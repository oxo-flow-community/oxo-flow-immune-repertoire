#!/usr/bin/env python3
"""compare_tools.py — cross-method validation metrics for one sample:
MiXCR vs TRUST4 on the unified clonotype table (or the two raw outputs).

Metrics (PLAN.md §4.4):
  1. top-100 CDR3aa Jaccard similarity                     (pass > 0.5)
  2. V-gene usage Spearman rank correlation                 (reported)
  3. per-clone clone_fraction Pearson r (log10, matched by CDR3aa, per chain)
     (pass > 0.7). Read counts are NOT comparable across tools: MiXCR counts
     read pairs, TRUST4 counts reads, and the ratio is locus-dependent. The
     per-chain clone_fraction is dimensionless and is the valid quantity.

Stdlib only (no numpy/scipy): ranks and Pearson implemented directly.
"""
import argparse
import csv
import math
import re
import statistics
import sys
from collections import Counter

# Comparable subset: complete, in-frame CDR3 amino-acid sequences
# (IMGT-anchored: leading C, trailing F/W, no stop/gap). Both tools emit
# partial or out-of-frame junctions for some contigs; comparing only the
# complete subset is the standard practice for cross-method agreement.
COMPLETE_CDR3 = re.compile(r"^C[A-Z]*[FW]$")


def is_complete_cdr3(aa: str) -> bool:
    return bool(aa) and "*" not in aa and bool(COMPLETE_CDR3.match(aa))


def read_unified(path):
    by_tool = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            by_tool.setdefault(r["tool"], []).append(r)
    return by_tool


def top_cdr3aa(rows, n=100):
    """Top-n UNIQUE complete CDR3aa by best read count.

    Independent verification (2026-11-15) showed first-100-rows-then-set can
    yield <100 uniques when the same junction appears on multiple contigs
    (OV trust4: 8,957 duplicate-aa rows). Dedupe-first is the standard
    clonotype-set definition.
    """
    best = {}
    for r in rows:
        aa = r.get("cdr3_aa", "")
        if not is_complete_cdr3(aa):
            continue
        c = float(r.get("read_count") or 0)
        if aa not in best or c > best[aa]:
            best[aa] = c
    ranked = sorted(best.items(), key=lambda kv: kv[1], reverse=True)
    return {aa for aa, _ in ranked[:n]}


def jaccard(a, b):
    if not a and not b:
        return 0.0
    inter = len(a & b)
    return inter / (len(a | b) or 1)


def v_gene_usage(rows):
    usage = Counter()
    for r in rows:
        v = r.get("v_gene", "")
        if v:
            usage[v] += float(r.get("read_count") or 0)
    return usage


def ranks(values):
    """Average ranks for ties, 1-based."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    rk = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            rk[order[k]] = avg
        i = j + 1
    return rk


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return 0.0
    mx = sum(xs) / n
    my = sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return 0.0
    return cov / math.sqrt(vx * vy)


def spearman(a_counter, b_counter):
    genes = sorted(set(a_counter) | set(b_counter))
    if len(genes) < 2:
        return 0.0
    xs = [a_counter.get(g, 0.0) for g in genes]
    ys = [b_counter.get(g, 0.0) for g in genes]
    return pearson(ranks(xs), ranks(ys))


def fraction_map(rows):
    """Sum clone_fraction per CDR3aa (per-chain normalized, dimensionless)."""
    m = {}
    for r in rows:
        aa = r.get("cdr3_aa", "")
        if not is_complete_cdr3(aa):
            continue
        try:
            frac = float(r.get("clone_fraction") or "")
        except (TypeError, ValueError):
            continue
        m[aa] = m.get(aa, 0.0) + frac
    return m


def by_chain(rows):
    d = {}
    for r in rows:
        d.setdefault(r.get("chain", ""), []).append(r)
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--unified", required=True, help="*.unified.tsv from standardize_clonotypes.py")
    ap.add_argument("--topn", type=int, default=100)
    args = ap.parse_args()

    by_tool = read_unified(args.unified)
    mx = by_tool.get("mixcr", [])
    tr = by_tool.get("trust4", [])
    if not mx or not tr:
        print("ERROR: unified file must contain both mixcr and trust4 rows", file=sys.stderr)
        sys.exit(1)

    jac = jaccard(top_cdr3aa(mx, args.topn), top_cdr3aa(tr, args.topn))
    v_rho = spearman(v_gene_usage(mx), v_gene_usage(tr))

    common = sorted(set(fraction_map(mx)) & set(fraction_map(tr)))
    per_chain = []
    mx_chains, tr_chains = by_chain(mx), by_chain(tr)
    for chain in sorted(set(mx_chains) & set(tr_chains)):
        fm = fraction_map(mx_chains[chain])
        ft = fraction_map(tr_chains[chain])
        shared = sorted(set(fm) & set(ft))
        if len(shared) >= 2:
            per_chain.append((chain, len(shared), pearson(
                [math.log10(fm[c] + 1e-12) for c in shared],
                [math.log10(ft[c] + 1e-12) for c in shared])))
    frac_r = statistics.median([r for _, _, r in per_chain]) if per_chain else 0.0

    verdict_j = "PASS" if jac > 0.5 else "FAIL"
    verdict_r = "PASS" if frac_r > 0.7 else "FAIL"
    print(f"sample: {args.unified}")
    print(f"clonotypes: mixcr={len(mx)} trust4={len(tr)} common_cdr3aa={len(common)}")
    print(f"top{args.topn}_cdr3aa_jaccard: {jac:.4f}  ({verdict_j}, threshold 0.5)")
    print(f"v_gene_usage_spearman:       {v_rho:.4f}")
    print(f"clone_fraction_pearson_median: {frac_r:.4f}  "
          f"({verdict_r}, threshold 0.7, chains={len(per_chain)})")
    for chain, n, r in per_chain:
        print(f"  {chain}: n={n} r={r:.4f}")


if __name__ == "__main__":
    main()
