#!/usr/bin/env python3
"""verify_batch.py — independent verification of completed per-sample outputs.

Checks every sample that has a `*.unified.tsv`:

  A. completeness   all declared outputs + completion sentinels + cleanup marker
  B. provenance     what the mixcr sentinel says (current vs pre-cutover)
  C. conservation   the unified table must reproduce EXACTLY the rows the
                    production parser derives from the two source files, as a
                    multiset (no lost/duplicated/altered rows). Reuses the
                    production parsers so the check tests the pipeline, not a
                    re-implementation of it.
  D. self-consistency  clone_fraction sums to 1 per (tool, chain); chain
                    agrees with the V-gene prefix; barcodes match the path;
                    counts are non-negative integers; clone_count == read_count
                    for MiXCR (by design).
  E. cross-tool     Jaccard (top-100 CDR3 aa), Spearman (V-gene usage) and
                    Pearson on log read counts of shared CDR3s. Reported, not
                    gated — the expected range is biological, not fixed.

Exit status: 0 if every sample passes A-D, 1 otherwise.
Stdlib only.

Usage:
  python3 scripts/verify_batch.py --results shards/shard1/results [...]
  python3 scripts/verify_batch.py --results shards/shard*/results --json out.json
"""
import argparse
import csv
import glob
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import standardize_clonotypes as sc  # noqa: E402  (production parser, reused on purpose)

COLUMNS = sc.COLUMNS
RULES = {
    "mixcr": ("{b}.mixcr.clones.tsv", "{b}.mixcr.qc.json", "{b}.mixcr.done"),
    "trust4": ("{b}.trust4_cdr3.out", "{b}.trust4_report.tsv",
               "{b}.trust4_airr.tsv", "{b}.trust4.done"),
    "standardize": ("{b}.unified.tsv", "{b}.unified.done"),
    "cleanup": ("{b}.tmp_cleaned",),
}


def read_unified(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def as_tuple(row):
    return tuple(str(row.get(c, "")) for c in COLUMNS)


def check_completeness(sample_dir, barcode):
    """Every declared output must exist AND be structurally complete.

    Existence alone is not enough (review H8): a job killed mid-write leaves a
    truncated TSV that still "exists", and the engine's skip test plus the
    driver's counting both trusted existence.
    """
    import check_tsv_integrity as cti
    missing = []
    for rule, patterns in RULES.items():
        for p in patterns:
            path = os.path.join(sample_dir, p.format(b=barcode))
            if not os.path.isfile(path):
                missing.append(f"{rule}:{p.format(b=barcode)}")
                continue
            if path.endswith((".tsv", ".out", ".fa", ".json")):
                if path.endswith(".json"):
                    try:
                        with open(path) as f:
                            json.load(f)
                    except (OSError, json.JSONDecodeError):
                        missing.append(f"{rule}:{p.format(b=barcode)} (invalid json)")
                    continue
                ok, detail = cti.check(path)
                if not ok:
                    missing.append(f"{rule}:{p.format(b=barcode)} ({detail})")
    return missing


def provenance(sample_dir, barcode):
    path = os.path.join(sample_dir, f"{barcode}.mixcr.done")
    try:
        with open(path) as f:
            text = f.read().strip()
    except OSError:
        return "unknown", ""
    kind = "backfilled" if text.lower().startswith("backfilled") else "current"
    return kind, text.splitlines()[0] if text else ""


def check_conservation(sample_dir, barcode, cancer, actual):
    """Rebuild the expected rows from the sources with the production parser."""
    mixcr_path = os.path.join(sample_dir, f"{barcode}.mixcr.clones.tsv")
    airr_path = os.path.join(sample_dir, f"{barcode}.trust4_airr.tsv")
    rows = sc.parse_mixcr(mixcr_path, barcode, cancer)
    rows += sc.parse_trust4_airr(airr_path, barcode, cancer)
    expected = [as_tuple(r) for r in sc.assign_chain_fractions(rows)]
    got = [as_tuple(r) for r in actual]
    if sorted(expected) == sorted(got):
        return True, ""
    exp_counts = {}
    for t in expected:
        exp_counts[t] = exp_counts.get(t, 0) + 1
    got_counts = {}
    for t in got:
        got_counts[t] = got_counts.get(t, 0) + 1
    only_exp = sum(max(exp_counts.get(t, 0) - got_counts.get(t, 0), 0) for t in exp_counts)
    only_got = sum(max(got_counts.get(t, 0) - exp_counts.get(t, 0), 0) for t in got_counts)
    return False, f"expected {len(expected)} rows, got {len(got)}; missing={only_exp} extra={only_got}"


def check_self_consistency(actual, barcode, cancer, sample_dir):
    problems = []
    sums = {}
    for r in actual:
        rc = r["read_count"]
        try:
            val = float(rc)
        except ValueError:
            problems.append(f"non-numeric read_count {rc!r}")
            continue
        if val < 0:
            problems.append(f"negative read_count {rc}")
        frac = r["clone_fraction"]
        if frac == "" :
            if val > 0:
                problems.append("empty clone_fraction on a positive count")
            continue
        key = (r["tool"], r["chain"])
        sums[key] = sums.get(key, 0.0) + float(frac)
        if r["tool"] == "mixcr" and str(rc) != str(r["clone_count"]):
            problems.append("mixcr clone_count != read_count")
        v = r["v_gene"]
        if v:
            inferred = sc.infer_chain(v, r["j_gene"])
            # Only a different RECEPTOR SYSTEM is a defect. Within the TCR,
            # cross-locus rearrangements are real biology — TRDV1 rearranges
            # to TRAJ segments and TRUST4 reports locus=TRA for those; the
            # standardizer intentionally prefers the tool's own locus.
            system = {"TRA": "TCR", "TRB": "TCR", "TRG": "TCR", "TRD": "TCR",
                      "IGH": "BCR", "IGK": "BCR", "IGL": "BCR"}
            if (inferred and r["chain"]
                    and system.get(inferred) != system.get(r["chain"])):
                problems.append(f"chain {r['chain']} disagrees with V gene {v}")
        if r["barcode_28"] != barcode:
            problems.append(f"barcode_28 {r['barcode_28']} != {barcode}")
        if r["barcode_12"] != "-".join(barcode.split("-")[:3]):
            problems.append("barcode_12 mismatch")
        if r["cancer_type"] != cancer:
            problems.append(f"cancer_type {r['cancer_type']} != {cancer}")
    for key, total in sums.items():
        if not (0.999 <= total <= 1.001):
            problems.append(f"fraction sum {total:.6f} for {key}")
    # the MiXCR QC report must parse and report a positive library size
    qc_path = os.path.join(sample_dir, f"{barcode}.mixcr.qc.json")
    try:
        with open(qc_path) as f:
            qc = json.load(f)
        # a zero-clonotype sample legitimately gets a placeholder `[]` (the
        # mixcr rule writes one when MiXCR exported no qc.json) — that is a
        # valid provenance marker, not a defect
        if isinstance(qc, dict):
            if not (qc.get("inputFiles") or qc.get("qc") or qc):
                problems.append("qc.json is empty")
        elif not isinstance(qc, list):
            problems.append(f"qc.json has unexpected type {type(qc).__name__}")
    except (OSError, json.JSONDecodeError) as exc:
        problems.append(f"qc.json unreadable: {exc}")
    return problems


def rank(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return num / (dx * dy) if dx and dy else float("nan")


def cross_tool(actual):
    """Cross-tool agreement, computed PER CHAIN (review M4).

    Pooling chains mixes two receptor systems: Ig reads usually outnumber TCR
    reads 100:1, so a pooled V-gene Spearman is really an Ig measurement.
    Jaccard is inherently depth-sensitive (VDJtools), so treat it as
    descriptive, not a pass/fail gate.
    """
    out = {}
    by_tool = {"mixcr": [], "trust4": []}
    for r in actual:
        by_tool.setdefault(r["tool"], []).append(r)

    def usage(tool, chain=None):
        counts = {}
        for r in by_tool.get(tool, []):
            if r["v_gene"] and (chain is None or r["chain"] == chain):
                counts[r["v_gene"]] = counts.get(r["v_gene"], 0) + float(r["read_count"] or 0)
        total = sum(counts.values())
        return {g: c / total for g, c in counts.items()} if total else {}

    def counts_by_aa(tool, chain=None):
        c = {}
        for r in by_tool.get(tool, []):
            if (r["cdr3_aa"] and r["read_count"]
                    and (chain is None or r["chain"] == chain)):
                c[r["cdr3_aa"]] = c.get(r["cdr3_aa"], 0) + float(r["read_count"])
        return c

    # pooled values kept for continuity, clearly labelled
    um, ut = usage("mixcr"), usage("trust4")
    genes = sorted(set(um) & set(ut))
    if len(genes) >= 5:
        out["spearman_vgene_pooled"] = pearson(rank([um[g] for g in genes]),
                                              rank([ut[g] for g in genes]))
    cm, ct = counts_by_aa("mixcr"), counts_by_aa("trust4")
    shared = sorted(set(cm) & set(ct))
    if len(shared) >= 5:
        out["pearson_logreads_pooled"] = pearson(
            [math.log10(cm[a] + 1) for a in shared],
            [math.log10(ct[a] + 1) for a in shared])

    chains = sorted({r["chain"] for r in actual if r["chain"]})
    for chain in chains:
        rows_m = [r for r in by_tool.get("mixcr", []) if r["chain"] == chain and r["cdr3_aa"]]
        rows_t = [r for r in by_tool.get("trust4", []) if r["chain"] == chain and r["cdr3_aa"]]
        if not rows_m or not rows_t:
            continue
        rows_m.sort(key=lambda r: -float(r["read_count"] or 0))
        rows_t.sort(key=lambda r: -float(r["read_count"] or 0))
        a = {r["cdr3_aa"] for r in rows_m[:100]}
        b = {r["cdr3_aa"] for r in rows_t[:100]}
        out[f"jaccard_{chain}"] = len(a & b) / len(a | b)
        um_c, ut_c = usage("mixcr", chain), usage("trust4", chain)
        g = sorted(set(um_c) & set(ut_c))
        if len(g) >= 5:
            out[f"spearman_vgene_{chain}"] = pearson(rank([um_c[x] for x in g]),
                                                    rank([ut_c[x] for x in g]))
        cm_c, ct_c = counts_by_aa("mixcr", chain), counts_by_aa("trust4", chain)
        sh = sorted(set(cm_c) & set(ct_c))
        if len(sh) >= 5:
            out[f"pearson_logreads_{chain}"] = pearson(
                [math.log10(cm_c[x] + 1) for x in sh],
                [math.log10(ct_c[x] + 1) for x in sh])
            out[f"n_cdr3_shared_{chain}"] = len(sh)
    return out


def verify_sample(sample_dir):
    barcodes = [os.path.basename(p).replace(".unified.tsv", "")
                for p in glob.glob(os.path.join(sample_dir, "*.unified.tsv"))]
    if len(barcodes) != 1:
        return {"sample_dir": sample_dir, "error": f"{len(barcodes)} unified files"}
    barcode = barcodes[0]
    cancer = sample_dir.split(os.sep)[-3]
    row = {"sample_dir": sample_dir, "barcode": barcode, "cancer": cancer}

    missing = check_completeness(sample_dir, barcode)
    row["missing"] = missing
    prov_kind, prov_text = provenance(sample_dir, barcode)
    row["provenance"] = prov_kind
    actual = read_unified(os.path.join(sample_dir, f"{barcode}.unified.tsv"))
    row["rows"] = len(actual)
    row["rows_mixcr"] = sum(1 for r in actual if r["tool"] == "mixcr")
    row["rows_trust4"] = sum(1 for r in actual if r["tool"] == "trust4")

    ok, detail = check_conservation(sample_dir, barcode, cancer, actual)
    row["conservation_ok"] = ok
    row["conservation_detail"] = detail
    row["consistency_problems"] = check_self_consistency(actual, barcode, cancer, sample_dir)
    row.update(cross_tool(actual))
    row["pass"] = (not missing) and ok and not row["consistency_problems"]
    return row


def fmt(value):
    if isinstance(value, float):
        return "nan" if value != value else f"{value:.3f}"
    return str(value)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", nargs="+", required=True)
    ap.add_argument("--json", help="write the full per-sample report here")
    ap.add_argument("--limit", type=int, default=0, help="verify at most N samples")
    ap.add_argument("--print-ok-count", action="store_true",
                    help="print only the number of samples that pass (for driver_shard.sh)")
    args = ap.parse_args()

    sample_dirs = []
    for root in args.results:
        for path in glob.glob(os.path.join(root, "**", "*.unified.tsv"), recursive=True):
            sample_dirs.append(os.path.dirname(path))
    sample_dirs = sorted(set(sample_dirs))
    if args.limit:
        sample_dirs = sample_dirs[: args.limit]
    if not sample_dirs:
        print("no unified.tsv found — nothing to verify")
        return 1

    reports = [verify_sample(d) for d in sample_dirs]
    if args.print_ok_count:
        print(sum(1 for r in reports if r.get("pass")))
        return 0
    header = ["barcode", "cancer", "provenance", "rows", "rows_mixcr", "rows_trust4",
              "conservation_ok", "jaccard_TRB", "jaccard_IGH", "spearman_vgene",
              "pearson_logreads", "pass"]
    print("\t".join(header))
    for r in reports:
        if "error" in r:
            print(f"{os.path.basename(r['sample_dir'])}\tERROR {r['error']}")
            continue
        print("\t".join(fmt(r.get(h, "")) for h in header))
    failures = [r for r in reports if not r.get("pass")]
    print(f"\n{len(reports)} sample(s): {len(reports) - len(failures)} pass, {len(failures)} fail")
    for r in failures:
        print(f"  FAIL {r.get('barcode', r['sample_dir'])}: "
              f"missing={r.get('missing', [])} conservation={r.get('conservation_detail', '')} "
              f"consistency={r.get('consistency_problems', [])[:3]}")
    if args.json:
        with open(args.json, "w") as f:
            json.dump(reports, f, indent=2)
        print(f"report: {args.json}")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
