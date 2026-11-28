#!/usr/bin/env python3
"""standardize_clonotypes.py — normalize MiXCR clones.tsv and TRUST4
trust_cdr3.out into the unified clonotype schema (AGENTS.md §5.3).

Unified columns:
  barcode_28, barcode_15, barcode_12, cancer_type, tool, chain,
  cdr3_nt, cdr3_aa, v_gene, j_gene, c_gene, read_count, clone_fraction, clone_count,
  count_unit

chain is inferred from the V-gene prefix (TRAV/TRBV/TRGV/TRDV/IGHV/IGKV/IGLV)
or from the tool's own locus/chain column.

clone_fraction is the fraction of that (sample, tool, **chain**) read pool —
chains are kept separate on purpose: BCR (IGH/IGK/IGL) and TCR (TRA/TRB/...)
transcripts come from different cell populations, so a cross-chain denominator
would conflate two receptor systems. MiXCR's native readFraction uses the same
denominator; we recompute it from read_count so BOTH tools carry the column
(TRUST4 has no native fraction). Verified 2026-11-16: per-chain
sum(clone_fraction) == 1.0 for every MiXCR chain.

clone_count is the legacy MiXCR 3.x alias of readCount (MiXCR 4.x renamed
readCount/cloneCount), carried for compatibility with pre-4.x consumers; it
duplicates read_count by design. TRUST4 rows use consensus_count for both.

count_unit is resolved per sample for MiXCR from the readinfo file written by
the mixcr rule (paired-end -> read_pair, single-end -> read); TRUST4 always
counts reads. The output is written to a temp file and renamed into place, so
a failed re-run leaves the previous table intact.

Stdlib only — runs in the project pixi env.
"""
import argparse
import csv
import os
import re
import sys

COLUMNS = [
    "barcode_28", "barcode_15", "barcode_12", "cancer_type", "tool", "chain",
    "cdr3_nt", "cdr3_aa", "v_gene", "j_gene", "c_gene",
    "read_count", "clone_fraction", "clone_count", "count_unit",
]

# COUNTS ARE NOT COMPARABLE ACROSS TOOLS (review 2026-11-18, H10).
# MiXCR's readCount counts read PAIRS, TRUST4's consensus_count counts READS.
# Verified on TCGA-13-0886:
#   * MiXCR: sum(readCount) = 31,212, exactly the assemble report's
#     readsInClones; and align totalReadsProcessed = 138,433,853 = primary
#     BAM records / 2 — so MiXCR's counting unit really is the pair.
#   * TRUST4: sum(consensus_count) = 49,871, exactly the `#count` column of
#     trust4_report.tsv, and trust-airr.pl sources it from trust_cdr3.out's
#     abundance column. What we verified is that it is TRUST4's own native
#     abundance; that the upstream tool calls it a read count comes from its
#     documentation, which the conda package does not ship.
# The raw ratio is 1.598 overall and locus-dependent (1.37x IGL … 1.80x IGH;
# per-clone median 1.0, IQR [1,2]). Converting MiXCR pairs to reads (x2, valid
# only for paired-end) brings the two tools to 0.70–1.05 per chain — i.e. most
# of the apparent discrepancy IS the unit difference, not a disagreement.
# Never sum or compare read_count across tools; count_unit says so in-band.
# A fixed x2 conversion is NOT applied here because MiXCR's unit is one read
# for single-end libraries — the multiplier depends on read_layout, which the
# mixcr rule writes to {barcode}.mixcr.readinfo.tsv (column 3). The unit is
# therefore resolved per sample; `read_pair` is the safe default when the
# readinfo file is missing (all 400 instances completed by 2026-11-27 were
# paired-end).
COUNT_UNIT = {"mixcr": "read_pair", "trust4": "read"}


def count_unit_for_layout(layout: str) -> str:
    """MiXCR counts pairs on paired-end libraries and reads on single-end."""
    return "read" if (layout or "").strip().lower() == "single-end" else "read_pair"


def read_layout(mixcr_path: str, barcode: str) -> str:
    """Read the layout column written by the mixcr rule (readinfo.tsv).

    The QC table reads the same file; keeping the two in sync means the
    unified `count_unit` always matches the QC `read_layout` column.
    """
    path = os.path.join(os.path.dirname(os.path.abspath(mixcr_path)),
                        f"{barcode}.mixcr.readinfo.tsv")
    try:
        with open(path, newline="") as f:
            row = next(csv.reader(f, delimiter="\t"))
        return row[2].strip().lower()
    except (OSError, StopIteration, IndexError):
        return ""


V_PREFIX_TO_CHAIN = {
    "TRAV": "TRA", "TRBV": "TRB", "TRGV": "TRG", "TRDV": "TRD",
    "IGHV": "IGH", "IGKV": "IGK", "IGLV": "IGL",
}

# J-gene prefixes let infer_chain fall back when V is empty (e.g., TRUST4
# rows with partial annotation).
J_PREFIX_TO_CHAIN = {
    "TRAJ": "TRA", "TRBJ": "TRB", "TRGJ": "TRG", "TRDJ": "TRD",
    "IGHJ": "IGH", "IGKJ": "IGK", "IGLJ": "IGL",
}


def infer_chain(v_gene: str, j_gene: str = "") -> str:
    for gene in (v_gene, j_gene):
        if not gene:
            continue
        for table in (V_PREFIX_TO_CHAIN, J_PREFIX_TO_CHAIN):
            for prefix, chain in table.items():
                if gene.startswith(prefix):
                    return chain
    return ""


def strip_allele(gene: str) -> str:
    """Normalize gene names to gene level (drop *XX allele suffix).

    MiXCR writes placeholder alleles (*00) while TRUST4 reports IMGT alleles
    (*01, *02...). Gene-level comparison is the correct cross-method
    granularity; allele strings would fragment usage tables (verified
    2026-11-15: it produced a spurious uniform Spearman ~ -0.85).
    """
    return gene.split("*")[0].strip() if gene else gene


def first_hit(field: str) -> str:
    """MiXCR *HitsWithScore like 'TRBV28(100,25/262);TRBV20-1(55,...)' or a bare gene."""
    if not field:
        return ""
    first = field.split(";")[0].split(",")[0]
    return first.split("(")[0].strip()


def parse_mixcr(path, barcode, cancer, count_unit=COUNT_UNIT["mixcr"]):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            v = strip_allele(first_hit(r.get("allVHitsWithScore", "")))
            j = strip_allele(first_hit(r.get("allJHitsWithScore", "")))
            c = strip_allele(first_hit(r.get("allCHitsWithScore", "")))
            aa = (r.get("aaSeqCDR3") or "").strip()
            nt = (r.get("nSeqCDR3") or "").strip()
            if not aa and not nt:
                continue
            rows.append({
                "barcode_28": barcode,
                "barcode_15": "-".join(barcode.split("-")[:4]),
                "barcode_12": "-".join(barcode.split("-")[:3]),
                "cancer_type": cancer,
                "tool": "mixcr",
                # merged per-chain exports carry a chain column; fall back to
                # gene-prefix inference when absent
                "chain": (r.get("chain") or "").strip() or infer_chain(v, j),
                "cdr3_nt": nt,
                "cdr3_aa": aa,
                "v_gene": v,
                "j_gene": j,
                "c_gene": c,
                # MiXCR 4.x export uses readCount/readFraction (3.x used
                # cloneCount/cloneFraction) — accept both; counts normalized
                # to int ('4472.0' -> 4472) for downstream numeric handling.
                # clone_fraction is assigned uniformly for both tools by
                # assign_chain_fractions() (see module docstring).
                "read_count": safe_int(r.get("readCount") or r.get("cloneCount") or ""),
                "clone_fraction": "",
                "clone_count": safe_int(r.get("readCount") or r.get("cloneCount") or ""),
                "count_unit": count_unit,
            })
    return rows


def parse_trust4_airr(path, barcode, cancer):
    """Parse TRUST4's AIRR output (*_airr.tsv) — the standardized schema.

    TRUST4's trust_cdr3.out is NOT stable to parse positionally: v1.1.5
    emits 13 columns with a leading contig id (e.g. `assemble0`) that shifts
    every assumed position (verified 2026-11-14). The AIRR file has named
    columns: locus=chain, junction=CDR3 nt (IMGT junction incl. anchors),
    junction_aa=CDR3 aa, consensus_count=reads.
    """
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            aa = (r.get("junction_aa") or "").strip()
            nt = (r.get("junction") or "").strip()
            if (not aa or aa == "*") and (not nt or nt == "*"):
                continue
            v = strip_allele(first_hit(r.get("v_call", "")))
            j = strip_allele(first_hit(r.get("j_call", "")))
            c = strip_allele(first_hit(r.get("c_call", "")))
            chain = (r.get("locus") or "").strip() or infer_chain(v, j)
            count = (r.get("consensus_count") or "").strip()
            rows.append({
                "barcode_28": barcode,
                "barcode_15": "-".join(barcode.split("-")[:4]),
                "barcode_12": "-".join(barcode.split("-")[:3]),
                "cancer_type": cancer,
                "tool": "trust4",
                "chain": chain,
                "cdr3_nt": nt,
                "cdr3_aa": aa,
                "v_gene": v,
                "j_gene": j,
                "c_gene": c,
                "read_count": safe_int(count),
                "clone_fraction": "",
                "clone_count": safe_int(count),
                "count_unit": COUNT_UNIT["trust4"],
            })
    return rows


def safe_int(x):
    """Normalize numeric-looking counts to int ('4472.0' -> 4472); passthrough otherwise.

    OverflowError covers 'inf' (float('inf') raises on int()); NaN raises
    ValueError. Both are passed through unchanged rather than crashing the
    whole sample.
    """
    try:
        f = float(x)
        return int(f) if f == int(f) else f
    except (TypeError, ValueError, OverflowError):
        return x


def assign_chain_fractions(rows):
    """Return NEW rows with clone_fraction = read_count / (tool, chain) total.

    One denominator per (sample, tool, chain) — see the module docstring for
    why chains are not pooled. Rows whose read_count is not numeric keep an
    empty fraction rather than fabricating a number.
    """
    totals = {}
    for r in rows:
        rc = r["read_count"]
        if isinstance(rc, (int, float)) and rc > 0:
            key = (r["tool"], r["chain"])
            totals[key] = totals.get(key, 0) + rc
    out = []
    for r in rows:
        rc = r["read_count"]
        total = totals.get((r["tool"], r["chain"]), 0)
        frac = f"{rc / total:.6g}" if isinstance(rc, (int, float)) and total > 0 else ""
        out.append({**r, "clone_fraction": frac})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mixcr", required=True)
    ap.add_argument("--trust4-airr", required=True, help="TRUST4 *_airr.tsv (AIRR schema)")
    ap.add_argument("--barcode", required=True, help="full 28-char TCGA barcode")
    ap.add_argument("--cancer", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    if not re.match(r"^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[A-Z0-9]{3}-[A-Z0-9]{3}-[A-Z0-9]{4}-[A-Z0-9]{2}$", args.barcode):
        print(f"WARNING: barcode {args.barcode} does not match TCGA 28-char format", file=sys.stderr)

    layout = read_layout(args.mixcr, args.barcode)
    mixcr_unit = count_unit_for_layout(layout)
    rows = parse_mixcr(args.mixcr, args.barcode, args.cancer, mixcr_unit)
    rows += parse_trust4_airr(args.trust4_airr, args.barcode, args.cancer)
    rows = assign_chain_fractions(rows)

    # Write then rename: a crash mid-write must never leave a truncated
    # unified.tsv under the previous run's still-present sentinel. The cluster
    # backend does not roll back failed outputs (unlike the local executor).
    tmp = f"{args.output}.tmp.{os.getpid()}"
    try:
        with open(tmp, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        os.replace(tmp, args.output)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)

    n_mixcr = sum(1 for r in rows if r["tool"] == "mixcr")
    n_trust4 = sum(1 for r in rows if r["tool"] == "trust4")
    print(f"wrote {args.output}: {len(rows)} rows (mixcr={n_mixcr}, trust4={n_trust4}; "
          f"count_unit: mixcr={mixcr_unit}, trust4={COUNT_UNIT['trust4']}; "
          f"layout={layout or 'unknown'})")


if __name__ == "__main__":
    main()
