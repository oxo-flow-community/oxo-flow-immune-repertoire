#!/usr/bin/env python3
"""sample_metrics.py — per-sample repertoire QC metrics table.

Produces manifests/../qc/sample_metrics.tsv, one row per sample, with the
covariates needed downstream (AGENTS.md §7): library size, immune-content
fraction, diversity/clonality, chain composition. Stdlib only.

Metrics rationale:
  total_reads          MiXCR totalReadsProcessed = sequenced pairs entering align
  aligned_reads/frac   MiXCR successfully aligned (repertoire extraction yield)
  cdr3_reads_*         sum of clonotype read counts per tool (repertoire captured)
  cdr3_frac            cdr3_reads / total_reads — immune-content covariate; the
                       primary filter candidate (low => no repertoire to compare)
  clonotypes           productive clonotype count per tool (richness)
  clonality            1 - Shannon H / ln(S) on clonotype frequencies (0 = polyclonal,
                       1 = oligoclonal; standard in immuno-oncology literature)
  gini                 Gini coefficient of clonotype read distribution (skew)
  top_clone_frac       largest clonotype fraction (dominance / QC per AGENTS.md §7.1)
  chain fractions      IGH/IGK/IGL/TRA/TRB/(TRD/TRG) share of cdr3 reads per tool
  mean_cdr3_nt         mean CDR3 nucleotide length per tool (QC + biology)
"""
import argparse
import csv
import json
import math
import os
import sys
from collections import defaultdict

CHAIN_ORDER = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"]


def safe_float(x, default=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def gini(counts):
    xs = sorted(counts)
    n = len(xs)
    if n == 0 or sum(xs) == 0:
        return 0.0
    cum = sum((i + 1) * x for i, x in enumerate(xs))
    return (2 * cum) / (n * sum(xs)) - (n + 1) / n


def clonality(counts):
    """1 - Shannon H / ln(S) computed on clone fractions."""
    total = sum(counts)
    if total <= 0:
        return 0.0
    ps = [c / total for c in counts if c > 0]
    s = len(ps)
    if s <= 1:
        return 1.0
    h = -sum(p * math.log(p) for p in ps)
    return 1.0 - h / math.log(s)


def load_unified(path):
    rows = defaultdict(list)
    with open(path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows[r.get("tool") or ""].append(r)
    return rows


def tool_block(rows):
    if not rows:
        return {k: "" for k in (
            "clonotypes", "clonotypes_raw", "cdr3_reads", "clonality", "gini",
            "top_clone_frac", "mean_cdr3_nt", *[f"frac_{c}" for c in CHAIN_ORDER])}
    counts = [max(1, int(safe_float(r.get("read_count"), 0))) for r in rows]
    chains = defaultdict(float)
    lens = []
    for r in rows:
        chains[r.get("chain", "")] += max(1, int(safe_float(r.get("read_count"), 0)))
        if r.get("cdr3_nt"):
            lens.append(len(r["cdr3_nt"]))
    total = sum(counts)
    # top clone by read share — per-chain clone_fraction columns are normalized
    # WITHIN chain (per-chain exports), so their max can be 1.0 spuriously.
    top = max(counts) / total if total else 0.0
    # clonotypes = productive (complete in-frame CDR3aa) rows only; raw row
    # count kept as clonotypes_raw. Independent verification 2026-11-15 showed
    # raw counts overstate richness (OV trust4: 41,761 raw vs 27,738 complete).
    n_complete = sum(1 for r in rows if r.get("cdr3_aa") and "*" not in r["cdr3_aa"]
                     and r["cdr3_aa"].startswith("C") and r["cdr3_aa"][-1:] in ("F", "W"))
    block = {
        "clonotypes": n_complete,
        "clonotypes_raw": len(rows),
        "cdr3_reads": total,
        "clonality": round(clonality(counts), 4),
        "gini": round(gini(counts), 4),
        "top_clone_frac": round(top, 4),
        "mean_cdr3_nt": round(sum(lens) / len(lens), 1) if lens else "",
    }
    for c in CHAIN_ORDER:
        block[f"frac_{c}"] = round(chains.get(c, 0.0) / total, 4) if total else 0.0
    return block


def mixcr_total_reads(align_report_json):
    try:
        with open(align_report_json) as f:
            rep = json.load(f)
        # v4.7: {"type":"alignerReport", ...,"totalReadsProcessed":N}
        if isinstance(rep, dict):
            for key in ("totalReadsProcessed", "totalSequencingReads"):
                if key in rep:
                    return int(rep[key])
        # older layout: list of blocks with key/value pairs
        if isinstance(rep, list):
            for block in rep:
                if isinstance(block, dict) and block.get("key") == "totalReadsProcessed":
                    return int(block.get("value", 0))
    except (OSError, ValueError, TypeError):
        pass
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-root", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    meta = []
    with open(args.manifest, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            # oxo-flow metadata files key the sample column as 'sample';
            # manifests/all_samples.tsv uses 'sample_id' — accept both.
            if "sample" in r and "sample_id" not in r:
                r["sample_id"] = r["sample"]
            meta.append(r)
    by_sample = {r["sample_id"]: r for r in meta}

    # sample dirs are results/<cancer>/<barcode_12>/<sample_id>/ — exactly 3
    # levels below the results root (absolute-depth comparison; relative
    # separator counting breaks when results_root is relative).
    root_abs = os.path.abspath(args.results_root)
    base_depth = len(root_abs.split(os.sep))
    sample_dirs = []
    for root, _dirs, files in os.walk(root_abs):
        if len(os.path.abspath(root).split(os.sep)) - base_depth == 3 and any(
                x.endswith(".unified.tsv") for x in files):
            sample_dirs.append(root)
    dirs_by_uuid = {os.path.basename(d): d for d in sample_dirs}

    header = None
    out_rows = []
    for uuid, r in by_sample.items():
        d = dirs_by_uuid.get(uuid)
        if not d:
            continue
        bc = r["barcode_28"]
        unified = os.path.join(d, f"{bc}.unified.tsv")
        align_json = os.path.join(d, f"{bc}.mixcr.align.report.json")
        if not os.path.exists(unified):
            continue
        tools = load_unified(unified)
        mx = tool_block(tools.get("mixcr", []))
        tr = tool_block(tools.get("trust4", []))
        total_reads = mixcr_total_reads(align_json)
        row = {
            "sample_id": uuid,
            "barcode_28": bc,
            "cancer_type": r["cancer_type"],
            "sample_type": r["barcode_28"].split("-")[3][:2],
            "has_chimeric": r.get("has_chimeric", ""),
            "bam_size_gb": round(safe_float(r.get("file_size_bytes")) / 1e9, 2),
            "total_reads": total_reads,
            "mixcr_aligned_reads": "",
            "mixcr_aligned_frac": "",
        }
        # aligned reads from the align report txt (stable text layout)
        txt = align_json.replace(".json", ".txt")
        try:
            with open(txt) as f:
                for line in f:
                    if line.startswith("Total sequencing reads:"):
                        row["total_reads"] = int(line.split(":")[1])
                    elif line.startswith("Successfully aligned reads:"):
                        parts = line.split(":")[1].split("(")
                        row["mixcr_aligned_reads"] = int(parts[0])
                        # report prints a PERCENT ("0.23%") — store as fraction
                        row["mixcr_aligned_frac"] = round(float(parts[1].rstrip("%)\n")) / 100, 6)
        except (OSError, ValueError, IndexError):
            pass
        for k, v in mx.items():
            row[f"mixcr_{k}"] = v
        row["mixcr_cdr3_frac"] = (
            round(safe_float(mx.get("cdr3_reads")) / row["total_reads"], 6)
            if row["total_reads"] and mx.get("cdr3_reads") != "" else "")
        for k, v in tr.items():
            row[f"trust4_{k}"] = v
        out_rows.append(row)

    out_rows.sort(key=lambda r: (r["cancer_type"], r["barcode_28"]))
    header = list(out_rows[0].keys()) if out_rows else ["sample_id"]
    with open(args.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)
    print(f"wrote {args.output}: {len(out_rows)} samples x {len(header)} metrics")


if __name__ == "__main__":
    main()
