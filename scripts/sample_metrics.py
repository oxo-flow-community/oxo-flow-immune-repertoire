#!/usr/bin/env python3
"""sample_metrics.py — per-sample repertoire QC metrics table.

One row per completed sample. Columns are documented in
docs/metrics_reference.md; the literature behind each definition is in
docs/scientific_rationale.md (section 5).

What this version changed after the 2026-11-17 scientific review:

* **productive-only base** (M3): diversity/richness/top-clone are computed on
  productive, complete CDR3aa rows only (starts with C, ends with F/W, no
  stop codon); the unfiltered counts are kept as `*_raw`.
* **TCR and BCR groups separately** (H7): pooling every chain lets the Ig
  fraction (usually >95 % of reads) dominate what is reported as "clonality",
  which then cannot support a T-cell statement. Thorsson 2018 computed the two
  lineages separately [SRC:thorsson2018].
* **rarefied diversity** (H6): richness/clonality/Gini all depend on depth,
  and cdr3_reads spans 4 orders of magnitude in this cohort. We therefore also
  report metrics after subsampling to a common depth — the requirement stated
  in the VDJtools paper/tutorial [SRC:shugay2015]. The subsample is seeded so
  the table is reproducible.
* **read length / layout and count unit** (H2/H10): TCGA mixes 76 bp
  single-end and 48 bp paired-end libraries, and extractability depends on
  both [SRC:bolotin2017]; MiXCR counts pairs while TRUST4 counts reads, so
  `mixcr_count_unit`/`trust4_count_unit` travel with the numbers.
* **chain-specific productive reads** (2026-11-27 review): `*_IGH_cdr3_reads`
  etc. expose the per-chain read floors that the authoritative "has
  repertoire" thresholds (§2a) require; they are computed on productive rows,
  unlike the read-weighted `frac_*` columns.
* **trust4_cdr3_frac** uses a read denominator (MiXCR total pairs x2 for PE)
  rather than dividing TRUST4 reads by MiXCR pairs, which would inflate it.
* **generation from the sentinel** (H3): the sentinel written by the current
  rule carries its parameters ("mixcr-v2 genomic-only threads=4 xmx=40g");
  older sentinels are labelled `legacy` instead of being called current.
* **has_chimeric / chimeric_analyzed** (H1): `has_chimeric` only says the
  chimeric BAM exists. `chimeric_analyzed` is the processing covariate: yes
  only for pre-cutover samples with has_chimeric=1 (the C1 fix dropped the
  chimeric merge; current samples all read the genomic BAM only).

Stdlib only.
"""
import argparse
import bisect
import csv
import json
import math
import os
import random
import sys
from collections import defaultdict

CHAIN_ORDER = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"]
TCR_CHAINS = {"TRA", "TRB", "TRG", "TRD"}
BCR_CHAINS = {"IGH", "IGK", "IGL"}
RAREFY_SEED = 17


def safe_float(x, default=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def count_unit_for_layout(layout):
    """Unified-table unit for a MiXCR sample: pairs on PE, reads on SE.

    Must match standardize_clonotypes.py's count_unit so the QC table and the
    unified table never disagree.
    """
    return "read_pair" if layout == "paired-end" else "read" if layout == "single-end" else ""


def cdr3_fraction(reads, total_reads, layout):
    """TRUST4 CDR3 reads as a fraction of sequencing reads.

    `total_reads` is MiXCR's totalReadsProcessed, which counts read PAIRS;
    TRUST4's consensus_count counts reads. Convert the denominator by the
    library layout (x2 for paired-end, x1 for single-end) rather than dividing
    reads by pairs (which would inflate the fraction ~2x). Empty when the
    layout is unknown — a wrong denominator is worse than a missing value.
    """
    if reads == "" or not total_reads:
        return ""
    if layout == "paired-end":
        denom = total_reads * 2
    elif layout == "single-end":
        denom = total_reads
    else:
        return ""
    return round(safe_float(reads) / denom, 10)


def bam_size_gb(row):
    """file_size_bytes from the manifest, else stat the BAM itself.

    The shard metadata.tsv intentionally omits file_size_bytes; it does carry
    the BAM path. A stat over the local/NFS mount is one syscall per sample.
    """
    if row.get("file_size_bytes"):
        return round(safe_float(row["file_size_bytes"]) / 1e9, 2)
    path = row.get("bam_genomic") or ""
    if path:
        try:
            return round(os.stat(path).st_size / 1e9, 2)
        except OSError:
            return ""
    return ""


def gini(counts):
    xs = sorted(counts)
    n = len(xs)
    if n == 0 or sum(xs) == 0:
        return 0.0
    cum = sum((i + 1) * x for i, x in enumerate(xs))
    return (2 * cum) / (n * sum(xs)) - (n + 1) / n


def clonality(counts):
    """1 - Shannon H / ln(S) on clone frequencies [SRC:thorsson2018]."""
    total = sum(counts)
    if total <= 0:
        return 0.0
    ps = [c / total for c in counts if c > 0]
    s = len(ps)
    if s <= 1:
        return 1.0
    h = -sum(p * math.log(p) for p in ps)
    return 1.0 - h / math.log(s)


def rarefy(counts, depth, seed=RAREFY_SEED):
    """Subsample clone counts without replacement down to `depth` reads.

    Exact hypergeometric draw implemented as "choose `depth` read positions
    out of `total`, then bucket them by clone". Deterministic for a seed, so
    the table is reproducible; a no-op when the sample is already shallower
    than `depth` (the flag `rarefied` says which happened).
    """
    total = sum(counts)
    if total <= depth:
        return list(counts), False
    rng = random.Random(seed)
    picked = rng.sample(range(total), depth)
    cum = []
    run = 0
    for c in counts:
        run += c
        cum.append(run)
    out = [0] * len(counts)
    for pos in picked:
        out[bisect.bisect_right(cum, pos)] += 1
    return out, True


def is_productive(row):
    """Complete, in-frame CDR3aa: starts C, ends F/W, no stop/gap, frame-consistent nt.

    `_` is MiXCR's frameshift marker (e.g. `CALPAL_SKLTF`): such a clone is not
    productive even though it carries no stop codon. The nucleotide length must
    be codon-complete whenever `cdr3_nt` is present (independent review,
    2026-11-23; regression test test/sample_metrics_selftest.py).
    """
    aa = (row.get("cdr3_aa") or "").strip()
    nt = (row.get("cdr3_nt") or "").strip()
    if not aa or "*" in aa or "_" in aa:
        return False
    if not aa.startswith("C") or aa[-1:] not in ("F", "W"):
        return False
    return not nt or len(nt) % 3 == 0


def load_unified(path):
    rows = defaultdict(list)
    with open(path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows[r.get("tool") or ""].append(r)
    return rows


def diversity_block(counts):
    if not counts:
        return {"clonotypes": "", "cdr3_reads": "", "clonality": "", "gini": "",
                "top_clone_frac": "", "mean_cdr3_nt": ""}
    total = sum(counts)
    return {
        "clonotypes": len([c for c in counts if c > 0]),
        "cdr3_reads": total,
        "clonality": max(0.0, round(clonality(counts), 4)),  # avoid -0.0
        "gini": max(0.0, round(gini(counts), 4)),
        "top_clone_frac": round(max(counts) / total, 4) if total else 0.0,
        "mean_cdr3_nt": "",
    }


def tool_block(rows, depth):
    empty = {k: "" for k in (
        "clonotypes", "clonotypes_raw", "cdr3_reads", "clonality", "gini",
        "top_clone_frac", "mean_cdr3_nt", "clonality_rare", "gini_rare",
        "rarefied", "tcr_clonotypes", "tcr_cdr3_reads", "tcr_clonality",
        "tcr_gini", "tcr_top_clone_frac", "bcr_clonotypes", "bcr_cdr3_reads",
        "bcr_clonality", "bcr_gini", "bcr_top_clone_frac",
        *[f"frac_{c}" for c in CHAIN_ORDER],
        *[f"{c}_cdr3_reads" for c in CHAIN_ORDER])}
    if not rows:
        return empty

    productive = [r for r in rows if is_productive(r)]
    def counts_of(subset):
        return [max(1, int(safe_float(r.get("read_count"), 0))) for r in subset]

    all_counts = counts_of(rows)
    prod_counts = counts_of(productive)
    block = diversity_block(prod_counts)

    # read-weighted CDR3 length on productive rows (review M7: an unweighted
    # mean is dragged down by singletons)
    lens = [(len(r["cdr3_nt"]), max(1, int(safe_float(r.get("read_count"), 0))))
            for r in productive if r.get("cdr3_nt")]
    block["mean_cdr3_nt"] = (round(sum(l * w for l, w in lens) / sum(w for _, w in lens), 1)
                             if lens else "")
    block["clonotypes_raw"] = len(rows)

    rare, did = rarefy(prod_counts, depth) if depth else (prod_counts, False)
    div = diversity_block(rare)
    block["clonality_rare"] = div["clonality"]
    block["gini_rare"] = div["gini"]
    block["rarefied"] = "yes" if did else "no"

    # lineage blocks (H7)
    for label, chains in (("tcr", TCR_CHAINS), ("bcr", BCR_CHAINS)):
        sub = [r for r in productive if r.get("chain") in chains]
        d = diversity_block(counts_of(sub))
        for k in ("clonotypes", "cdr3_reads", "clonality", "gini", "top_clone_frac"):
            block[f"{label}_{k}"] = d[k]

    per_chain = defaultdict(int)
    for r in rows:
        per_chain[r.get("chain", "")] += max(1, int(safe_float(r.get("read_count"), 0)))
    total = sum(all_counts)
    for c in CHAIN_ORDER:
        block[f"frac_{c}"] = round(per_chain.get(c, 0) / total, 4) if total else 0.0

    # Chain-specific productive reads: the authoritative "has repertoire"
    # thresholds (docs/metrics_reference.md §2a) are per chain, and the
    # read-weighted *_frac columns above include non-productive rows.
    per_chain_prod = defaultdict(int)
    for r in productive:
        per_chain_prod[r.get("chain", "")] += max(1, int(safe_float(r.get("read_count"), 0)))
    for c in CHAIN_ORDER:
        block[f"{c}_cdr3_reads"] = per_chain_prod.get(c, 0)
    return block


def mixcr_total_reads(align_report_json):
    """Last JSON object wins: old re-runs appended blocks (review H9)."""
    try:
        with open(align_report_json) as f:
            text = f.read()
    except OSError:
        return ""
    dec = json.JSONDecoder()
    best = ""
    i = 0
    while i < len(text):
        while i < len(text) and text[i] not in "{[":
            i += 1
        if i >= len(text):
            break
        try:
            obj, end = dec.raw_decode(text, i)
        except ValueError:
            break
        i = end
        if isinstance(obj, dict):
            for key in ("totalReadsProcessed", "totalSequencingReads"):
                if key in obj:
                    best = int(obj[key])
        elif isinstance(obj, list):
            for block_ in obj:
                if isinstance(block_, dict) and block_.get("key") == "totalReadsProcessed":
                    best = int(block_.get("value", 0))
    return best


def generation_of(sample_dir, barcode):
    """Derive the MiXCR generation from what the sentinel actually says."""
    try:
        with open(os.path.join(sample_dir, f"{barcode}.mixcr.done")) as f:
            text = f.read().strip()
    except OSError:
        return "unknown"
    if text.startswith("backfilled"):
        return "pre-cutover"
    if text.startswith("mixcr-v2"):
        return "current"
    return "legacy"


def read_info(sample_dir, barcode):
    """read_length / read_layout written by the mixcr rule (review H2)."""
    path = os.path.join(sample_dir, f"{barcode}.mixcr.readinfo.tsv")
    try:
        with open(path, newline="") as f:
            row = next(csv.reader(f, delimiter="\t"))
        return row[1], row[2]
    except (OSError, StopIteration, IndexError):
        return "", ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-root", required=True, action="append",
                    help="may be repeated: one per shard / results tree")
    ap.add_argument("--manifest", required=True, action="append",
                    help="metadata file(s); may be repeated (one per shard). "
                         "Include has_chimeric for covariate analyses")
    ap.add_argument("--output", required=True)
    ap.add_argument("--rarefy-depth", type=int, default=1000,
                    help="subsample depth for *_rare diversity columns (0 disables)")
    ap.add_argument("--force", action="store_true",
                    help="overwrite an existing output (the QC workflow skips on existence)")
    args = ap.parse_args()

    if os.path.exists(args.output) and not args.force:
        print(f"{args.output} exists — use --force to overwrite", file=sys.stderr)
        return 1

    meta = []
    for path in args.manifest:
        with open(path, newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if "sample" in r and "sample_id" not in r:
                    r["sample_id"] = r["sample"]
                meta.append(r)
    by_sample = {r["sample_id"]: r for r in meta}

    dirs_by_uuid = {}
    for root in args.results_root:
        root_abs = os.path.abspath(root)
        base_depth = len(root_abs.split(os.sep))
        for dirpath, _dirs, files in os.walk(root_abs):
            if len(os.path.abspath(dirpath).split(os.sep)) - base_depth == 3 and any(
                    x.endswith(".unified.tsv") for x in files):
                dirs_by_uuid.setdefault(os.path.basename(dirpath), dirpath)

    out_rows = []
    for uuid, r in by_sample.items():
        d = dirs_by_uuid.get(uuid)
        if not d:
            continue
        bc = r["barcode_28"]
        unified = os.path.join(d, f"{bc}.unified.tsv")
        if not os.path.exists(unified):
            continue
        tools = load_unified(unified)
        mx = tool_block(tools.get("mixcr", []), args.rarefy_depth)
        tr = tool_block(tools.get("trust4", []), args.rarefy_depth)
        read_len, read_layout = read_info(d, bc)
        generation = generation_of(d, bc)
        row = {
            "sample_id": uuid,
            "barcode_28": bc,
            "cancer_type": r["cancer_type"],
            "sample_type": bc.split("-")[3][:2],
            "has_chimeric": r.get("has_chimeric", ""),
            "chimeric_analyzed": ("yes" if generation == "pre-cutover"
                                  and r.get("has_chimeric") == "1" else "no"),
            "read_length": read_len,
            "read_layout": read_layout,
            "bam_size_gb": bam_size_gb(r),
            "mixcr_generation": generation,
            "total_reads": mixcr_total_reads(os.path.join(d, f"{bc}.mixcr.align.report.json")),
            "mixcr_aligned_reads": "",
            "mixcr_aligned_frac": "",
            "mixcr_count_unit": count_unit_for_layout(read_layout),
            "trust4_count_unit": "read",
            "rarefy_depth": args.rarefy_depth,
        }
        txt = os.path.join(d, f"{bc}.mixcr.align.report.txt")
        try:
            with open(txt) as f:
                for line in f:
                    if line.startswith("Total sequencing reads:"):
                        row["total_reads"] = int(line.split(":")[1])
                    elif line.startswith("Successfully aligned reads:"):
                        parts = line.split(":")[1].split("(")
                        row["mixcr_aligned_reads"] = int(parts[0])
                        row["mixcr_aligned_frac"] = round(
                            float(parts[1].rstrip("%)\n")) / 100, 8)
        except (OSError, ValueError, IndexError):
            pass
        for k, v in mx.items():
            row[f"mixcr_{k}"] = v
        row["mixcr_cdr3_frac"] = (
            round(safe_float(mx.get("cdr3_reads")) / row["total_reads"], 10)
            if row["total_reads"] and mx.get("cdr3_reads") != "" else "")
        for k, v in tr.items():
            row[f"trust4_{k}"] = v
        row["trust4_cdr3_frac"] = cdr3_fraction(
            tr.get("cdr3_reads", ""), row["total_reads"], read_layout)
        out_rows.append(row)

    if not out_rows:
        print("no completed samples found", file=sys.stderr)
        return 1
    out_rows.sort(key=lambda r: (r["cancer_type"], r["barcode_28"]))
    with open(args.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()),
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)
    print(f"wrote {args.output}: {len(out_rows)} samples x {len(out_rows[0])} metrics")
    return 0


if __name__ == "__main__":
    sys.exit(main())
