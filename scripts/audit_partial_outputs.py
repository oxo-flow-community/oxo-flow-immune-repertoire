#!/usr/bin/env python3
"""audit_partial_outputs.py — find (and optionally remove) per-sample outputs
left behind by killed jobs, so they cannot be mistaken for completed work.

A rule is complete **only** when its LAST-written sentinel exists. The engine's
skip test is "all declared outputs exist", so a job killed after writing a
data file but before its sentinel would otherwise be re-run anyway — but any
partial intermediate (`.vdjca`, `.clna`, half-written `unified.tsv`) still
sits in the sample directory and can confuse a re-run or a later merge.
This tool makes that state explicit.

Per sample directory, the rules declare:
  mixcr        clones.tsv + qc.json                    sentinel .mixcr.done
  trust4       trust4_cdr3.out + trust4_report.tsv
               + trust4_airr.tsv                       sentinel .trust4.done
  standardize  unified.tsv                             sentinel .unified.done
  cleanup      tmp_cleaned                             (no sentinel)

Usage:
  python3 scripts/audit_partial_outputs.py --results shards/shard*/results
  python3 scripts/audit_partial_outputs.py --results ... --fix
Exit: 0 = no partial outputs (or all fixed), 1 = partial outputs found.
Stdlib only.
"""
import argparse
import glob
import os
import sys

# rule -> (sentinel, [outputs], [intermediates to remove on --fix])
RULES = {
    "mixcr": (".mixcr.done", ["{b}.mixcr.clones.tsv", "{b}.mixcr.qc.json"],
              ["{b}.mixcr.alignments.vdjca", "{b}.mixcr.partial_assemblies.clna",
               "{b}.mixcr.contigs.clns", "{b}.mixcr.clones_*.tsv",
               "{b}.mixcr.align.report.txt", "{b}.mixcr.align.report.json",
               "{b}.mixcr.align.list.tsv"]),
    "trust4": (".trust4.done", ["{b}.trust4_cdr3.out", "{b}.trust4_report.tsv",
                                "{b}.trust4_airr.tsv"],
               ["{b}.trust4_final.out", "{b}.trust4_annot.fa",
                "{b}.trust4_assembled.fa", "{b}.trust4_bam-report.tsv"]),
    "standardize": (".unified.done", ["{b}.unified.tsv"], []),
}


def audit_sample(d, barcode):
    findings = []
    for rule, (_sentinel, outputs, _inter) in RULES.items():
        sentinel = os.path.join(d, barcode + _sentinel)
        present = [o.format(b=barcode) for o in outputs
                   if os.path.isfile(os.path.join(d, o.format(b=barcode)))]
        if os.path.isfile(sentinel) or not present:
            continue
        findings.append((rule, present))
    return findings


def fix_sample(d, barcode, findings):
    removed = []
    for rule, _present in findings:
        _sentinel, outputs, intermediates = RULES[rule]
        patterns = list(outputs) + list(intermediates)
        for pat in patterns:
            for path in glob.glob(os.path.join(d, pat.format(b=barcode))):
                try:
                    os.remove(path)
                    removed.append(os.path.basename(path))
                except OSError as exc:
                    print(f"  WARN cannot remove {path}: {exc}", file=sys.stderr)
    return removed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", nargs="+", required=True)
    ap.add_argument("--fix", action="store_true",
                    help="delete the partial outputs so the rules re-run cleanly")
    args = ap.parse_args()

    dirs = set()
    for root in args.results:
        for pattern in ("*.mixcr.clones.tsv", "*.trust4_airr.tsv", "*.unified.tsv"):
            for path in glob.glob(os.path.join(root, "**", pattern), recursive=True):
                dirs.add(os.path.dirname(path))

    partial = 0
    fixed = 0
    for d in sorted(dirs):
        files = os.listdir(d)
        barcodes = {os.path.basename(f).split(".mixcr.")[0] for f in files
                    if ".mixcr." in f}
        barcodes |= {os.path.basename(f).split(".trust4")[0] for f in files
                     if ".trust4" in f}
        barcodes |= {os.path.basename(f).split(".unified.")[0] for f in files
                     if ".unified." in f}
        for barcode in sorted(barcodes):
            findings = audit_sample(d, barcode)
            if not findings:
                continue
            partial += 1
            rules = ",".join(r for r, _ in findings)
            print(f"PARTIAL {d}  [{rules}]  files={[f for _, p in findings for f in p]}")
            if args.fix:
                removed = fix_sample(d, barcode, findings)
                fixed += 1
                print(f"  removed {len(removed)} partial file(s)")

    print(f"\n{len(dirs)} sample dir(s) scanned, {partial} with partial outputs"
          + (f", {fixed} cleaned" if args.fix else ""))
    if partial and not args.fix:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
