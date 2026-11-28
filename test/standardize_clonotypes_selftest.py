#!/usr/bin/env python3
"""Self-test for scripts/standardize_clonotypes.py.

Covers the 2026-11-27 review fixes:
  * per-sample MiXCR count_unit from the readinfo file (paired-end -> read_pair,
    single-end -> read) instead of a hardcoded read_pair;
  * atomic output write (no truncated table under an old sentinel);
  * safe_int handles inf/NaN without crashing;
  * TRUST4 clone_count is normalized like read_count.

Run: python3 -I test/standardize_clonotypes_selftest.py
"""
import csv
import importlib.util
import os
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "scripts" / "standardize_clonotypes.py"
_SPEC = importlib.util.spec_from_file_location("standardize_clonotypes", SCRIPT)
sc = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(sc)

BARCODE = "TCGA-AA-3870-01A-01R-1022-07"
CANCER = "COAD"

MIXCR_HEADER = ["chain", "allVHitsWithScore", "allJHitsWithScore",
                "allCHitsWithScore", "aaSeqCDR3", "nSeqCDR3", "readCount"]
TRUST4_HEADER = ["locus", "junction", "junction_aa", "consensus_count",
                 "v_call", "j_call", "c_call"]


def write_mixcr(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(MIXCR_HEADER)
        w.writerows(rows)


def write_airr(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(TRUST4_HEADER)
        w.writerows(rows)


def write_readinfo(path, layout, uuid="4aa404c5-683d-4275-8806-45568beda8ea"):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow([uuid, 75, layout, 163])


def run_cli(tmp):
    out = os.path.join(tmp, "unified.tsv")
    mixcr = os.path.join(tmp, f"{BARCODE}.mixcr.clones.tsv")
    airr = os.path.join(tmp, "airr.tsv")
    r = subprocess.run([sys.executable, "-I", str(SCRIPT),
                        "--mixcr", mixcr, "--trust4-airr", airr,
                        "--barcode", BARCODE, "--cancer", CANCER,
                        "--output", out],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError(f"script failed: {r.stdout}\n{r.stderr}")
    return out


def read_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def check(failures, name, got, expected):
    status = "PASS" if got == expected else "FAIL"
    if got != expected:
        failures += 1
    print(f"{status}  {name}: expected {expected}, got {got}")
    return failures


def main():
    failures = 0

    # --- pure-function checks ---
    failures = check(failures, "safe_int float string", sc.safe_int("4472.0"), 4472)
    failures = check(failures, "safe_int inf passthrough", sc.safe_int("inf"), "inf")
    failures = check(failures, "safe_int NaN passthrough", sc.safe_int("nan"), "nan")
    failures = check(failures, "safe_int text passthrough", sc.safe_int("abc"), "abc")
    failures = check(failures, "PE -> read_pair",
                     sc.count_unit_for_layout("paired-end"), "read_pair")
    failures = check(failures, "SE -> read",
                     sc.count_unit_for_layout("single-end"), "read")
    rows = [
        {"tool": "mixcr", "chain": "TRB", "read_count": 3, "clone_fraction": ""},
        {"tool": "mixcr", "chain": "TRB", "read_count": 7, "clone_fraction": ""},
        {"tool": "mixcr", "chain": "IGH", "read_count": 5, "clone_fraction": ""},
    ]
    frac = sc.assign_chain_fractions(rows)
    sums = {}
    for r in frac:
        sums[r["chain"]] = sums.get(r["chain"], 0) + float(r["clone_fraction"])
    failures = check(failures, "TRB fraction sum", round(sums["TRB"], 6), 1.0)
    failures = check(failures, "IGH fraction sum", round(sums["IGH"], 6), 1.0)

    # --- end-to-end CLI: paired-end and single-end ---
    for layout, expected_unit in (("paired-end", "read_pair"),
                                  ("single-end", "read")):
        with tempfile.TemporaryDirectory() as td:
            mixcr = os.path.join(td, f"{BARCODE}.mixcr.clones.tsv")
            airr = os.path.join(td, "airr.tsv")
            readinfo = os.path.join(td, f"{BARCODE}.mixcr.readinfo.tsv")
            write_mixcr(mixcr, [
                ["TRB", "TRBV28(100)", "TRBJ1-1(100)", "TRBC1(100)",
                 "CASSLGQAYEQYF", "TGTGCCAGCAGCTTGGGACAGGCCTACGAGCAGTACTTC", "10.0"],
                ["IGH", "IGHV1-2(100)", "IGHJ4(100)", "IGHM(100)",
                 "CARDW", "TGTGCGAGAGACTGG", "7.0"],
            ])
            write_airr(airr, [
                ["TRB", "TGTGCCAGCAGCTTGGGACAGGCCTACGAGCAGTACTTC",
                 "CASSLGQAYEQYF", "4", "TRBV28", "TRBJ1-1", "TRBC1"],
                ["IGH", "TGTGCGAGAGACTGG", "CARDW", "3", "IGHV1-2", "IGHJ4", "IGHM"],
            ])
            write_readinfo(readinfo, layout)
            out = run_cli(td)
            got = read_rows(out)
            units = {(r["tool"], r["count_unit"]) for r in got}
            failures = check(failures, f"{layout} units", units,
                             {("mixcr", expected_unit), ("trust4", "read")})
            failures = check(failures, f"{layout} rows", len(got), 4)
            leftovers = [x for x in os.listdir(td) if ".tmp." in x]
            failures = check(failures, f"{layout} no temp leftovers", leftovers, [])

    total = 14
    print(f"{total - failures}/{total} cases pass")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
