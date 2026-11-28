#!/usr/bin/env python3
"""Self-test for scripts/sample_metrics.py::is_productive.

Regression for the defect found by independent review on 2026-11-23: the
filter rejected stop codons (`*`) but accepted MiXCR frameshift gaps (`_`)
and did not check the nucleotide frame, so an out-of-frame clone such as
`CALPAL_SKLTF` (35-nt CDR3) was counted as productive in clonality/gini.

Run: python3 -I test/sample_metrics_selftest.py
"""
import importlib.util
import os
import pathlib
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parent.parent
_SPEC = importlib.util.spec_from_file_location(
    "sample_metrics", ROOT / "scripts" / "sample_metrics.py")
sm = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(sm)

INTACT_NT = "A" * 36       # codon-complete CDR3 nucleotide
FRAMESHIFT_NT = "A" * 35   # one base over a codon


def row(aa, nt=""):
    return {"cdr3_aa": aa, "cdr3_nt": nt}


CASES = [
    # name, row, expected
    ("in-frame MiXCR clone", row("CARDSTIVATGGGRSGMDVW", INTACT_NT), True),
    ("TRUST4 productive row", row("CASSLGQAYEQYF", INTACT_NT), True),
    ("no cdr3_nt but clean aa", row("CARDW", ""), True),
    ("stop codon", row("CARD*VW", INTACT_NT), False),
    ("frameshift gap, no stop", row("CALPAL_SKLTF", FRAMESHIFT_NT), False),
    ("gap inside an otherwise clean aa", row("CAR_DW", INTACT_NT), False),
    ("clean aa but out-of-frame nt", row("CARDW", FRAMESHIFT_NT), False),
    ("empty aa", row("", INTACT_NT), False),
    ("does not start with C", row("XARDW", INTACT_NT), False),
    ("does not end with F/W", row("CARDV", INTACT_NT), False),
]


def check(failures, name, got, expected):
    status = "PASS" if got == expected else "FAIL"
    if got != expected:
        failures += 1
    print(f"{status}  {name}: expected {expected}, got {got}")
    return failures


def helper_cases():
    """The 2026-11-27 additions: per-chain reads, units, fractions, BAM size."""
    failures = 0
    rows = [
        {"chain": "TRB", "cdr3_aa": "CASSLGQAYEQYF", "cdr3_nt": INTACT_NT,
         "read_count": "10"},
        {"chain": "TRB", "cdr3_aa": "CALPAL_SKLTF", "cdr3_nt": FRAMESHIFT_NT,
         "read_count": "5"},
        {"chain": "IGH", "cdr3_aa": "CARDW", "cdr3_nt": INTACT_NT,
         "read_count": "7"},
        {"chain": "IGH", "cdr3_aa": "", "cdr3_nt": "", "read_count": "2"},
    ]
    block = sm.tool_block(rows, 1000)
    # productive-only per-chain reads
    failures = check(failures, "TRB productive reads", block["TRB_cdr3_reads"], 10)
    failures = check(failures, "IGH productive reads", block["IGH_cdr3_reads"], 7)
    # read-weighted fraction still includes every row (15/24)
    failures = check(failures, "TRB all-row fraction", round(block["frac_TRB"], 4),
                     round(15 / 24, 4))
    failures = check(failures, "PE unit", sm.count_unit_for_layout("paired-end"),
                     "read_pair")
    failures = check(failures, "SE unit", sm.count_unit_for_layout("single-end"), "read")
    # unknown layout falls back to read_pair, matching standardize
    failures = check(failures, "unknown unit fallback",
                     sm.count_unit_for_layout(""), "read_pair")
    failures = check(failures, "PE trust4 fraction",
                     sm.cdr3_fraction(100, 1000, "paired-end"), 0.05)
    failures = check(failures, "SE trust4 fraction",
                     sm.cdr3_fraction(100, 1000, "single-end"), 0.1)
    failures = check(failures, "unknown-layout fraction",
                     sm.cdr3_fraction(100, 1000, ""), "")
    failures = check(failures, "missing reads fraction",
                     sm.cdr3_fraction("", 1000, "paired-end"), "")
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "x.bam")
        with open(path, "wb") as f:
            f.truncate(123_456_789)  # sparse file: stat size without disk use
        failures = check(failures, "bam size from stat",
                         sm.bam_size_gb({"bam_genomic": path}), 0.12)
        failures = check(failures, "bam size from manifest",
                         sm.bam_size_gb({"file_size_bytes": "2000000000"}), 2.0)
        failures = check(failures, "bam size NA falls back to stat",
                         sm.bam_size_gb({"file_size_bytes": "NA",
                                         "bam_genomic": path}), 0.12)
    return failures


def main():
    failures = 0
    for name, sample_row, expected in CASES:
        failures = check(failures, name, sm.is_productive(sample_row), expected)
    failures += helper_cases()
    total = len(CASES) + 13
    print(f"{total - failures}/{total} cases pass")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
