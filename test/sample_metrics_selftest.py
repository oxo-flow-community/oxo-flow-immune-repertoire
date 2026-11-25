#!/usr/bin/env python3
"""Self-test for scripts/sample_metrics.py::is_productive.

Regression for the defect found by independent review on 2026-11-23: the
filter rejected stop codons (`*`) but accepted MiXCR frameshift gaps (`_`)
and did not check the nucleotide frame, so an out-of-frame clone such as
`CALPAL_SKLTF` (35-nt CDR3) was counted as productive in clonality/gini.

Run: python3 -I test/sample_metrics_selftest.py
"""
import importlib.util
import pathlib
import sys

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


def main():
    failures = 0
    for name, sample_row, expected in CASES:
        got = sm.is_productive(sample_row)
        status = "PASS" if got == expected else "FAIL"
        if got != expected:
            failures += 1
        print(f"{status}  {name}: expected {expected}, got {got}")
    print(f"{len(CASES) - failures}/{len(CASES)} cases pass")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
