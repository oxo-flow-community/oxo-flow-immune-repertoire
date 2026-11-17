#!/usr/bin/env python3
"""check_tsv_integrity.py — is this TSV structurally complete?

A killed job leaves a file that exists but is truncated: the last line is cut
mid-record, or a later line has fewer columns than the header. "File exists"
is therefore not evidence of completeness (scientific review H8), and both
backfill_sentinels.sh and driver_shard.sh used existence as their test.

Checks: header present, every row has the header's column count, the file ends
with a newline, and (optionally) the row count reaches a minimum.

Exit: 0 complete, 1 incomplete, 2 unreadable.
Usage: check_tsv_integrity.py <file.tsv> [--min-rows N] [--quiet]
"""
import argparse
import sys


def check(path, min_rows=0):
    try:
        with open(path, "rb") as f:
            raw = f.read()
    except OSError as exc:
        return False, f"unreadable: {exc}"
    if not raw:
        return False, "empty file"
    if not raw.endswith(b"\n"):
        return False, "truncated: no trailing newline"
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        return False, f"invalid utf-8: {exc}"
    lines = text.splitlines()
    header_cols = lines[0].count("\t")
    for i, line in enumerate(lines, 1):
        cols = line.count("\t")
        if cols != header_cols:
            return False, f"line {i}: {cols + 1} columns, header has {header_cols + 1}"
    rows = len(lines) - 1
    if rows < min_rows:
        return False, f"{rows} rows < required {min_rows}"
    return True, f"{rows} rows x {header_cols + 1} columns"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--min-rows", type=int, default=0)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    ok, detail = check(args.path, args.min_rows)
    if not args.quiet:
        print(("OK   " if ok else "BAD  ") + f"{args.path}: {detail}")
    return 0 if ok else (2 if detail.startswith("unreadable") else 1)


if __name__ == "__main__":
    sys.exit(main())
