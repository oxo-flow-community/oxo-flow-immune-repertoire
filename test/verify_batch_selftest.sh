#!/usr/bin/env bash
# verify_batch_selftest.sh — regression test for scripts/verify_batch.py.
#
# Builds one synthetic sample with the production standardizer, confirms the
# verifier passes, then mutates the output in three ways that MUST be caught:
#   (1) a dropped row        -> conservation failure
#   (2) an altered fraction  -> conservation + fraction-sum failure
#   (3) a missing sentinel   -> completeness failure
# Exits non-zero if any expectation is not met.
set -euo pipefail

BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "${WORK:?}"' EXIT
cd "$WORK"

B=TCGA-AA-0001-01A-01R-1022-07
D=results/COAD/TCGA-AA-0001/11111111-2222-3333-4444-555555555555
mkdir -p "$D"

printf 'cloneId\tallVHitsWithScore\tallJHitsWithScore\tallCHitsWithScore\taaSeqCDR3\tnSeqCDR3\treadCount\tchain\n1\tTRBV28(100,25/262)\tTRBJ2-7(100)\tTRBC2(100)\tCASSLGQAYEQYF\tTGTGCCAGCAGCTTAGGACAGGCCTATGAACAGTACTTC\t4472.0\tTRB\n2\tIGHV1-2(100)\tIGHJ4(100)\tIGHM(100)\tCARDRGYSYGYFDY\tTGTGCGAGAGATAGGGGGTATAGCTATGGTTACTTTGACTAC\t120.0\tIGH\n' > "$D/$B.mixcr.clones.tsv"

printf 'sequence_id\tlocus\tjunction\tjunction_aa\tv_call\tj_call\tc_call\tconsensus_count\nassemble0\tTRB\tTGTGCCAGCAGCTTAGGACAGGCCTATGAACAGTACTTC\tCASSLGQAYEQYF\tTRBV28*01\tTRBJ2-7*01\tTRBC2*01\t3100\nassemble1\tIGH\tTGTGCGAGAGATAGGGGGTATAGCTATGGTTACTTTGACTAC\tCARDRGYSYGYFDY\tIGHV1-2*01\tIGHJ4*01\tIGHM*01\t60\nassemble2\tTRB\tTGTGCCAGCAGCTTAGGACAGGCCTATGAACAGTACTTA\tCASSLGQAYEQYF\tTRBV28*01\tTRBJ2-7*01\tTRBC2*01\t12\n' > "$D/$B.trust4_airr.tsv"

printf '{"inputFiles":{"x":1}}\n' > "$D/$B.mixcr.qc.json"
# header-only stubs, like the rule writes for a zero-clonotype sample
# (an *empty* file is now a completeness failure — review H8)
printf '#count\tfrequency\tCDR3nt\tCDR3aa\tV\tD\tJ\tC\n' > "$D/$B.trust4_cdr3.out"
printf 'count\tfrequency\tCDR3nt\tCDR3aa\tV\tD\tJ\tC\n' > "$D/$B.trust4_report.tsv"
echo "current $(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$D/$B.mixcr.done"
echo ok > "$D/$B.trust4.done"
: > "$D/$B.tmp_cleaned"

standardize() {
  python3 "$BASE/scripts/standardize_clonotypes.py" \
    --mixcr "$D/$B.mixcr.clones.tsv" --trust4-airr "$D/$B.trust4_airr.tsv" \
    --barcode "$B" --cancer COAD --output "$D/$B.unified.tsv" >/dev/null
  echo ok > "$D/$B.unified.done"
}

expect() {
  local want=$1 label=$2
  set +e
  python3 "$BASE/scripts/verify_batch.py" --results results >/dev/null 2>&1
  local rc=$?
  set -e
  if [ "$rc" -ne "$want" ]; then
    echo "FAIL: $label (verify_batch rc=$rc, want $want)"
    exit 1
  fi
  echo "ok: $label"
}

standardize
expect 0 "clean sample passes"

grep -v IGHM "$D/$B.unified.tsv" > "$D/u.tmp" && mv "$D/u.tmp" "$D/$B.unified.tsv"
expect 1 "dropped row is caught"
standardize

python3 - "$D/$B.unified.tsv" <<'PY'
import sys
p = sys.argv[1]
lines = open(p).read().splitlines()
f = lines[1].split("\t")
f[12] = "0.9"
lines[1] = "\t".join(f)
open(p, "w").write("\n".join(lines) + "\n")
PY
expect 1 "altered fraction is caught"
standardize

rm "${D:?}/${B:?}.mixcr.done"
expect 1 "missing sentinel is caught"

echo "all verify_batch self-tests passed"
