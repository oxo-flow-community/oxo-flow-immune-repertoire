#!/usr/bin/env bash
# backfill_sentinels.sh — create completion sentinels for samples finished by
# the PRE-sentinel workflow, so switching to the sentinel workflow does not
# re-run work that is already complete.
#
# Usage:
#   backfill_sentinels.sh <results_root>            # dry run (reports only)
#   backfill_sentinels.sh <results_root> --apply    # write the sentinel files
#
# A sentinel is written only when the corresponding outputs are demonstrably
# complete:
#   mixcr      : <barcode>.mixcr.clones.tsv AND <barcode>.mixcr.qc.json exist
#                (qc.json may be `[]` for a zero-clonotype sample — that is a
#                valid result, so only existence is required)
#   trust4     : cdr3.out + report.tsv + airr.tsv all exist (the first two may
#                be zero-length stubs for zero-clonotype samples by design)
#   standardize: <barcode>.unified.tsv carries the canonical header
#
# Never touches samples that are mid-run: a rule only writes its sentinel here
# when all of its declared outputs are present.
set -uo pipefail

ROOT=${1:?usage: backfill_sentinels.sh <results_root> [--apply]}
APPLY=${2:-}
[ -d "$ROOT" ] || { echo "no such directory: $ROOT" >&2; exit 1; }

n_mixcr=0; n_trust4=0; n_unified=0; n_incomplete=0
while IFS= read -r f; do
  base=${f%.mixcr.clones.tsv}
  if [ -e "${base}.mixcr.qc.json" ]; then
    if [ ! -e "${base}.mixcr.done" ]; then
      n_mixcr=$((n_mixcr + 1))
      [ "$APPLY" = "--apply" ] && echo "backfilled mixcr $(date -Is)" > "${base}.mixcr.done"
    fi
  else
    n_incomplete=$((n_incomplete + 1))
    echo "SKIP mixcr (qc.json missing): $f" >&2
  fi
done < <(find "$ROOT" -name '*.mixcr.clones.tsv' -type f 2>/dev/null)

while IFS= read -r f; do
  base=${f%.trust4_airr.tsv}
  if [ -e "${base}.trust4_cdr3.out" ] && [ -e "${base}.trust4_report.tsv" ]; then
    if [ ! -e "${base}.trust4.done" ]; then
      n_trust4=$((n_trust4 + 1))
      [ "$APPLY" = "--apply" ] && echo "backfilled trust4 $(date -Is)" > "${base}.trust4.done"
    fi
  else
    n_incomplete=$((n_incomplete + 1))
    echo "SKIP trust4 (sibling outputs missing): $f" >&2
  fi
done < <(find "$ROOT" -name '*.trust4_airr.tsv' -type f 2>/dev/null)

while IFS= read -r f; do
  if head -1 "$f" 2>/dev/null | grep -q '^barcode_28'; then
    base=${f%.unified.tsv}
    if [ ! -e "${base}.unified.done" ]; then
      n_unified=$((n_unified + 1))
      [ "$APPLY" = "--apply" ] && echo "backfilled standardize $(date -Is)" > "${base}.unified.done"
    fi
  else
    n_incomplete=$((n_incomplete + 1))
    echo "SKIP unified (no canonical header): $f" >&2
  fi
done < <(find "$ROOT" -name '*.unified.tsv' -type f 2>/dev/null)

mode="DRY-RUN"
[ "$APPLY" = "--apply" ] && mode="APPLIED"
echo "$mode root=$ROOT mixcr=$n_mixcr trust4=$n_trust4 unified=$n_unified incomplete_or_skipped=$n_incomplete"
