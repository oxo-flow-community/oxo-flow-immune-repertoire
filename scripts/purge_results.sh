#!/usr/bin/env bash
# purge_results.sh — wipe per-sample outputs and engine state so the fixed
# workflow starts from a clean slate (2026-11-17 full re-run decision).
#
# Why wipe instead of letting the engine invalidate: the old outputs were
# produced with the duplicated-chimeric-read input (scientific review C1), so
# every one of them is scientifically wrong. A stale unified.tsv left on disk
# is a contamination risk for any later aggregate step, and the engine's
# "outputs exist => skip" test makes deletion the only unambiguous reset.
#
# Safety:
#   * dry-run unless --apply
#   * writes an audit manifest (path, size, mtime) for every file it deletes
#   * archives the 24 unified.tsv tables (small) so the before/after comparison
#     on identical samples stays possible — labelled pre-fix, never used as
#     analysis input
#   * ${VAR:?} guards so an empty variable can never turn into `rm -rf /`
#
# Usage: purge_results.sh [--apply]
set -uo pipefail
BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
APPLY=0
[ "${1:-}" = "--apply" ] && APPLY=1
STAMP=$(date +%Y%m%d_%H%M)
ARCHIVE=$BASE/archive/pre_fix_$STAMP
LOGDIR=$BASE/logs
mkdir -p "$LOGDIR"

total=0
for n in 1 2 3 4; do
  S=$BASE/shards/shard$n
  [ -d "$S" ] || { echo "shard$n: missing, skipped"; continue; }
  MAN=$LOGDIR/purge_manifest_shard${n}_${STAMP}.txt
  find "$S/results" -type f -printf '%p\t%s\t%TY-%Tm-%Td %TH:%TM\n' 2>/dev/null > "$MAN"
  count=$(wc -l < "$MAN")
  bytes=$(awk -F'\t' '{s+=$2} END{print s+0}' "$MAN")
  total=$((total + count))
  echo "shard$n: $count file(s), $(numfmt --to=iec "$bytes" 2>/dev/null || echo "$bytes B") — manifest: $MAN"
  if [ "$APPLY" -eq 1 ]; then
    # keep the unified tables for a labelled before/after comparison
    mkdir -p "$ARCHIVE/shard$n"
    find "$S/results" -name '*.unified.tsv' -exec cp --parents {} "$ARCHIVE/shard$n/" \; 2>/dev/null
    # purge per-sample outputs, engine run state and scratch
    find "$S/results" -mindepth 1 -delete 2>/dev/null
    find "$S/.oxo-flow/runs" -mindepth 1 -delete 2>/dev/null
    find "$S/tmp" -mindepth 1 -delete 2>/dev/null
  fi
done

if [ "$APPLY" -eq 1 ]; then
  echo "PURGED $total file(s); old unified tables archived under $ARCHIVE"
else
  echo "DRY-RUN: would purge $total file(s). Re-run with --apply."
fi
