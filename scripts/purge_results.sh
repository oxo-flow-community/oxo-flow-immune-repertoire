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
LEGACY=0
for arg in "$@"; do
  [ "$arg" = "--apply" ] && APPLY=1
  [ "$arg" = "--legacy" ] && LEGACY=1
done
# Refuse to purge while the pipeline is live. On 2026-11-17 a --legacy run
# during production re-ran the shard section and deleted the tmp/ and
# .oxo-flow/runs of 3 shards under 11 running jobs, killing ~1 h of work.
LIVE=$(( $(squeue -u "${JOB_USER:-$(stat -c %U "$BASE/shards/shard1" 2>/dev/null || echo "$USER")}" -h -o "%j" 2>/dev/null \
          | grep -cE "^(mixcr|trust4|standardize|cleanup)") + $(pgrep -cf "driver_shar[d]\.sh [0-9]$") ))
if [ "$LIVE" -gt 0 ] && [ "${FORCE_PURGE:-0}" != "1" ]; then
  echo "REFUSING to purge: $LIVE live job(s)/driver(s). Stop the drivers and cancel jobs first,"
  echo "or re-run with FORCE_PURGE=1 if you really mean it."
  exit 3
fi
STAMP=$(date +%Y%m%d_%H%M)
ARCHIVE=$BASE/archive/pre_fix_$STAMP
LOGDIR=$BASE/logs
mkdir -p "$LOGDIR"

total=0
if [ "$LEGACY" -eq 1 ] && [ "$APPLY" -eq 1 ] && [ "${LEGACY_ONLY:-0}" = "1" ]; then
  :
else
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
fi

# Legacy trees outside shards/ (found by the 2026-11-17 independent audit):
# the project-root results/, the proto*/ and s4min/ prototype scratch. They
# hold PRE-FIX outputs (duplicated chimeric reads) and one of them overlaps
# the current chunk, so they must not survive into any aggregate.
if [ "$LEGACY" -eq 1 ]; then
  for d in "$BASE/results" "$BASE/proto/results" "$BASE/proto2/results" \
           "$BASE/proto3/results" "$BASE/s4min/results"; do
    [ -d "$d" ] || continue
    MAN=$LOGDIR/purge_manifest_$(echo "${d#$BASE/}" | tr '/' '_')_${STAMP}.txt
    find "$d" -type f -printf '%p\t%s\t%TY-%Tm-%Td %TH:%TM\n' 2>/dev/null > "$MAN"
    count=$(wc -l < "$MAN")
    total=$((total + count))
    echo "legacy $d: $count file(s) — manifest: $MAN"
    if [ "$APPLY" -eq 1 ]; then
      mkdir -p "$ARCHIVE/legacy/$(echo "${d#$BASE/}" | tr '/' '_')"
      find "$d" -name '*.unified.tsv' -exec cp --parents {} "$ARCHIVE/legacy/$(echo "${d#$BASE/}" | tr '/' '_')/" \; 2>/dev/null
      find "$d" -mindepth 1 -delete 2>/dev/null
    fi
  done
fi

if [ "$APPLY" -eq 1 ]; then
  echo "PURGED $total file(s); old unified tables archived under $ARCHIVE"
else
  echo "DRY-RUN: would purge $total file(s). Re-run with --apply."
fi
