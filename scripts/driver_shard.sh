#!/usr/bin/env bash
# driver_shard.sh — production driver for one TCGA-mixcr shard using the
# oxo-flow Slurm **cluster backend** (replaces supervisor.sh / sbatch wrapper).
#
# Usage:  driver_shard.sh <shardN>          (N = 1..4)
# Launch: setsid nohup driver_shard.sh 1 >/dev/null 2>&1 &
#
# Behaviour:
#   * one instance per shard (flock); duplicate launches exit immediately
#   * runs `oxo-flow run ... --profile slurm -k` so a failing sample never
#     cancels the rest of the shard (keep-going)
#   * a restart re-runs every rule whose declared outputs are missing/stale
#     (validated T3e: recorded PENDING state is re-queued, not adopted).
#     Safety against "partial output mistaken for completion" (T3b–T3d,
#     docs/oxoflow_issue_restart_silent_loss.md) comes from the workflow's
#     LAST-written sentinel outputs, so output-existence implies completion.
#   * leftover jobs from the previous run are cancelled and the driver waits
#     for them to leave the queue before re-submitting the same rules.
#   * an output-count guard turns "rc=0 but samples missing" into a retry
#   * exponential backoff between restarts (5 min -> 30 min cap)
#   * all output appended to logs/driver_shardN.log
#
# CPU policy: per-rule resources are 8 threads and the shard profile caps
# max_submitted at 2 => 16 threads/node (hard policy: <= 20).
set -uo pipefail

SHARD=${1:?usage: driver_shard.sh <1|2|3|4>}
BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
DIR=$BASE/shards/shard$SHARD
LOG=$BASE/logs/driver_shard$SHARD.log
LOCK=$BASE/logs/.driver_shard$SHARD.lock
OXO=oxo-flow

[ -d "$DIR" ] || { echo "missing shard dir: $DIR" >&2; exit 1; }
[ -x "$OXO" ] || { echo "missing oxo-flow binary: $OXO" >&2; exit 1; }

exec 9>"$LOCK"
flock -n 9 || { echo "driver already running for shard$SHARD" >&2; exit 0; }

cd "$DIR" || exit 1
# oxo-flow must be on PATH; override OXO=... if it lives elsewhere
export PATH="$HOME/.cargo/bin:$PATH"

# Expected number of samples in this shard (groups.tsv: name<TAB>comma-list).
EXPECTED=$(awk -F'\t' 'NR>1 { s += split($2, a, ",") } END { print s + 0 }' manifests/groups.tsv)

# Cancel EVERY live rule job that belongs to this shard, not merely the ids
# recorded in the newest run directory. A driver restart (or a workflow
# reconfiguration) leaves the previous generation's Slurm jobs running; if
# they are not cancelled they overlap the new generation and double the
# node's CPU footprint (incident 2026-11-16: 32 CPUs/node vs the 20 policy).
# The job's Slurm WorkDir is the shard directory, which makes the match
# unambiguous even when several shards run concurrently.
cancel_leftovers() {
  local id wd n=0
  for id in $(squeue -u "$USER" -h -o "%i %j" 2>/dev/null \
      | awk '$2 ~ /^(mixcr|trust4|standardize|cleanup)/ {print $1}'); do
    wd=$(scontrol show job "$id" 2>/dev/null | sed -n 's/.*WorkDir=\([^ ]*\).*/\1/p' | head -1)
    [ "$wd" = "$PWD" ] || continue
    scancel "$id" 2>/dev/null && n=$((n + 1))
  done
  [ "$n" -gt 0 ] && echo "$(date -Is) cancelled $n leftover rule job(s) from a previous generation" >> "$LOG"
  for _ in $(seq 1 30); do
    [ -z "$(squeue -u "$USER" -h -o "%j" 2>/dev/null | grep -E '^(mixcr|trust4|standardize|cleanup)')" ] && break
    sleep 2
  done
}

BACKOFF=300
MAX_BACKOFF=1800
while true; do
  # Cancel jobs the previous run left behind, then wait for them to leave the
  # queue. Killed rules are re-submitted automatically (T3e: recorded PENDING
  # state is NOT adopted as complete), so an orphan job would race the fresh
  # submission and both would write the same MiXCR prefix.
  cancel_leftovers

  echo "$(date -Is) driver shard$SHARD start (expected_samples=$EXPECTED)" >> "$LOG"
  # -k / --retry / --max-threads are NOT honoured on the cluster path (the
  # engine warns about each): the cluster driver already keeps dispatching
  # independent jobs until the queue drains, and retry/backoff is this
  # wrapper's. Passing them only adds warnings to the log.
  "$OXO" run tcga-mixcr.oxoflow --profile slurm >> "$LOG" 2>&1
  rc=$?

  ACTUAL=$(find results -name '*.unified.tsv' 2>/dev/null | wc -l)
  echo "$(date -Is) driver shard$SHARD exit rc=$rc unified=${ACTUAL}/${EXPECTED}" >> "$LOG"
  if [ "$rc" -eq 0 ] && [ "$ACTUAL" -lt "$EXPECTED" ]; then
    echo "$(date -Is) OUTPUT GAP: $((EXPECTED - ACTUAL)) sample(s) without unified.tsv — retrying" >> "$LOG"
    rc=42
  fi
  if [ "$rc" -eq 0 ]; then
    echo "$(date -Is) driver shard$SHARD finished cleanly" >> "$LOG"
    exit 0
  fi

  echo "$(date -Is) driver shard$SHARD restart in ${BACKOFF}s" >> "$LOG"
  sleep "$BACKOFF"
  BACKOFF=$(( BACKOFF * 2 ))
  [ "$BACKOFF" -gt "$MAX_BACKOFF" ] && BACKOFF=$MAX_BACKOFF
done
