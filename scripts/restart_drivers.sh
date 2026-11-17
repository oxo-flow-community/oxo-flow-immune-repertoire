#!/usr/bin/env bash
# restart_drivers.sh — stop the running shard drivers and relaunch them so a
# profile change (walltime, max_submitted) takes effect.
#
# Why this is a script and not a one-liner: the external SSH path (nps) drops
# long/interactive sessions at random. A multi-step command pasted inline can
# die half-way with no record of what ran (AGENTS.md §2.2). The script logs
# every action so the state is recoverable after a dropped connection.
#
# Order matters:
#   1. kill the bash drivers first (else they would notice their child dying
#      and start a retry/backoff branch)
#   2. kill the oxo-flow run children — SIGKILL, the driver does not always
#      honour SIGTERM (AGENTS.md §5.4)
#   3. cancel every rule job (matched by job name; the drivers match by
#      WorkDir, but here all shards are being restarted)
#   4. relaunch all four drivers detached
#
# Usage: restart_drivers.sh [--dry-run]
set -uo pipefail

BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
LOG=$BASE/logs/restart_drivers.log
DRY=0
[ "${1:-}" = "--dry-run" ] && DRY=1

log() { echo "$(date -Is) $*" | tee -a "$LOG"; }

cd "$BASE" || exit 1
JOB_USER=${JOB_USER:-$(stat -c %U shards/shard1 2>/dev/null || echo "$USER")}

# Guard: refuse to relaunch with a walltime too short for the big samples.
short=$(grep -L 'walltime *= *"48h"' shards/shard*/profiles/slurm.toml)
if [ -n "$short" ]; then
  log "ABORT: profile(s) without walltime 48h: $short"
  exit 2
fi

driver_pids=$(pgrep -f 'driver_shar[d]\.sh [0-9]$' || true)
oxo_pids=$(pgrep -f 'oxo-flow[ ]run tcga-mixcr\.oxoflow' || true)
log "found drivers: ${driver_pids:-none}; oxo-flow runs: ${oxo_pids:-none}"

if [ "$DRY" -eq 1 ]; then
  log "dry-run: would kill the above, scancel rule jobs, and relaunch shards 1-4"
  exit 0
fi

if [ -n "$driver_pids" ]; then
  kill -9 $driver_pids 2>/dev/null && log "killed driver(s): $driver_pids"
fi
sleep 2
if [ -n "$oxo_pids" ]; then
  kill -9 $oxo_pids 2>/dev/null && log "killed oxo-flow run(s): $oxo_pids"
fi
sleep 3

left=$(pgrep -f 'driver_shar[d]\.sh|oxo-flow[ ]run tcga-mixcr\.oxoflow' || true)
if [ -n "$left" ]; then
  log "ABORT: processes survived SIGKILL: $left"
  exit 3
fi
log "all driver processes stopped"

n=0
for id in $(squeue -u "$JOB_USER" -h -o "%i %j" \
    | awk '$2 ~ /^(mixcr|trust4|standardize|cleanup)/ {print $1}'); do
  scancel "$id" && n=$((n + 1))
done
log "cancelled $n rule job(s)"

for _ in $(seq 1 30); do
  [ "$(squeue -u "$JOB_USER" -h -o "%j" | grep -cE '^(mixcr|trust4|standardize|cleanup)')" -eq 0 ] && break
  sleep 2
done
still=$(squeue -u "$JOB_USER" -h -o "%j" | grep -cE '^(mixcr|trust4|standardize|cleanup)')
if [ "$still" -ne 0 ]; then
  log "ABORT: $still rule job(s) still in the queue after 60s"
  exit 4
fi
log "queue clear of rule jobs; relaunching drivers"

for s in 1 2 3 4; do
  setsid nohup bash "$BASE/scripts/driver_shard.sh" "$s" \
    >> "$BASE/logs/driver_shard${s}.out" 2>&1 < /dev/null &
  log "launched driver_shard.sh $s (pid $!)"
done

sleep 20
log "now running: $(pgrep -cf 'driver_shar[d]\.sh [0-9]$') driver(s), $(squeue -u "$JOB_USER" -h -o "%j" | grep -cE '^(mixcr|trust4|standardize|cleanup)') rule job(s)"
