#!/usr/bin/env bash
# watchdog_drivers.sh — keep one driver alive per shard without a human (or an
# SSH session) watching.
#
# The drivers are detached daemons (`setsid`, ppid=1), so they survive a
# closed workstation session; this watchdog covers the remaining failure mode:
# a driver that dies (OOM, node reboot, MAX_CYCLES exhausted) leaves the shard
# idle forever. Run it from cron every 10 minutes on <node>.
#
# Two consecutive misses are required before relaunching, so the relaunch does
# not race a normal restart; a real relaunch cancels the dead generation's
# leftovers (driver_shard.sh does that itself).
#
# Install (root on <node>, as wsx owns the files):
#   echo '*/10 * * * * $HOME/Projects/TCGA-mixcr/scripts/watchdog_drivers.sh' | crontab -
# Remove: crontab -e (delete the line).
set -uo pipefail
BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
LOG=$BASE/logs/watchdog.log
STATE=$BASE/logs/.watchdog_miss
JOB_USER=$(stat -c %U "$BASE/shards/shard1" 2>/dev/null || echo wsx)
mkdir -p "$BASE/logs" "$STATE"

log() { echo "$(date -Is) $*" >> "$LOG"; }

# heartbeat: one line per run so "is the watchdog actually being executed?" can
# be answered from the log (cron does not log user-crontab CMD lines on this
# host — verified 2026-11-17 with a trivial '* * * * * date' entry that ran
# without any journal entry).
log "check drivers=$(pgrep -cf 'driver_shar[d]\.sh [0-9]$')"

for n in 1 2 3 4; do
  if pgrep -f "driver_shar[d]\.sh $n$" >/dev/null 2>&1; then
    rm -f "$STATE/shard$n"
    continue
  fi
  # first miss: record and wait for the next round
  if [ ! -f "$STATE/shard$n" ]; then
    date -Is > "$STATE/shard$n"
    log "shard$n driver missing (first observation)"
    continue
  fi
  jobs=$(squeue -u "$JOB_USER" -h -o "%j" 2>/dev/null | grep -cE '^(mixcr|trust4|standardize|cleanup)$')
  log "shard$n driver still missing after two checks (queue=$jobs) — relaunching"
  cd "$BASE" || exit 1
  sudo -u "$JOB_USER" setsid nohup bash "$BASE/scripts/driver_shard.sh" "$n" \
    >> "$BASE/logs/driver_shard${n}.out" 2>&1 < /dev/null &
  rm -f "$STATE/shard$n"
  sleep 2
done
