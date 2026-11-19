#!/usr/bin/env bash
# check_stray_jobs.sh — find leftover interactive allocations in the queue.
#
# WHY (2026-11-19): two root jobs named `bash` (1833, 1834 — 8 + 16 CPUs,
# 32 + 48 GB, 6-day time limit) sat on <node> for 2 days 3 hours. They were
# leftovers from the 2026-11-17 C1 A/B verification: an interactive allocation
# whose shell had exited, so the cgroup held no processes while the allocation
# stayed reserved. Harmless thermally (<node>'s load was 22.7 while our own jobs
# already accounted for ~19.5 cores, so those 24 CPUs were idle) but they held
# a third of <node>'s CPU and 80 GB that other jobs could have used, and they
# made `squeue` look like someone else was competing for the cluster.
#
# RULE: an interactive allocation is for a session, not for a step. Release it
# (`exit`, `scancel`) as soon as the step is done; after any `salloc`/`srun`
# work, run this script and confirm the queue holds only pipeline jobs.
#
# Usage: check_stray_jobs.sh [--fix]
#   --fix   scancel the stray jobs (only ever cancel your own)
# Exit: 0 clean, 1 strays found.
set -uo pipefail
FIX=0
[ "${1:-}" = "--fix" ] && FIX=1

# job names that are only ever suspicious when they linger: a bare shell, or
# the default `srun`/"interactive" names. Pipeline jobs are named after rules.
STRAY='^(bash|sh|zsh|fish|srun|interactive|sleep)$'
# a legitimate short interactive step finishes well inside this window
MAX_MIN=${MAX_MIN:-60}

now=$(date +%s)
strays=0
while read -r jid user name secs; do
  [ -z "$jid" ] && continue
  name=${name%% *}
  echo "$name" | grep -qE "$STRAY" || continue
  if [ "${secs:-0}" -gt $((MAX_MIN * 60)) ]; then
    hrs=$(awk -v s="$secs" 'BEGIN{printf "%.1f", s/3600}')
    echo "STRAY  job=$jid user=$user name=$name running=${hrs}h (> ${MAX_MIN}m) — an interactive allocation left behind?"
    strays=$((strays + 1))
    if [ "$FIX" -eq 1 ]; then
      scancel "$jid" && echo "       cancelled $jid"
    fi
  fi
done < <(squeue -h -o "%i %u %j %S" 2>/dev/null | awk '$4 ~ /^[0-9]+$/ {print $1, $2, $3, $4}')

# squeue's %S is a start time, so re-derive elapsed from %M instead — the awk
# above is a fallback; this is the authoritative pass.
if [ "$strays" -eq 0 ]; then
  while read -r jid user name elapsed; do
    [ -z "$jid" ] && continue
    echo "$name" | grep -qE "$STRAY" || continue
    case "$elapsed" in
      *-*) secs=$(( 10#${elapsed%%-*} * 86400 )); rest=${elapsed#*-}
           IFS=: read -r h m s <<< "$rest"; secs=$(( secs + 10#${h:-0}*3600 + 10#${m:-0}*60 + 10#${s:-0} )) ;;
      *)   IFS=: read -r h m s <<< "$elapsed"
           secs=$(( 10#${h:-0}*3600 + 10#${m:-0}*60 + 10#${s:-0} )) ;;
    esac
    if [ "$secs" -gt $((MAX_MIN * 60)) ]; then
      hrs=$(awk -v s="$secs" 'BEGIN{printf "%.1f", s/3600}')
      echo "STRAY  job=$jid user=$user name=$name running=${hrs}h (> ${MAX_MIN}m)"
      strays=$((strays + 1))
      [ "$FIX" -eq 1 ] && { scancel "$jid" && echo "       cancelled $jid"; }
    fi
  done < <(squeue -h -o "%i %u %j %M" 2>/dev/null)
fi

if [ "$strays" -eq 0 ]; then
  echo "clean: no leftover interactive allocations"
  exit 0
fi
echo "$strays stray allocation(s) found$([ "$FIX" -eq 1 ] || echo ' — re-run with --fix to cancel')"
exit 1
