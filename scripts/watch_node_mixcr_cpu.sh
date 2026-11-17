#!/usr/bin/env bash
# watch_node_mixcr_cpu.sh — detect a wedged MiXCR job by CPU-time progress.
#
# File-write heuristics give false positives: during alignment MiXCR buffers in
# memory and writes its output at the end, so tmp/ can be quiet for an hour
# while the job is healthy (observed 2026-11-17). CPU time is the honest
# signal — a running alignment burns seconds of CPU per minute.
#
# Matches pids ACROSS the two snapshots (a pid that exited must not be
# compared against a different process — the first version summed by position
# and reported a negative delta).
#
# Run ON the node (root or the job owner):
#   bash watch_node_mixcr_cpu.sh [sampling_seconds] [min_ticks]
# Exit 0 = advancing, 3 = stalled, 4 = no MiXCR process (nothing to judge).
set -uo pipefail
SAMPLE=${1:-45}
MIN_TICKS=${2:-30}       # 30 ticks = 0.3 s CPU over the sample window

pids=$(pgrep -f 'mixcr\.jar' 2>/dev/null | head -20)
[ -z "$pids" ] && { echo "idle (no mixcr jar process)"; exit 4; }

declare -A before
for p in $pids; do
  [ -r "/proc/$p/stat" ] || continue
  before[$p]=$(awk '{print $14 + $15}' "/proc/$p/stat")
done

sleep "$SAMPLE"

total=0
alive=0
for p in "${!before[@]}"; do
  [ -r "/proc/$p/stat" ] || continue
  now=$(awk '{print $14 + $15}' "/proc/$p/stat")
  delta=$(( now - ${before[$p]} ))
  [ "$delta" -gt 0 ] && total=$(( total + delta ))
  alive=$(( alive + 1 ))
done

echo "alive=$alive cpu_ticks_in_${SAMPLE}s=$total min=$MIN_TICKS"
if [ "$alive" -eq 0 ]; then echo "no surviving mixcr process"; exit 4; fi
if [ "$total" -lt "$MIN_TICKS" ]; then echo "STALLED"; exit 3; fi
echo "advancing"
exit 0
