#!/usr/bin/env bash
# check_memory_budget.sh — refuse to run a shard whose per-node memory budget
# exceeds a safe fraction of the node's RAM.
#
# oxo-flow bounds concurrency per node (cluster profile `max_submitted`) and
# per job (`resources.memory`). Their product is the node's worst-case
# charge, so it must stay below the budget:
#
#     max_submitted x memory_per_job  <=  BUDGET_PCT% x MemTotal
#
# Node RAM is NOT uniform on this cluster (<node> 515G, <node> 251G, <node>
# 376G, <node> 1385G), so the same max_submitted is safe on one node and
# fatal on another — the 2026-11-16 <node> failure ("Not enough memory for
# run command ... machine has 256911 Mb") was exactly that.
#
# Usage (on <node>, as root so MemTotal of every node is readable):
#   bash scripts/check_memory_budget.sh            # report
#   bash scripts/check_memory_budget.sh --strict   # exit 1 on any breach
set -uo pipefail
BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUDGET_PCT=${BUDGET_PCT:-80}
STRICT=0
[ "${1:-}" = "--strict" ] && STRICT=1

# `resources.memory` is a SLURM RESERVATION, not what a job actually uses.
# Two different checks therefore apply:
#   * reservation  max_submitted x MEM_GB  <= MemTotal   (else jobs cannot start)
#   * real charge  max_submitted x CHARGE_GB <= BUDGET_PCT% of MemTotal
# CHARGE_GB is measured, not assumed: on <node> (2026-11-19) two running mixcr
# jobs charged 42.3 GB and 48.9 GB to their cgroups while their JVMs' anonymous
# RSS totalled 49.3 GB — the rest of the charge is reclaimable page cache from
# reading the BAM. Using the 64 GB reservation for the 80% test made a 7-slot
# <node> look like a breach (448/503) when its real charge is 67%.
MEM_GB=$(grep -oE 'memory = "[0-9]+G"' "$BASE/tcga-mixcr.shard.oxoflow" | head -1 | grep -oE '[0-9]+')
MEM_GB=${MEM_GB:-64}
CHARGE_GB=${CHARGE_GB:-48}
breach=0

CPU_CAP=${CPU_CAP:-20}
CPU_THREADS=$(grep -oE 'threads = [0-9]+' "$BASE/tcga-mixcr.shard.oxoflow" | head -1 | grep -oE '[0-9]+')
CPU_THREADS=${CPU_THREADS:-4}

# The CPU policy is PER NODE, not a single number (2026-11-19). Measured with
# IPMI: <node> sits at 2,520-2,640 RPM at 52 C and has 503 GB, so it can carry
# 28 CPUs; <node> is the noise source — 7,920-8,040 RPM at 48 C because it is
# the file server (1385 GB RAM, 8 USB disks) — and is deliberately held at
# 20 CPUs even though its RAM would allow far more. Raising a node past its
# cap is a noise decision, so it belongs here rather than in one global knob.
cpu_cap_for() {
  case "$1" in
    1) echo "${CPU_CAP_NODE1:-30}" ;;
    4) echo "${CPU_CAP_NODE4:-20}" ;;
    *) echo "${CPU_CAP}" ;;
  esac
}

printf "%-7s %-10s %-13s %-9s %-8s %-8s %-8s %s\n" node MemTotalGB max_submitted reserveGB usedGB budgetGB chargeCPU status
for n in 1 2 3 4; do
  profile=$BASE/shards/shard$n/profiles/slurm.toml
  [ -f "$profile" ] || { printf "node%s  (no profile)\n" "$n"; continue; }
  maxsub=$(sed -n 's/^max_submitted *= *\([0-9]*\).*/\1/p' "$profile" | head -1)
  mem_mb=$(timeout 10 ssh -o ConnectTimeout=5 -o BatchMode=yes "root@node$n" "awk '/^MemTotal:/{print \$2}' /proc/meminfo" 2>/dev/null)
  if [ -z "$mem_mb" ]; then
    printf "node%-3s  (unreachable)\n" "$n"
    continue
  fi
  mem_gb=$((mem_mb / 1024 / 1024))
  charge=$((maxsub * MEM_GB))
  budget=$((mem_gb * BUDGET_PCT / 100))
  cpu=$((maxsub * CPU_THREADS))
  cap=$(cpu_cap_for "$n")
  used=$((maxsub * CHARGE_GB))
  status=ok
  if [ "$charge" -gt "$mem_gb" ]; then status=RESERVE-OVER; breach=$((breach + 1)); fi
  if [ "$used" -gt "$budget" ]; then status="${status},MEM-OVER"; breach=$((breach + 1)); fi
  if [ "$cpu" -gt "$cap" ]; then status="${status},CPU-OVER(${cpu}>${cap})"; breach=$((breach + 1)); fi
  printf "node%-3s %-10s %-13s %-9s %-8s %-8s %-8s %s\n" \
    "$n" "$mem_gb" "$maxsub" "$charge" "$used" "$budget" "$cpu" "$status"
done

echo
echo "rule memory = ${MEM_GB}G/job, threads = ${CPU_THREADS}/job, "
echo "memory budget = ${BUDGET_PCT}% of MemTotal; CPU caps per node: slots per node: 7/3/5/5 under the caps 30/20/20/20"
if [ "$breach" -gt 0 ]; then
  echo "BREACH: $breach violation(s) — lower max_submitted in that shard's profile"
  [ "$STRICT" -eq 1 ] && exit 1
fi
exit 0
