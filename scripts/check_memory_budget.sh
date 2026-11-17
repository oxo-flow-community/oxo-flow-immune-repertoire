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
# Node RAM is NOT uniform on this cluster (e.g. 515 GB, 251 GB,
# 376 GB, 1385 GB), so the same max_submitted is safe on one node and
# fatal on another — the 2026-11-16 node2 failure ("Not enough memory for
# run command ... machine has 256911 Mb") was exactly that.
#
# Usage (on node4, as root so MemTotal of every node is readable):
#   bash scripts/check_memory_budget.sh            # report
#   bash scripts/check_memory_budget.sh --strict   # exit 1 on any breach
set -uo pipefail
BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUDGET_PCT=${BUDGET_PCT:-80}
STRICT=0
[ "${1:-}" = "--strict" ] && STRICT=1

# The largest per-rule memory in the workflow decides the per-job charge.
MEM_GB=$(grep -oE 'memory = "[0-9]+G"' "$BASE/tcga-mixcr.shard.oxoflow" | head -1 | grep -oE '[0-9]+')
MEM_GB=${MEM_GB:-64}
breach=0

CPU_CAP=${CPU_CAP:-20}
CPU_THREADS=$(grep -oE 'threads = [0-9]+' "$BASE/tcga-mixcr.shard.oxoflow" | head -1 | grep -oE '[0-9]+')
CPU_THREADS=${CPU_THREADS:-4}

printf "%-7s %-10s %-13s %-9s %-8s %-8s %s\n" node MemTotalGB max_submitted chargeGB budgetGB chargeCPU status
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
  status=ok
  if [ "$charge" -gt "$budget" ]; then status=MEM-OVER; breach=$((breach + 1)); fi
  if [ "$cpu" -gt "$CPU_CAP" ]; then status="${status},CPU-OVER"; breach=$((breach + 1)); fi
  printf "node%-3s %-10s %-13s %-9s %-8s %-8s %s\n" "$n" "$mem_gb" "$maxsub" "$charge" "$budget" "$cpu" "$status"
done

echo
echo "rule memory = ${MEM_GB}G/job, threads = ${CPU_THREADS}/job, "
echo "memory budget = ${BUDGET_PCT}% of MemTotal, CPU cap = ${CPU_CAP}/node (project policy)"
if [ "$breach" -gt 0 ]; then
  echo "BREACH: $breach violation(s) — lower max_submitted in that shard's profile"
  [ "$STRICT" -eq 1 ] && exit 1
fi
exit 0
