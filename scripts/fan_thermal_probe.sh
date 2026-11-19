#!/usr/bin/env bash
# fan_thermal_probe.sh — per-node CPU temperature, FAN SPEED and our-job CPU use.
#
# Why: the <=20 CPUs/node policy exists because the machine room complained
# about fan noise (2026-11-15: 4 nodes at load ~50, package 70 C, fans at full
# speed). Cooling can only be judged on the two numbers that matter — package
# temperature and fan RPM — and `sensors` exposes no fans on these hosts, so we
# read them over IPMI (`ipmitool`, available on <node>-4 without a BMC password
# for local/root use; `sensors` covers the temperatures).
#
# It also measures how much of the allowed CPU the jobs actually consume, so a
# slot that is idle-but-reserved is visible instead of implied.
#
# Usage: fan_thermal_probe.sh [-q]
#   (default) human table per node, plus a one-line summary
#   -q        a single machine-readable line per run (for logging/trending)
# Exit: 0 always (a probe must not break a monitoring loop); the caller decides.
set -uo pipefail
QUIET=0
[ "${1:-}" = "-q" ] && QUIET=1

# Cluster-agnostic: the node list and which one is local come from the
# environment, so the published copy carries no site names.
NODES=${NODES:?set your node names, e.g. NODES="n1 n2 n3"}
LOCAL_NODE=${LOCAL_NODE:-$(hostname -s)}

# The node that runs this script is probed locally; the others are reached by
# ssh (they must trust root for the ipmitool call).
probe_local() {
  temp=$(sensors 2>/dev/null | awk '/Package id 0/{gsub(/[+°C]/,"",$4); print $4; exit}')
  fans=$(ipmitool sensor 2>/dev/null | awk -F'|' '/FAN[0-9]+ Speed/{gsub(/ /,"",$2); printf "%s ", $2}')
  fanmax=$(echo "$fans" | tr ' ' '\n' | sort -rn | head -1)
  load=$(awk '{print $1}' /proc/loadavg)
  cores=$(nproc)
  ourcpus=$(ps -eo pcpu,args --no-headers 2>/dev/null | grep -E "[j]ava .*-Xmx|[b]am-extractor|[s]amtools view -@" \
            | awk '{s+=$1} END{printf "%.0f", s}')
  echo "$temp|$fanmax|$load|$cores|${ourcpus:-0}"
}
probe_remote() {
  ssh -o ConnectTimeout=8 "root@$1" '
    temp=$(sensors 2>/dev/null | awk "/Package id 0/{gsub(/[+°C]/,\"\",\$4); print \$4; exit}")
    fan=$(ipmitool sensor 2>/dev/null | awk -F"|" "/FAN[0-9]+ Speed/{gsub(/ /,\"\",\$2); print \$2}" | sort -rn | head -1)
    load=$(awk "{print \$1}" /proc/loadavg)
    cores=$(nproc)
    ourcpus=$(ps -eo pcpu,args --no-headers 2>/dev/null | grep -E "[j]ava .*-Xmx|[b]am-extractor|[s]amtools view -@" \
              | awk "{s+=\$1} END{printf \"%.0f\", s}")
    echo "$temp|$fan|$load|$cores|${ourcpus:-0}"' 2>/dev/null
}

if [ "$QUIET" -eq 0 ]; then
  printf "%-7s %-7s %-9s %-8s %-10s %s\n" NODE TEMP fan_RPM load cores "our_CPU%"
fi
for n in $NODES; do
  if [ "$n" = "$LOCAL_NODE" ]; then out=$(probe_local); else out=$(probe_remote "$n"); fi
  if [ -z "$out" ]; then
    [ "$QUIET" -eq 0 ] && printf "%-7s %s\n" "$n" "unreachable"
    continue
  fi
  IFS='|' read -r temp fan load cores ourcpus <<< "$out"
  if [ "$QUIET" -eq 1 ]; then
    echo "$n temp=${temp:-NA} fan=${fan:-NA} load=${load:-NA} cores=${cores:-NA} ourcpu=${ourcpus:-NA}"
  else
    printf "%-7s %-7s %-9s %-8s %-10s %s\n" "$n" "${temp:-?}C" "${fan:-?}" "$load" "$cores" "$ourcpus"
  fi
done
exit 0
