#!/bin/bash
# cluster_health_probe.sh — run as root on node4 (the only node that can ssh
# to node1-3). Prints one line per node:
#   nodeN temp=<°C> load1=<1-min load> cpus=<CPUs held by project jobs>
#
# Used by the temperature/load monitor: the hard project policy is <=20 CPUs
# per node; package temperature must stay well under the 88 °C high / 98 °C
# critical thresholds (the 2026-11-15 incident hit 70 °C+ with fans at full
# speed when a run saturated all 64 cores).
BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
# Node list is site-specific: NODES="node4 node1 node2 node3" ./cluster_health_probe.sh
NODES=${NODES:-"node4 node1 node2 node3"}
for n in $NODES; do
  if [ "$n" = "node4" ]; then
    out=$(bash "$BASE/scripts/node_health_probe.sh")
  else
    out=$(ssh -o ConnectTimeout=5 -o BatchMode=yes "root@$n" 'bash -s' < "$BASE/scripts/node_health_probe.sh" 2>/dev/null)
  fi
  echo "$n ${out:-unreachable}"
done
