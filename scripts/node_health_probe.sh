#!/bin/bash
# Payload executed ON each node (fed over ssh stdin by cluster_health_probe.sh).
# Prints one line: temp=<pkg °C> load1=<1-min> cpus=<CPUs allocated to project jobs>
t=$(sensors 2>/dev/null | awk '/Package id 0/{gsub(/[+°C]/,"",$4); print $4; exit}')
[ -z "$t" ] && t=$(sensors 2>/dev/null | awk '/^Core 0/{gsub(/[+°C]/,"",$4); print $4; exit}')
l=$(cut -d' ' -f1 /proc/loadavg)
c=$(squeue -h -u wsx -w "$(hostname -s)" -o "%C" 2>/dev/null | paste -sd+ - | bc 2>/dev/null)
echo "temp=${t:-na} load1=$l cpus=${c:-0}"
