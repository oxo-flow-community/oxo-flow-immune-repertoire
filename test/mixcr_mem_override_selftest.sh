#!/usr/bin/env bash
# mixcr_mem_override_selftest.sh — prove that MAX_MEMORY_OVERRIDE pins the MiXCR
# step-JVM heap independently of MemFree, and that the env file carries it.
#
# Mechanism (share/mixcr-4.7.0-0/mixcr):
#     targetXmx = min(defaultXmx=12000, MemFree_MB - 100)
#     if MAX_MEMORY_OVERRIDE is set, maxMb := MAX_MEMORY_OVERRIDE
# Page cache is excluded from MemFree, so N concurrent BAM readers drive the
# computed heap toward zero; live on 2026-11-18 a step JVM got ~1950m on a
# 384 GB node and MiXCR restarted the pipeline from `align`.
#
# The test shadows `free` with a fake reporting 2 GB free, then reads the JVM's
# real MaxHeapSize through -XX:+PrintFlagsFinal. No MiXCR run, no license, no I/O.
#
# Usage: bash test/mixcr_mem_override_selftest.sh
set -uo pipefail
BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
L=$(ls -d "$BASE"/shards/shard*/.pixi/envs/*/share/mixcr-*/mixcr 2>/dev/null | head -1)
[ -n "$L" ] || L=$(ls -d "$BASE"/.pixi/envs/*/share/mixcr-*/mixcr 2>/dev/null | head -1)
if [ -z "$L" ] || [ ! -f "$L" ]; then
  echo "SKIP: no MiXCR launcher found (run on the cluster after pixi install)"
  exit 0
fi

FAKE=$(mktemp -d)
trap 'rm -rf "$FAKE"' EXIT
cat > "$FAKE/free" <<'FAKE_FREE'
#!/bin/sh
# MemFree = 2000 MB: a node whose page cache has eaten everything
echo "              total        used        free      shared  buff/cache   available"
echo "Mem:         384927      300000        2000         244      179280      256682"
FAKE_FREE
chmod +x "$FAKE/free"

heap_mb() {  # $1 = extra env assignment ("" for none) -> MaxHeapSize in MB
  env $1 PATH="$FAKE:$PATH" "$L" -XX:+PrintFlagsFinal --help 2>&1 \
    | awk '/MaxHeapSize/{printf "%d", $4/1048576; exit}'
}

before=$(heap_mb "")
after=$(heap_mb "MAX_MEMORY_OVERRIDE=12000")
envfile=${ENVF:-$HOME/.config/tcga-mixcr/env}
infile=$(grep -qx 'export MAX_MEMORY_OVERRIDE=12000' "$envfile" 2>/dev/null && echo yes || echo no)

echo "MemFree simulated        : 2000 MB"
echo "heap without override    : ${before} MB   (expected ~1900: the collapse)"
echo "heap with override       : ${after} MB   (expected 12000: pinned)"
echo "env file exports it      : $infile ($envfile)"

rc=0
[ "$after" -eq 12000 ] || { echo "FAIL: override does not pin the heap"; rc=1; }
[ "$before" -lt 3000 ] || { echo "FAIL: could not reproduce the collapse"; rc=1; }
[ "$infile" = yes ] || { echo "FAIL: $envfile does not export MAX_MEMORY_OVERRIDE (run scripts/set_mixcr_mem_override.sh)"; rc=1; }
[ "$rc" -eq 0 ] && echo "PASS: step-JVM heap is pinned and the env file carries the override"
exit "$rc"
