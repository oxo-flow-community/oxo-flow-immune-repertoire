#!/usr/bin/env bash
# probe_mixcr_memory.sh — two checks on a live node:
#
#   A. does a small alignment run at all?  (node-level pressure check)
#   B. does an explicit `-Xmx` *argument* beat JAVA_TOOL_OPTIONS?  The MiXCR
#      launcher appends its own -Xmx (min(12 GB, free-memory cap)) when no
#      -Xmx argument is present, and that later JVM option overrides the env
#      var — so JAVA_TOOL_OPTIONS=-Xmx40g never took effect (verified
#      2026-11-17: "This run used approximately 1394m" is Runtime.maxMemory).
#      With JAVA_TOOL_OPTIONS=-Xmx16m (far too small to boot MiXCR) and an
#      explicit -Xmx2g argument, a successful run proves the argument wins.
#
# Uses the shard's pixi env directly (no pixi binary on PATH under srun):
#   srun -n 1 -p workq -w <node> -c 1 --mem=4G bash \
#        $HOME/Projects/TCGA-mixcr/scripts/probe_mixcr_memory.sh
set -uo pipefail
S=${SHARD_DIR:-$HOME/Projects/TCGA-mixcr/shards/shard3}
cd "$S" || exit 1
PATH="$S/.pixi/envs/default/bin:$PATH"
command -v mixcr >/dev/null || { echo "mixcr not found in $S/.pixi/envs/default/bin"; exit 3; }

# BAM is site-specific: pass it in (BAM=/path/to/sample.bam) rather than
# relying on the hardcoded pattern this project uses.
BAM=${BAM:-$(ls ${BAM_GLOB:-*.bam} 2>/dev/null | head -1)}
[ -n "$BAM" ] || { echo "set BAM=/path/to/sample.bam (or BAM_GLOB) first"; exit 3; }
echo "host=$(hostname -s) bam=$BAM"
grep -E "MemFree|MemAvailable" /proc/meminfo

mkdir -p tmp/probe
samtools view -h "$BAM" 2>/dev/null | head -4000 | samtools view -b -o tmp/probe/p.bam - 2>/dev/null
ENVF=$HOME/.config/tcga-mixcr/env
[ -f "$ENVF" ] || ENVF=$HOME/.config/tcga-mixcr/env
# shellcheck disable=SC1090
. "$ENVF" 2>/dev/null || echo "license env not found ($ENVF)"

echo "--- A: plain alignment"
mixcr -Xmx2g align -f --preset rna-seq --species hsa -t 1 \
  tmp/probe/p.bam tmp/probe/a.vdjca > tmp/probe/a.log 2>&1
echo "A-RC=$? $(tail -1 tmp/probe/a.log)"

echo "--- B: -Xmx as argument vs JAVA_TOOL_OPTIONS (env set far too small)"
JAVA_TOOL_OPTIONS="-Xmx16m" mixcr -Xmx2g align -f --preset rna-seq --species hsa -t 1 \
  tmp/probe/p.bam tmp/probe/b.vdjca > tmp/probe/b.log 2>&1
rc=$?
echo "B-RC=$rc $(tail -1 tmp/probe/b.log)"
if [ "$rc" -eq 0 ]; then
  echo "VERDICT: explicit -Xmx argument wins over JAVA_TOOL_OPTIONS"
else
  echo "VERDICT: JAVA_TOOL_OPTIONS still overrides — the workflow fix is not enough"
fi
rm -rf tmp/probe
exit "$rc"
