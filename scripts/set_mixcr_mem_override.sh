#!/usr/bin/env bash
# set_mixcr_mem_override.sh — pin MiXCR's self-computed step-JVM heap through the
# ENVIRONMENT, never by editing the launcher script.
#
# THE BUG (2026-11-18)
#   `mixcr analyze` runs each pipeline step as a child JVM whose command line has
#   no -Xmx, so share/mixcr-4.7.0-0/mixcr computes
#       targetXmx = min(defaultXmx=12000, MemFree_MB - 100)
#   MemFree excludes page cache — exactly what concurrent BAM readers create — so
#   on a busy node the step heap collapsed to ~1.95 GB:
#       Not enough memory for run command. This run used approximately 1950m of
#       memory, this machine has 384927 Mb in total
#   The step died and MiXCR restarted the whole pipeline from `align` (6-7 align
#   invocations per sample measured, i.e. two full passes where one is needed).
#
# THE FIX
#   The launcher already honours an environment override BEFORE capping:
#       if [[ -n ${MAX_MEMORY_OVERRIDE+x} ]]; then maxMb="${MAX_MEMORY_OVERRIDE}"; fi
#       targetXmx=defaultXmx; calculateMaxMb
#       if [[ $targetXmx -gt $maxMb ]]; then targetXmx=$maxMb; fi
#   With MAX_MEMORY_OVERRIDE=12000 every launcher-spawned JVM gets exactly
#   -Xmx12000m regardless of MemFree. 12 GB is the launcher's own default and is
#   empirically sufficient — every step that ever completed ran under a <=12 GB
#   child heap; only the top-level `mixcr -Xmx40g analyze` needs 40 GB, and its
#   explicit argument still wins. Deliberately NOT 40 GB: the parent analyze JVM
#   already holds ~27 GB RSS of the 64 GB cgroup, so a 40 GB child risks a
#   cgroup OOM.
#
#   It goes in ~/.config/tcga-mixcr/env because the workflow's `shell_prelude`
#   already sources that file for every rule. That means NO workflow edit — and
#   therefore no rule-definition change, which would invalidate completed
#   instances and re-run every sample the rule had already produced.
#
# WHY NOT PATCH THE LAUNCHER SCRIPT (learned the hard way, same day)
#   Editing `share/mixcr-4.7.0-0/mixcr` while bash is executing it shifts byte
#   offsets; bash re-reads scripts by offset and misparses the tail
#   (`syntax error near unexpected token ']]'`). Under `set -euo pipefail` that
#   killed a rule AFTER MiXCR had printed "Analysis finished successfully." —
#   the science was done, the merge and sentinel never ran, and the sample would
#   have been re-aligned from scratch. Never alter a script that running jobs are
#   executing.
#
# Usage: set_mixcr_mem_override.sh [--check]
set -uo pipefail
ENVF=${ENVF:-$HOME/.config/tcga-mixcr/env}
LINE='export MAX_MEMORY_OVERRIDE=12000'
CHECK=0
[ "${1:-}" = "--check" ] && CHECK=1

[ -f "$ENVF" ] || { echo "missing $ENVF (license env — see AGENTS.md §8.1)"; exit 2; }

if grep -qx "$LINE" "$ENVF"; then
  echo "already set: $LINE ($ENVF)"
  exit 0
fi
if [ "$CHECK" -eq 1 ]; then
  echo "NOT set: $LINE in $ENVF"
  exit 1
fi
{
  echo "# MiXCR step JVMs have no -Xmx of their own and the launcher derives the"
  echo "# heap from MemFree, which page cache drives to ~0 on busy nodes. 12000 is"
  echo "# the launcher's own default ceiling — see scripts/set_mixcr_mem_override.sh."
  echo "$LINE"
} >> "$ENVF"
chmod 600 "$ENVF"
if grep -qx "$LINE" "$ENVF"; then
  echo "set: $LINE -> $ENVF"
  exit 0
fi
echo "FAILED to append to $ENVF"
exit 3
