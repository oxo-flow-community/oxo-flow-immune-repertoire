#!/usr/bin/env bash
# backfill_mixcr_merge.sh — finish a MiXCR rule whose `mixcr analyze` completed
# but whose post-processing did not run.
#
# WHY THIS EXISTS (2026-11-18 incident)
#   The rule shell runs `mixcr ... | tee -a log` under `set -euo pipefail`.
#   MiXCR itself printed "Analysis finished successfully." and wrote all six
#   per-chain `*.clones_*.tsv` tables, but the launcher script then died with
#   `syntax error near unexpected token ']]'` because the launcher file had been
#   rewritten while bash was still executing it (bash re-reads scripts by byte
#   offset). `pipefail` turned that into a nonzero pipeline status, `set -e`
#   aborted the rule before the merge and the sentinel, and the engine recorded
#   FAILED. Re-running the rule costs ~10 h of alignment for work that already
#   succeeded, so instead: redo exactly the two steps that were skipped.
#
# WHAT IT DOES — byte-for-byte the same operations as the rule's shell:
#   1. merge `*.clones_*.tsv` into `<barcode>.mixcr.clones.tsv` with a trailing
#      `chain` column (an empty chain suffix becomes `ambiguous`)
#   2. write the completion sentinel `<barcode>.mixcr.done` last
# It refuses to act unless measurable evidence says the analyze step finished:
# at least one per-chain table, a non-empty qc.json, and a readinfo.tsv.
#
# Usage: backfill_mixcr_merge.sh <sample-output-dir> [--apply]
#        (dry-run by default; prints every check it makes)
set -uo pipefail
DIR=${1:?usage: backfill_mixcr_merge.sh <sample-output-dir> [--apply]}
APPLY=0
[ "${2:-}" = "--apply" ] && APPLY=1
[ -d "$DIR" ] || { echo "not a directory: $DIR"; exit 2; }

shopt -s nullglob
parts=( "$DIR"/*.clones_*.tsv )
# find -print -quit, not `ls` under nullglob: with no match `ls` lists the CWD
# and the NOTE below then describes an unrelated file (2026-11-28 audit L1).
clones=$(find "$DIR" -maxdepth 1 -name '*.mixcr.clones.tsv' -print -quit 2>/dev/null)
qc=$(find "$DIR" -maxdepth 1 -name '*.mixcr.qc.json' -print -quit 2>/dev/null)
barcode=""
[ -n "$qc" ] && barcode=$(basename "$qc" .mixcr.qc.json)
[ -n "$barcode" ] || { echo "REFUSE: no *.mixcr.qc.json — analyze never produced its report"; exit 3; }
[ "${#parts[@]}" -gt 0 ] || { echo "REFUSE: no *.clones_*.tsv — zero-clonotype sample, use the rule itself"; exit 3; }
[ -s "$DIR/$barcode.mixcr.readinfo.tsv" ] || { echo "REFUSE: readinfo empty/missing"; exit 3; }
# The sentinel this script writes claims "genomic-only". Verify that from the
# align report's inputFiles instead of trusting the call site: a pre-cutover
# sample would have aligned a samtools-merged BAM and must not be relabelled.
report=$DIR/$barcode.mixcr.align.report.json
bam_input=$(grep -o '"[^"]*\.bam"' "$report" 2>/dev/null | head -1)
case "$bam_input" in
  "") echo "REFUSE: cannot verify the align input from $report"; exit 6 ;;
  *merged*|*chimeric*) echo "REFUSE: align input is not genomic-only ($bam_input)"; exit 6 ;;
esac
echo "sample      : $barcode"
echo "align input : $bam_input (genomic-only)"
echo "chain tables: ${#parts[@]} -> $(for p in "${parts[@]}"; do printf '%s ' "${p##*.clones_}"; done)"
if [ -n "$clones" ]; then
  echo "NOTE: $clones already exists ($(wc -l < "$clones") lines) — it will be rewritten from the per-chain tables"
fi
want_rows=$(for p in "${parts[@]}"; do echo $(( $(wc -l < "$p") - 1 )); done | awk '{s+=$1} END{print s}')
echo "expected data rows: $want_rows"
[ "$APPLY" -eq 1 ] || { echo "DRY-RUN: no files written. Re-run with --apply."; exit 0; }

PREFIX=$DIR/$barcode.mixcr
tmpf=$PREFIX.clones.tsv.tmp.$$
{ head -1 "${parts[0]}" | sed 's/$/\tchain/'
  for f in "${parts[@]}"; do
    ch=$(basename "$f" .tsv); ch=${ch##*.clones_}
    [ -z "$ch" ] && ch="ambiguous"
    tail -n +2 "$f" | sed "s/$/\t$ch/"
  done; } > "$tmpf"
got_rows=$(( $(wc -l < "$tmpf") - 1 ))
if [ "$got_rows" -ne "$want_rows" ]; then
  echo "FAIL: merged $got_rows rows, expected $want_rows — leaving the original in place"
  rm -f "$tmpf"; exit 4
fi
mv "$tmpf" "$PREFIX.clones.tsv"
echo "WROTE $PREFIX.clones.tsv ($got_rows rows, $(head -1 "$PREFIX.clones.tsv" | awk -F'\t' '{print NF}') cols)"

# verification before the sentinel: chain values legal, per-chain cloneFraction ~1
awk -F'\t' -v bc="$barcode" '
  NR==1 { for (i=1;i<=NF;i++) {
            if ($i=="chain") ci=i
            # MiXCR 4.x exports readFraction; 3.x called it cloneFraction
            if ($i=="readFraction" || $i=="cloneFraction") fi=i
          }
          if (!ci) { print "FAIL: no chain column"; exit 9 }
          if (!fi) { print "FAIL: no readFraction/cloneFraction column"; exit 9 }
          hdr=1; next }
  hdr { n++; if ($ci=="") { e++ } c[$ci]++; s[$ci]+=($fi+0) }
  END {
    if (!hdr) { print "FAIL: empty table"; exit 9 }
    if (e) { printf "FAIL: %d rows with empty chain\n", e; exit 9 }
    bad=0
    for (k in c) { d = s[k]-1; if (d<0) d=-d; if (d>0.001) { printf "FAIL: chain %s fraction sum=%.6f\n", k, s[k]; bad=1 } }
    if (bad) exit 9
    printf "OK: %d rows, chains=%d, every chain readFraction sums to 1\n", n, length(c)
  }' "$PREFIX.clones.tsv" || { echo "verification failed — sentinel NOT written"; exit 5; }

echo "mixcr-v2 genomic-only threads=4 xmx=40g $(date -Is) [backfilled: merge after launcher fault]" > "$PREFIX.done"
echo "WROTE $PREFIX.done"
