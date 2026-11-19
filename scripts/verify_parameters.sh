#!/usr/bin/env bash
# verify_parameters.sh — re-derive, from the INSTALLED tools, every parameter
# claim this project makes. Nothing here trusts a web page or a written note:
# each check reads the tool's own artifact and compares it to what we pinned.
#
# Evidence hierarchy used (strongest first):
#   1. the tool's own output        `mixcr exportPreset` — a byte-for-byte
#                                   reproducible expansion of the preset we pass
#   2. the installed source         run-trust4 (perl) and BamExtractor.cpp —
#                                   shows what the tool actually does with flags
#   3. official documentation       mixcr.com preset reference, TRUST4 README
#   4. peer-reviewed literature     DOIs, checked against Crossref metadata
#   5. our own experiments          used to confirm *or overturn* 1-4
#                                   (docs/scientific_rationale.md §7)
#
# Usage: verify_parameters.sh
# Exit: 0 all checks pass, 1 any failure, 2 missing environment.
set -uo pipefail
BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
E=$BASE/.pixi/envs/default/bin
PINNED=$BASE/docs/parameters/mixcr_rna-seq_v4.7.0.yaml
fail=0
ok()   { printf "  PASS  %s\n" "$*"; }
bad()  { printf "  FAIL  %s\n" "$*"; fail=$((fail+1)); }
skip() { printf "  SKIP  %s\n" "$*"; }

[ -x "$E/mixcr" ] || { echo "no MiXCR at $E/mixcr — run on the cluster after pixi install"; exit 2; }
if [ -f "$HOME/.config/tcga-mixcr/env" ]; then . "$HOME/.config/tcga-mixcr/env"; fi

echo "== 1. MiXCR rna-seq preset: what the tool expands our flags to =="
# NB: exportPreset refuses to write over an existing file, so do NOT use mktemp
# (which creates it) — build a name that does not exist yet.
tmp=/tmp/preset_verify_$$_$RANDOM.yaml
rm -f "$tmp"
if "$E/mixcr" exportPreset --preset-name rna-seq -s hsa --assemble-longest-contigs "$tmp" >/dev/null 2>&1; then
  if cmp -s "$tmp" "$PINNED"; then
    ok "exportPreset output is byte-identical to docs/parameters/$(basename "$PINNED") ($(md5sum "$tmp" | cut -c1-12))"
  else
    bad "exportPreset output differs from the pinned file — the pinned copy is stale"
    diff <(head -30 "$tmp") <(head -30 "$PINNED") | head -10 | sed 's/^/        /'
  fi
  # the pipeline the preset actually runs (the YAML key is `pipeline:`, not
  # `steps:` — the preset file is the tool's own export, so trust its shape)
  steps=$(awk '/^pipeline:/{f=1;next} f&&/^[a-z]/{f=0} f&&/^  - /{printf "%s ", $2}' "$tmp")
  echo "        pipeline: $steps"
  for s in align assemblePartial extend assemble assembleContigs exportClones; do
    case " $steps " in *" $s "*) ;; *) bad "step '$s' missing from the preset";; esac
  done
  grep -q 'trimmingQualityThreshold: 10' "$tmp" && ok "align trimmingQualityThreshold=10" || bad "trimmingQualityThreshold not 10"
  grep -q 'repeats: 2'                  "$tmp" && ok "assemblePartial repeats=2"          || bad "assemblePartial repeats not 2"
  grep -q 'clnaOutput: true'            "$tmp" && ok "assemble clnaOutput=true"           || bad "clnaOutput not true"
else
  bad "exportPreset failed (is the license env sourced?)"
fi
rm -f "$tmp"

echo "== 2. TRUST4: what run-trust4 does with our flags (installed source) =="
RT=$E/run-trust4
if [ -f "$RT" ]; then
  grep -q 'bam-extractor -b "\.\$bamFiles\[0\]' "$RT" \
    && ok 'only the FIRST -b reaches bam-extractor ($bamFiles[0], line ~402) — a second -b is silently ignored' \
    || bad 'could not confirm the single -b behaviour in run-trust4'
  grep -qi 'the data is from bulk,non-UMI-based TCR-seq or BCR-seq' "$RT" \
    && ok '--repseq is documented as targeting bulk non-UMI TCR/BCR-seq, so we correctly do not pass it' \
    || bad '--repseq wording not found'
  grep -qi 'the flag in BAM for the unmapped read-pair is nonconcordant' "$RT" \
    && ok '--abnormalUnmapFlag exists with the documented meaning we rely on' \
    || bad '--abnormalUnmapFlag wording not found'
else
  skip "run-trust4 not installed here"
fi
if [ -f "$BASE/docs/parameters/trust4_bam_extractor_notes.txt" ]; then
  ok "bam-extractor behaviour notes pinned (see docs/parameters/trust4_bam_extractor_notes.txt)"
fi

echo "== 3. Pinned reference FASTA digests =="
# PINNED.txt is a human-readable note ("name   md5 <digest>"), not `md5sum -c`
# format, so parse it explicitly.
if [ -f "$BASE/reference/trust4/PINNED.txt" ]; then
  while read -r name _ digest; do
    case "$name" in ''|'#'*) continue;; esac
    [ -n "${digest:-}" ] || continue
    actual=$(md5sum "$BASE/reference/trust4/$name" 2>/dev/null | cut -d' ' -f1)
    if [ "$actual" = "$digest" ]; then
      ok "$name md5 $digest matches the pin"
    else
      bad "$name md5 is ${actual:-missing}, pinned $digest"
    fi
  done < "$BASE/reference/trust4/PINNED.txt"
else
  skip "reference/trust4/PINNED.txt missing"
fi

echo
if [ "$fail" -eq 0 ]; then
  echo "ALL PARAMETER CHECKS PASSED"
  exit 0
fi
echo "$fail check(s) FAILED"
exit 1
