#!/usr/bin/env bash
# driver_shard.sh — production driver for one TCGA-mixcr shard using the
# oxo-flow Slurm **cluster backend** (replaces supervisor.sh / sbatch wrapper).
#
# Usage:  driver_shard.sh <shardN>          (N = 1..4)
# Launch: setsid nohup driver_shard.sh 1 >/dev/null 2>&1 &
#
# Scheduling model — CHUNKED sample selection.
#   oxo-flow's `priority` is a *global* rule ordering (documented as
#   "effective = declared + rounds waited" aging): a single run over 2,780
#   samples is phase-ordered — every mixcr, then every trust4 — and the first
#   `unified.tsv` only appears once both producer phases are done. Giving the
#   downstream rules a higher priority does not fix it, it inverts it: every
#   trust4 then outranks every mixcr and the long pole starves (observed
#   2026-11-16: 20/20 slots on trust4, no mixcr for two hours).
#
#   The engine's own mechanism for batch-complete execution is sample
#   selection: `--samples @sheet.tsv` replaces the sample set of that run.
#   This driver therefore walks the shard in chunks sized by CHUNK_SIZES and
#   every chunk runs end-to-end (mixcr -> trust4 -> standardize -> cleanup)
#   before the next one starts. That yields results per batch, keeps planning
#   proportional to the chunk (seconds, not the ~12 min an 11k-instance plan
#   costs) and confines a bad sample to its own chunk.
#
# Behaviour:
#   * one instance per shard (flock); duplicate launches exit immediately
#   * chunks/chunk_NNN.tsv generated from manifests/groups.tsv, regenerated
#     only when the sample list or the chunk plan changes (md5 stamp)
#   * CHUNK_SIZES is a staircase, not a constant: chunk i uses the i-th size
#     (the last size repeats). A small first chunk yields a fully validated
#     end-to-end result within ~1 day; later chunks are as large as the
#     straight-line-chunk tradeoff allows. Rationale: inside a chunk the
#     engine still prefers every ready `mixcr` over any `trust4` (rule
#     priority), so a chunk only produces unified output once its producer
#     phase has drained — first-result latency is ~chunk_size/slots waves.
#   * per chunk: `run --profile slurm --samples @chunk`, up to CHUNK_RETRIES
#     retries with backoff; a chunk that keeps failing is logged and the
#     driver moves on — one bad sample never blocks the shard
#   * after a full pass the driver compares the per-sample output count with
#     the shard's sample count and repeats (up to MAX_CYCLES) for stragglers
#   * before every attempt: cancel leftover jobs from a previous generation
#     (matched by Slurm WorkDir) and sweep stale per-sample scratch
#   * all output appended to logs/driver_shardN.log
#
# CPU policy: per-rule resources are 4 threads, the profile caps
# max_submitted at 5 => 20 CPUs/node (hard policy: <= 20).
set -uo pipefail

SHARD=${1:?usage: driver_shard.sh <1|2|3|4>}
BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
DIR=$BASE/shards/shard$SHARD
LOG=$BASE/logs/driver_shard$SHARD.log
LOCK=$BASE/logs/.driver_shard$SHARD.lock
OXO=${OXO:-$(command -v oxo-flow || echo $HOME/.cargo/bin/oxo-flow)}
CHUNK_SIZES=${CHUNK_SIZES:-"10 40 100"}
CHUNK_RETRIES=${CHUNK_RETRIES:-2}
MAX_CYCLES=${MAX_CYCLES:-3}
STALE_MIN=${STALE_MIN:-360}
# Slurm account owning the jobs: the owner of the shard directory, NOT $USER —
# the drivers are sometimes started from a root shell on node4 (bioinfo-node4),
# and `squeue -u root` sees none of wsx's jobs. That bug let a restart leave
# the previous generation running while the new one submitted on top of it
# (2026-11-17).
JOB_USER=${JOB_USER:-$(stat -c %U "$DIR" 2>/dev/null || echo "$USER")}

[ -d "$DIR" ] || { echo "missing shard dir: $DIR" >&2; exit 1; }
[ -x "$OXO" ] || { echo "missing oxo-flow binary: $OXO" >&2; exit 1; }

# Guard: the shard must run the workflow that is in the repository. On
# 2026-11-17 shards 2-4 silently kept an OLD workflow after a fix (the cutover
# had not been re-run), so the "fixed" re-run still merged the chimeric BAM on
# three of four shards. Refuse to start on a mismatch.
REPO_WF=$BASE/tcga-mixcr.shard.oxoflow
SHARD_WF=$DIR/tcga-mixcr.oxoflow
if [ -f "$REPO_WF" ] && [ -f "$SHARD_WF" ]; then
  r=$(md5sum "$REPO_WF" | cut -d' ' -f1)
  s_=$(md5sum "$SHARD_WF" | cut -d' ' -f1)
  if [ "$r" != "$s_" ]; then
    echo "FATAL: $SHARD_WF (${s_}) differs from $REPO_WF (${r}) — run scripts/cutover_cluster_backend.sh --apply first" >&2
    echo "$(date -Is) FATAL workflow mismatch shard=${s_} repo=${r}" >> "$LOG"
    exit 5
  fi
fi

exec 9>"$LOCK"
flock -n 9 || { echo "driver already running for shard$SHARD" >&2; exit 0; }

cd "$DIR" || exit 1
export PATH="$HOME/.cargo/bin:$PATH"

EXPECTED=$(awk -F'\t' 'NR>1 { s += split($2, a, ",") } END { print s + 0 }' manifests/groups.tsv)

# Cancel EVERY live rule job that belongs to this shard, not merely the ids
# recorded in the newest run directory: a driver restart (or a workflow
# reconfiguration) leaves the previous generation's Slurm jobs running, and
# uncancelled they double the node's CPU footprint (incident 2026-11-16:
# 32 CPUs/node against the 20 policy). The job's Slurm WorkDir is the shard
# directory, which disambiguates concurrently running shards.
cancel_leftovers() {
  local id wd n=0
  for id in $(squeue -u "$JOB_USER" -h -o "%i %j" 2>/dev/null \
      | awk '$2 ~ /^(mixcr|trust4|standardize|cleanup)/ {print $1}'); do
    wd=$(scontrol show job "$id" 2>/dev/null | sed -n 's/.*WorkDir=\([^ ]*\).*/\1/p' | head -1)
    [ "$wd" = "$PWD" ] || continue
    scancel "$id" 2>/dev/null && n=$((n + 1))
  done
  [ "$n" -gt 0 ] && echo "$(date -Is) cancelled $n leftover rule job(s) from a previous generation" >> "$LOG"
  for _ in $(seq 1 30); do
    [ -z "$(squeue -u "$JOB_USER" -h -o "%j" 2>/dev/null | grep -E '^(mixcr|trust4|standardize|cleanup)')" ] && break
    sleep 2
  done
}

# Reclaim per-sample scratch no live rule can own: a running rule keeps
# writing inside its own tmp/<uuid>/<rule>/ directory, so a whole subtree
# untouched for STALE_MIN minutes is orphaned (jobs killed before their own
# rm -rf, or debris from the pre-cutover layout).
sweep_stale_tmp() {
  local d n=0
  [ -d tmp ] || return 0
  for d in tmp/*/; do
    [ -d "$d" ] || continue
    if [ -n "$(find "$d" -type f -mmin -"$STALE_MIN" 2>/dev/null | head -1)" ]; then
      continue
    fi
    rm -rf "$d" && n=$((n + 1))
  done
  [ "$n" -gt 0 ] && echo "$(date -Is) swept $n stale tmp director(ies) older than ${STALE_MIN}m" >> "$LOG"
}

# Deterministic chunks from the shard's sample list; regenerated only when
# groups.tsv or the chunk plan changes (stamp covers both — a plan change
# with an unchanged sample list previously left stale chunks behind).
prepare_chunks() {
  local stamp
  stamp="$(md5sum manifests/groups.tsv | cut -d' ' -f1) sizes=$CHUNK_SIZES"
  if [ -f chunks/.stamp ] && [ "$(cat chunks/.stamp)" = "$stamp" ]; then
    return 0
  fi
  rm -rf chunks
  mkdir -p chunks
  awk -F'\t' -v sizes="$CHUNK_SIZES" '
    BEGIN { ns = split(sizes, sz, " ") }
    NR > 1 {
      n = split($2, a, ",")
      for (i = 1; i <= n; i++) {
        if (left <= 0) {
          k++
          ci = (k <= ns) ? k : ns
          left = sz[ci]
        }
        buf[k] = buf[k] (buf[k] == "" ? "" : ",") a[i]
        left--
      }
    }
    END {
      for (j = 1; j <= k; j++)
        printf "name\tsamples\nall\t%s\n", buf[j] > sprintf("chunks/chunk_%03d.tsv", j)
    }' manifests/groups.tsv
  echo "$stamp" > chunks/.stamp
  echo "$(date -Is) prepared $(ls chunks/chunk_*.tsv | wc -l) chunk(s) with sizes '$CHUNK_SIZES'" >> "$LOG"
}

# One pass over every chunk. A failing chunk is retried; after CHUNK_RETRIES
# it is logged and skipped so the rest of the shard proceeds.
run_cycle() {
  local cycle=$1 chunk rc tries
  for chunk in chunks/chunk_*.tsv; do
    [ -s "$chunk" ] || continue
    tries=0
    while :; do
      echo "$(date -Is) cycle $cycle chunk $(basename "$chunk") start (try $((tries + 1)))" >> "$LOG"
      "$OXO" run tcga-mixcr.oxoflow --profile slurm --samples "@$chunk" >> "$LOG" 2>&1
      rc=$?
      [ "$rc" -eq 0 ] && break
      tries=$((tries + 1))
      if [ "$tries" -gt "$CHUNK_RETRIES" ]; then
        echo "$(date -Is) chunk $(basename "$chunk") FAILED rc=$rc after $tries retries — continuing with the next chunk" >> "$LOG"
        break
      fi
      cancel_leftovers
      sweep_stale_tmp
      sleep $((60 * tries))
    done
  done
}

prepare_chunks
cycle=0
while [ "$cycle" -lt "$MAX_CYCLES" ]; do
  cycle=$((cycle + 1))
  cancel_leftovers
  sweep_stale_tmp
  echo "$(date -Is) driver shard$SHARD cycle $cycle start (expected_samples=$EXPECTED)" >> "$LOG"
  run_cycle "$cycle"
  # Count only VERIFIED completions (review H8): a truncated unified.tsv used
  # to satisfy "file exists" here and in the engine's skip test.
  ACTUAL=$(python3 -I scripts/verify_batch.py --results results --print-ok-count 2>/dev/null || echo 0)
  [ -z "$ACTUAL" ] && ACTUAL=0
  echo "$(date -Is) driver shard$SHARD cycle $cycle done — unified=${ACTUAL}/${EXPECTED}" >> "$LOG"
  if [ "$ACTUAL" -ge "$EXPECTED" ]; then
    echo "$(date -Is) driver shard$SHARD finished cleanly" >> "$LOG"
    exit 0
  fi
  cancel_leftovers
  sweep_stale_tmp
  BACKOFF=$((300 * cycle))
  echo "$(date -Is) $((EXPECTED - ACTUAL)) sample(s) still missing — next cycle in ${BACKOFF}s" >> "$LOG"
  sleep "$BACKOFF"
done

echo "$(date -Is) driver shard$SHARD exhausted $MAX_CYCLES cycle(s); unified=$(find results -name '*.unified.tsv' | wc -l)/$EXPECTED — manual review needed" >> "$LOG"
exit 42
