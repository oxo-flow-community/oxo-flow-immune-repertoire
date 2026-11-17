# T1 — cluster-backend end-to-end validation (2026-11-16)

> Scope: prove that the oxo-flow Slurm cluster backend can run the production
> workflow (mixcr → trust4 → standardize → cleanup) on real TCGA samples, with
> correct outputs, correct resource use, and correct restart behaviour.
> Hosts: node4 (Slurm `workq`), 12 samples, profile `max_submitted=2`,
> per-rule `threads=8 / memory=32G` → 16 CPUs (policy ≤ 20).

## What ran

Final run directory `proto/.oxo-flow/runs/2026-11-16T04-06-21-530`:

| Rule | Instances planned | Skipped (sentinel) | Submitted | Completed | Failed |
|---|---:|---:|---:|---:|---:|
| mixcr | 12 | 12 | 0 | – | – |
| trust4 | 12 | 1 | 11 | 11 | 0 |
| standardize | 12 | 0 | 12 | 12 | 0 |
| cleanup | 12 | 0 | 12 | 12 | 0 |
| **total** | **48** | **13** | **35** | **35** | **0** |

Note: the 12 MiXCR instances were produced by the earlier pre-sentinel
workflow (jobs exit 0, 20:14–22:45) and their completion sentinels were
backfilled so they were not re-run. The MiXCR **rule body** under the new
workflow is exercised by the shard4 canary, not by this T1 run.

## Independent verification (read-only, separate agent)

| Check | Result | Evidence |
|---|---|---|
| Sentinel completeness | PASS | 48/48 sentinels (12 × 4 rule types) |
| unified header | PASS | 14 columns, canonical order, 12/12 |
| per-(tool,chain) sum(clone_fraction) | PASS | max deviation 1.4e-06 (%.6g rounding) |
| row conservation | PASS | 12/12 exact; full multiset comparison of (chain, nt, aa, v, j, read_count) matches source files |
| spot value check | PASS | 24/24 rows identical to source |
| cleanup | PASS | 12/12 `.tmp_cleaned`; `tmp/` = 0 B |
| event stream | PASS | SUBMITTED 35 = COMPLETED 35, FAILED 0; ≤2 concurrent tasks (16 threads) |
| job resources | PASS | `job.sh`: `--cpus-per-task=8 --mem=32G --nodelist=<node>` |
| directory hygiene | PASS | no >100 MB files outside `results/`; `/dev/shm` empty |

Known limitations (all addressed):

1. MiXCR rule body not executed in this run (see note above) — covered by the
   shard4 canary.
2. TRUST4 can leave `*_toassemble.fq` / `*_assembly.fa` (100–200 MB each) in
   `results/` when a run is interrupted → added `rm -f` to the trust4 rule.
3. Prototype debris (`READ/` leftovers, 2 orphan `.fq`, ~340 MB) removed.

## Restart behaviour (tests T3–T3f, synthetic)

| Test | Scenario | Result |
|---|---|---|
| T3 | kill driver, jobs already finished | restart submitted 10/12, **zero duplicate work**, 12/12 outputs |
| T3b | kill driver **and** scancel in-flight jobs | 2 samples silently skipped (partial output counted as complete) |
| T3c | + delete checkpoint only | still 10/12 — insufficient |
| T3d | + move checkpoint AND run directory | still 10/12 — root cause was the partial output file |
| T3e | kill+cancel with **zero outputs**, no state moved | restart submitted **12/12** — recorded PENDING state is re-queued, not adopted |
| T3f | sentinel pattern (payload first, sentinel last) | restart submitted **12/12**, killed samples re-ran, 12/12 sentinels |

**Conclusion**: the engine skips a rule iff all declared outputs exist, so every
rule must declare a LAST-written sentinel. The driver additionally cancels
leftover jobs and waits for them to leave the queue before re-submitting.

## Gates run

`oxo-flow validate | dry-run | lint` on both workflow shapes:
4 rules / 48 instances (12 samples) — 0 errors.
