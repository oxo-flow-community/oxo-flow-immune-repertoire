> Filed as https://github.com/Traitome/oxo-flow/issues/1083

# Cluster backend: partially-written outputs are treated as completion — rules are silently skipped after a driver/job crash

**Environment**: oxo-flow v0.30.0, commit `384db709`, Slurm (`workq`, 4 nodes),
reported from a production rollout of an 11,123-sample RNA-seq pipeline
(MiXCR + TRUST4) using `oxo-flow run --profile <slurm>`.

**Severity**: high for any long-running cluster job. A rule whose Slurm job was
killed can be permanently skipped while the run still exits 0.

---

## 1. Summary of findings

| # | Finding | Severity |
|---|---|---|
| F1 | A rule is skipped when its **declared outputs exist**, so a killed job that wrote a partial output is never re-run. Combined with F2 this silently drops samples. | high |
| F2 | Per-rule **completion state is not persisted during a cluster run**; the checkpoint still shows `completed_rules: []` after a crash. The restart therefore has no trustworthy completion record and falls back to F1's heuristic. | high |
| F3 | `--provenance` is accepted but not honoured on the cluster path (engine warns). | low |
| F4 | Planning cost at scale: ~12 min for 10–11k rule instances, paid again on every driver restart; no incremental planning. | medium |
| F5 | Workflow-shape errors surface as raw parse failures: deploying a workflow whose `metadata_file` path does not exist yields `parse error ... No such file or directory` with no hint about the intended shape (`sample_groups_file`, etc.). | low |

F1+F2 are one defect chain and are the main subject below. F3–F5 are smaller
and independent; evidence is included so they can be triaged separately.

---

## 2. F1 + F2 — silent skip after a crash

### 2.1 Reproduction (synthetic, no real compute)

Workflow (12 samples from `manifests/groups.tsv`, one rule per sample):

```toml
[workflow]
name = "cluster-synth"
sample_groups_file = "manifests/groups.tsv"

[defaults]
shell_prelude = "set -euo pipefail"

[[rules]]
name = "emit"
resources = { threads = 1, memory = "1G" }
output = ["out/{sample}/done.txt"]
target = true
# NOTE: this rule creates its declared output at the START (payload first,
# "done" line appended at the end) — the same shape as real tools that stream
# into their final path.
shell = "mkdir -p out/{sample} && echo \"start $(date -Is) job=${SLURM_JOB_ID:-na}\" > {output[0]} && sleep 45 && echo \"done $(date -Is)\" >> {output[0]}"
```

Profile: `backend=slurm`, `partition=workq`, `max_submitted=2`,
`extra_args=["--nodelist=node4","--export=ALL"]`.

Steps:

1. Start the run; wait until 2 jobs are RUNNING.
2. `kill -9 <driver pid>` **and** `scancel` the two running job ids.
3. Restart the identical command.

### 2.2 Observed results (four variants)

| Variant | Restart preparation | Jobs submitted on restart | Outcome |
|---|---|---:|---|
| A | none (jobs had already finished) | 10 / 12 | correct: the 2 completed rules were skipped, no duplicate work |
| B | none (2 jobs were killed mid-run) | 10 / 12 | **2 samples never ran again**; run exited 0 |
| C | `rm .oxo-flow/checkpoint.json` | 10 / 12 | still skipped |
| D | moved `checkpoint.json` **and** `.oxo-flow/runs/` aside | 10 / 12 | still skipped |
| E | killed with **zero output files**, state untouched | **12 / 12** | proves recorded `PENDING` state is re-queued, not adopted |
| F | rule declares a **sentinel output written LAST** | **12 / 12** | killed samples re-ran; mitigation validated |

Decisive evidence from variant D — the killed job's partial output was treated
as a completed result:

```
$ cat out/S12/done.txt
start 2026-11-15T22:36:55-05:00 host=node4 job=1259 task=2
                                       # no "done" line, yet the rule was skipped

$ cat out/S11/done.txt                 # also killed, but its file had been
start ... job=1261                     # removed before the restart -> re-ran
done  2026-11-15T22:37:54-05:00 ...
```

Run directories / job scripts for each variant are on the reporting side
(`proto3/logs/t3*.log`); the workflow above is the complete reproduction.

### 2.3 Root cause (source references)

1. **Skip test is output existence only** —
   `crates/oxo-flow-cli/src/commands/run_preview.rs:178` `rule_outputs_exist_fresh`
   (and `oxo-flow-core/src/executor/checkpoint.rs:1986`
   `all_outputs_exist = expanded_outputs.iter().all(|o| workdir.join(o).exists())`).
   Nothing checks that a rule *finished*; a partially written file passes.

2. **Completion is folded into the checkpoint only after the whole run returns** —
   `crates/oxo-flow-cli/src/commands/run_cluster.rs`:
   - the submit hook (`on_submit`, ~line 570) persists a `running` placeholder
     per submission;
   - `record_outcome(...)` → `mark_completed(...)` (line ~341) /
     `save_to_file_async` (line ~386) is called from a loop **after**
     `driver.run(...)` returns (line ~642).
   After a crash the checkpoint therefore contains
   `completed_rules: []` and `rule_runs: {<rule>: running}` — no truth about
   what finished, forcing the output-existence heuristic.

3. `DriverOptions::on_checkpoint` exists but is passed as `None` by the
   cluster path, so there is no incremental fold even though the driver
   supports it.

### 2.4 Suggested fixes (in order of preference)

1. **Fold outcomes incrementally on the cluster path** — drive
   `DriverOptions::on_checkpoint` from the same `record_outcome` logic so
   `completed_rules` is persisted as each job settles. A restart then has
   authoritative data and never needs the heuristic.
2. **Do not trust output existence for rules recorded as `running`/`pending`
   at startup.** Verify the scheduler: if the recorded job id is unknown and
   no completion record exists, re-queue the rule (the jobs are idempotent by
   contract for this use case).
3. **Support/document a completion marker**: let a rule declare a "last
   written" sentinel, and require it (or a `checksum` entry) for the skip
   test. We adopted this as the user-side workaround and it closed the
   failure mode completely (variant F).
4. Log skipped-without-a-record rules in the terminal summary instead of
   silently counting them as complete.

---

## 3. F3 — `--provenance` is a no-op on the cluster path

`oxo-flow run --profile slurm --provenance` prints a warning that provenance
is not supported on the cluster path and produces no provenance record.
Observed on every cluster run; local runs are unaffected.

Suggested: either implement (the run directory already has job scripts,
status and events) or reject the flag with an error instead of a warning.

---

## 4. F4 — planning cost at production scale

`Cluster: submitting 11124 jobs to slurm (max 2 in flight)` — the planning
phase for 2,781 samples × 4 rules took ~12 minutes at ~1–3 cores before the
first submission. The same cost is paid on every driver restart (we do not
rotate state, so this is not avoidable by us).

Measurements: 2,781 samples × 3 rules ≈ 12 min; 2,781 × 4 rules ≈ 12 min
(same order). Single-node, 1.4 TB RAM, no memory pressure.

Suggested: incremental/persisted plan, or submission streaming (start
submitting the first ready rules while later instances are still planned).

---

## 5. F5 — workflow-shape parse error has no guidance

Deploying a workflow that references `metadata_file = "manifests/x.tsv"` where
the file is absent fails at parse time with:

```
parse error in ./manifests/test_samples_metadata.tsv: failed to read metadata
for metadata file: No such file or directory (os error 2)
```

There is no hint that the file is missing (vs. malformed), nor which key was
used. In our case the repo carries two shapes of the same workflow (test shape
with inline `[[sample_groups]]`, shard shape with `sample_groups_file`), and a
mistaken deployment produced exactly this error.

Suggested: name the missing path and the key that referenced it in the error
message; mention the alternative `sample_groups_file` form.

---

## 6. Workaround we ship meanwhile (user side, no engine change)

- every rule declares a sentinel written as its **last** action
  (`.mixcr.done`, `.trust4.done`, `.unified.done`), so "all outputs exist"
  implies "rule finished";
- the driver cancels leftover Slurm jobs from the previous run and waits for
  them to leave the queue before re-submitting;
- an output-count guard turns `rc=0` with missing per-sample outputs into a
  retry.

This has run 48 real rule instances with 0 failures, and a 11,124-job
production submission is in progress at the time of filing.

Related, already filed and closed upstream: #1078 (expansion semantics /
placeholder validation / lifecycle gaps), #1080 (CPU-oversubscription guard).
This issue is the cluster-path counterpart found after those fixes were
requested.
