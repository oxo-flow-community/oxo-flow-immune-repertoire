# `cluster status` is not actionable during a long run: job id + state only, unscoped, no `--json`, and no way to see per-rule progress

**Environment**: oxo-flow v0.30.0 (`384db709`), Slurm (`workq`, 4 nodes).
**Context**: a production pipeline runs 4 independent workflow directories in
parallel (one per node), each expanding 2,780 samples × 4 rules into ~11,120
rule instances (`mixcr` → `trust4` → `standardize` → `cleanup`), driven by
`oxo-flow run --profile slurm` under a per-shard wrapper. Jobs run for 2–4
hours each; a full pass is measured in weeks.

## Observed

```text
$ oxo-flow cluster status -v
Cluster: Executing 'squeue -u wsx --noheader -o %i|%T'...
Cluster: 20 queued jobs
  1400_1: running
  1405_1: running
  ...
1400_1  running
1405_1  running
...
```

Same with a single id:

```text
$ oxo-flow cluster status 1432
Cluster: Executing 'squeue -j 1432 --noheader -o %i|%T'...
Cluster: 1 job:
  1432: running
1432  running
```

## Findings

| # | Finding |
|---|---|
| F1 | **Only `job id` + scheduler state.** No rule name, sample id, node, elapsed/limit, CPU/mem reservation, or last log line. With four shards and four rule kinds, an operator cannot tell which sample a job belongs to, how far a rule has progressed, or where the long tail is. |
| F2 | **Unscoped.** Without job ids the command lists every job the user owns (`squeue -u <user>`), so jobs from different workflow directories are indistinguishable — the same command run from any shard shows the same list. |
| F3 | **Duplicated output in a terminal.** The human listing goes to stderr and the machine list (`<id>\t<state>`) to stdout by the codebase convention (`cluster.rs`: "human output stays on stderr"). In an interactive TTY both are visible and the whole list appears twice. `--quiet` suppresses the stderr half, but neither the help text nor the duplication itself hints at that. |
| F4 | **`--json` is rejected** for `cluster status` (`oxo-flow cluster does not support --json yet`), so automation must scrape the tab list. |
| F5 | **No alternative progress view.** `oxo-flow status` reads `.oxo-flow/checkpoint.json`, which during a cluster run still says `Completed: 0` (completions are folded only after the whole run returns — see #1083), so neither command answers "how many samples are done?". |

## What the data already available looks like

The run directory already holds everything needed to answer F1/F2/F5:

- `.oxo-flow/runs/<ts>/events.jsonl` —
  `{"job":"1400_1","rule":"mixcr_all_<uuid>","t":"SUBMITTED|COMPLETED|CANCELLED|SKIPPED","ts":...}`
- `.oxo-flow/runs/<ts>/jobs/<rule_instance>/status.json` —
  `{"state","job_id","command","submitted_at","finished_at","exit_code"}`
- `.oxo-flow/runs/<ts>/jobs/<rule_instance>/job.id`, `job.sh`
- the workflow itself: per-rule resources, declared inputs/outputs, and the
  expanded DAG (planned instance count)

## Proposal

1. **Scope to a run.** Default to the newest run directory in the current
   workflow (or `--run-dir <path>`); add `--all` for today's "every job of the
   user" view. Map scheduler job ids back to rule instances through the
   `job.id` files, so `cluster status <id>` explains *which rule of which
   sample* the id is.
2. **One table, per instance**, e.g.:

   ```text
   tcga-mixcr @ shards/shard1 · run 2026-11-16T07-12
   planned 11,120   completed 34   running 16   pending 11,070   failed 0

   rule          sample       state    node   elapsed/limit  cpu  mem   progress
   mixcr        0056507f…    running  node1  2:24:11/12h     4    48G   alignment 93%
   mixcr        008e5b39…    running  node1  0:40:06/12h     4    48G   assembly
   trust4       01ea7cb3…    running  node3  0:12:02/12h     4    24G   extracting reads
   standardize  12183b9e…    pending  –      –              1    4G    waiting for trust4
   ```

   with an aggregate header per rule (`mixcr 12/2780 done, 4 running, 2764
   pending`) and failure counts pulled from `events.jsonl`.
3. **`--json`**: emit the same document (per-instance records + aggregates).
   The run-dir files already serialize this; the command only has to join them
   with `squeue`/`sacct` output.
4. **Job naming**: include the rule instance in the scheduler job name (today
   every job is named `mixcr`, `trust4`, …), or expose a `--job-name`
   template, so `squeue` itself is readable and multi-shard hosts can be
   separated without the run directory.
5. **TTY detection** (optional): when stdout is a terminal, skip the stderr
   copy of the list (or vice versa) so the duplication in F3 disappears
   without requiring `--quiet`.

## Impact

For a multi-week production run the operator question is "how many samples are
done, which stage is each sample in, what is failing, and where is the long
tail". Today the answer requires scripting `squeue` + `events.jsonl` +
`status.json` by hand (which is what we do), and the engine's own two status
commands both fail to answer it — one shows ids only, the other shows zeros.

Happy to test a patch on our live 4-shard deployment; we can provide the run
directory layout and a captured `squeue`/`sacct` sample.
