# `oxo-flow cluster *` audit on a real site: status/logs are not run-aware, `logs` needs `sacct`, and three submission paths diverge

**Environment**: oxo-flow v0.30.0 (`384db709`), Slurm `workq`, 4 nodes.
**Deployment**: the pipeline (`tcga-mixcr`) runs from four workflow directories
(one per node) via `oxo-flow run --profile slurm` under a per-shard driver
wrapper; each directory expands 2,780 samples × 4 rules ≈ 11.1k rule
instances (~44k cluster-wide), jobs last 10 min – 3 h, a full pass runs for
weeks.

This is a consolidated report for the whole `cluster` command family; it
supersedes the status-only report filed earlier as #1085.

## What the family offers today

```text
$ oxo-flow cluster --help
Commands:
  submit  Generate and submit job scripts for a workflow to the cluster
  status  Show the status of submitted cluster jobs ...
  cancel  Cancel submitted cluster jobs
  logs    Fetch logs for a submitted cluster job
```

Three different submission/monitoring models coexist:

| path | state model | monitoring that understands it |
|---|---|---|
| `run --profile slurm` (driver-integrated) | run directory: `events.jsonl`, `jobs/<instance>/status.json`, `job.id`, checkpoint | none |
| `cluster submit` (standalone script generator) | `cluster_scripts/` + scheduler dependencies | none |
| `cluster status/cancel/logs` | raw scheduler job ids | — |

The monitoring commands are anchored to the loosest model (raw ids), so they
cannot answer any question about a `run --profile` deployment.

## Findings

### C1 — `status` returns `id` + scheduler state only, unscoped, and prints the list twice

```text
$ oxo-flow cluster status -v
Cluster: Executing 'squeue -u wsx --noheader -o %i|%T'...
Cluster: 20 queued jobs
  1400_1: running
  ...
1400_1  running
...
```

- no rule name / sample / node / elapsed / CPU / memory / progress;
- no scoping to a workflow or run directory — with four shards every
  invocation lists **all** of the user's jobs identically;
- the human listing (stderr) and the machine listing (stdout) are both shown
  in a TTY, so the list appears twice (`--quiet` suppresses one copy, but the
  help does not say so).

### C2 — `status` cannot map a job id back to a rule instance

`oxo-flow cluster status 1432` says `1432 running`. The run directory already
contains `jobs/<rule_instance>/job.id` files, so the mapping exists locally;
the command just does not use it.

### C3 — `logs` requires Slurm accounting and fails outright without it

```text
$ oxo-flow cluster logs 1432
Error: failed to fetch job logs
Caused by:
    config error: 'sacct' exited 1: Slurm accounting storage is disabled
```

Many small Slurm sites (including this one) run without `slurmdbd`. The
command has no fallback to (a) `scontrol show job`, or (b) the logs the engine
itself writes into the run directory (`logs/<rule>.<jobid>.out/.err`), which
are available for every submitted job. It also accepts only scheduler ids,
not rule-instance names.

### C4 — `cancel` is id-only

`cluster cancel 1432 1433` works, but there is no way to cancel "this run" or
"the rule instance X" from the workflow side, and no dry-run/preview of what
would be cancelled.

### C5 — no `--json` anywhere in the family

`cluster status --json` → `Error: 'oxo-flow cluster' does not support --json
yet`. Automation must scrape a tab-separated list, and `logs` cannot be parsed
at all.

### C6 — the companion status command reports zeros during a cluster run

`oxo-flow status` reads `.oxo-flow/checkpoint.json`, which during a cluster
run still says `Completed: 0` because completions are folded only after the
whole run returns (see the restart-correctness issue). So neither `status`
nor `cluster status` answers "how many samples are done?" while the run is
live.

### C7 — job names carry no instance identity

Every submitted job is named after the rule (`mixcr`, `trust4`, …), so
`squeue` itself is unreadable on a multi-shard host and two runs are
indistinguishable without inspecting `scontrol show job` per id.

## Data that already exists (so this is a presentation problem, not a data problem)

```
.oxo-flow/runs/<ts>/events.jsonl
  {"job":"1400_1","rule":"trust4_all_<uuid>","t":"SUBMITTED|COMPLETED|FAILED|CANCELLED|SKIPPED","ts":...}
.oxo-flow/runs/<ts>/jobs/<rule_instance>/status.json
  {"state","job_id","command","submitted_at","finished_at","exit_code"}
.oxo-flow/runs/<ts>/jobs/<rule_instance>/job.id, job.sh
logs/<rule>.<jobid>.out / .err        ← captured stdout/stderr, present for every job
<workflow> per-rule resources + expanded DAG (planned instance count)
```

## Proposal (one coherent run-aware surface)

1. **Run-aware by default.** `cluster status|logs|cancel` should resolve the
   newest run directory of the workflow in the current directory (or
   `--run-dir`), map scheduler ids ↔ rule instances through `job.id`, and add
   `--all` for today's raw behaviour.
2. **Status table** per instance with aggregates:

   ```text
   tcga-mixcr @ shards/shard1 · run 2026-11-16T12:06
   planned 11,120   completed 34   running 16   pending 11,070   failed 1

   rule          sample       state    node   elapsed/limit  cpu  mem   progress
   mixcr         0056507f…    running  node1  2:24:11/12h     4    48G   alignment 93%
   trust4       01ea7cb3…    running  node3  0:12:02/12h     4    24G   extracting reads
   standardize  12183b9e…    pending  –      –              1    4G    waiting for trust4
   ```

3. **`logs` fallback chain**: run-dir capture → `scontrol`/`qstat -f` →
   `sacct`, with `--tail N`; accept rule-instance names as well as ids.
4. **`cancel`** should accept `--run`, `--rule <name>`, and instance patterns,
   with a `--dry-run` preview.
5. **`--json`** for all four subcommands with the same document shape as the
   table.
6. **TTY-aware narration** so the status list is not printed twice.
7. **Job names**: include the rule instance (or allow a `--job-name`
   template), so `squeue` is legible on multi-workflow hosts.

## Workaround in production

We render our own dashboard from the run directories
(`scripts/shards_status.py`: planned/submitted/completed/failed per shard and
per rule, in-flight jobs with node/elapsed/CPU, plus node temperatures), and
read logs from the run directory. It works, but every operator of a
long-running `--profile slurm` deployment has to rebuild the same view.

Happy to test a patch against the live 4-shard deployment.
