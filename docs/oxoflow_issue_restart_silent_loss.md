# oxo-flow issue draft — cluster restart can silently skip rules whose jobs died

> Discovered 2026-11-16 while validating the Slurm cluster backend for an
> 11k-sample production pipeline. Engine: oxo-flow v0.30.0 (commit 384db709).
> Tests: `tests/oxoflow_cluster/` (synthetic, 12 samples, `sleep 45` rules).
> Status: root cause isolated — output-existence skipping. A state-based
> mechanism was ruled out by T3e (below).

## Summary

On `oxo-flow run --profile <slurm>`, if the driver process dies and its
in-flight Slurm jobs also die (node failure, walltime, `scancel`, maintenance),
the restarted driver may **not** resubmit those rules. The run then completes
with missing outputs and exits successfully. Nothing in the log names the
skipped rules.

## Reproduction

Profile: `backend=slurm`, `max_submitted=2`, `--nodelist=<node>`.

1. Start the workflow; wait until 2 jobs are RUNNING.
2. `kill -9 <driver pid>` and `scancel` the two RUNNING job ids.
3. Restart the identical command.

Observed across four variants (T3b/T3c/T3d, 2026-11-16): the restart submits
**10 of 12** jobs. The two cancelled samples never run.

## Root cause (isolated): output existence is treated as completion

`crates/oxo-flow-cli/src/commands/run_preview.rs:178 rule_outputs_exist_fresh`
skips a rule when every declared output exists (and is newer than its inputs).

A rule that creates its declared output **early** (or a tool that streams into
the final path) leaves a partial file when killed. The restart then treats the
rule as complete:

```
$ cat proto3/out/S12/done.txt        # job 1259 was scancelled mid-run
start 2026-11-15T22:36:55-05:00 host=node4 job=1259 task=2
                                      # <- no "done" line; file still counts
```

The complementary sample (`S11`, also killed) was re-run and completed because
its partial file had been removed — confirming the skip is driven by file
existence, not by scheduler state.

## State-based suppression is NOT the cause (T3e)

A second candidate was the recorded `PENDING` state in
`runs/<ts>/jobs/<rule>/status.json`. T3e isolated it: a rule whose output is
written only at the end, killed while RUNNING with **zero** output files, then
restarted with **no state removed** — the restart submitted **12/12** and all
completed. The recorded `PENDING` entries are re-queued, not adopted as
complete. Only output existence suppresses a rule.

## Impact

A long production driver is exposed through either path: a killed job can be
silently skipped, and the workflow reports success. Detecting this requires an
external expected-vs-actual output accounting.

## Suggested fixes

1. Verify **job liveness** when adopting `PENDING`/`RUNNING` records at
   startup: if the scheduler no longer knows the job and the outputs are
   absent *or older than the job's start*, re-queue the rule.
2. Persist `mark_completed` incrementally on the cluster path (the
   `DriverOptions::on_checkpoint` hook exists; `run_cluster.rs` passes
   `None`), so a restart has truthful data instead of submit-time
   placeholders.
3. Document that output existence is not a completion proof, and support a
   completion marker (e.g. a declared sentinel written last), which the skip
   logic can require.
4. Never report success when a rule was dropped; at minimum count and print
   skipped-without-output rules at the end of the run.

## Workaround used in production (pipeline-side, no engine change)

- Every rule declares a **sentinel output** written as its last action, so
  "all outputs exist" implies completion.
- The driver moves `.oxo-flow/checkpoint.json` and `.oxo-flow/runs/` aside on
  restart, cancels leftover jobs recorded in the old run directory, waits for
  the scheduler to drop them, then re-plans from scratch.
- An output-count guard turns "rc=0 but samples missing" into a retry.

Cost: a full plan per restart (noticeable at 2.8k samples).
