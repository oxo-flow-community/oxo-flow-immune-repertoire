# Cutover runbook — sbatch supervisor → oxo-flow Slurm cluster backend

> **STATUS 2026-11-16 06:0x**: cutover EXECUTED for all four shards. The old
> `supervisor.sh` is stopped; each shard runs `scripts/driver_shard.sh N`
> (flock, retry loop, output-gap guard) with the sentinel workflow and
> `profiles/slurm.toml` (`--nodelist=nodeN`, `max_submitted=2`, 8 threads →
> 16 CPUs/node). Rule definitions are frozen: editing them invalidates every
> completed instance and forces a re-run.

> Scope: production Phase 5 (shards 1–4 on node1–4). Written 2026-11-16 after
> the prototype and tests T1–T4. Every step is reversible; the rollback is at
> the end.

## Preconditions (do not start otherwise)

- [ ] T1 prototype finished 72/72 rule instances (mixcr + trust4 + standardize)
      with the sentinel workflow.
- [ ] No `FAILED` events in any shard's run directory.
- [ ] `/data` has ≥ 5 TB free (`df -h /data`), because re-planned rules can
      briefly duplicate intermediates.
- [ ] Cluster load policy check: each node ≤ 20 CPUs for this project.

## 1. Freeze the old pipeline

```bash
ssh <login-node>
# stop the supervisor so it cannot resubmit the sbatch-based shard jobs
pkill -f "supervisor[.]sh"
pgrep -af "supervisor[.]sh" || echo "supervisor stopped"
# note which shard jobs are running, then cancel them
squeue -u wsx -o "%i %j %T %R"
scancel -n tcga-shard1 -n tcga-shard2 -n tcga-shard3 -n tcga-shard4
```

`oxo-flow` writes its local-mode checkpoint after each completed rule, so the
work already done is recorded in `shards/shardN/.oxo-flow/checkpoint.json`.

## 2. Prepare all shards (dry run, then apply)

```bash
cd $HOME/Projects/TCGA-mixcr
bash scripts/cutover_cluster_backend.sh            # inspect the report
bash scripts/cutover_cluster_backend.sh --apply
```

This deploys the sentinel workflow, writes `profiles/slurm.toml`, backfills
completion sentinels for finished samples (so nothing is re-run) and runs the
validate/lint gates per shard.

## 3. Start the drivers

```bash
for n in 1 2 3 4; do
  setsid nohup $HOME/Projects/TCGA-mixcr/scripts/driver_shard.sh $n \
      >/dev/null 2>&1 < /dev/null &
  sleep 3
done
pgrep -af "driver_shard[.]sh"
```

## 4. Verify (first 30 minutes)

```bash
# jobs are per-rule, 8 CPUs each, 2 in flight per node -> 16 CPUs/node
squeue -u wsx -o "%.10i %.28j %.8T %.10M %R"
scontrol show job <id> | grep -E "NumCPUs|TRES|NodeList"

# each node shows only its own shard's jobs
for n in 1 2 3 4; do echo "node$n: $(squeue -h -w node$n | wc -l) jobs"; done

# sentinels appear as rules finish
find shards/shard1/results -name "*.mixcr.done" | wc -l

# drivers are alive and logging
tail -5 logs/driver_shard1.log
```

Expected: 2 jobs per node, `NumCPUs=8` each, `--nodelist` fixed to that node,
no job on a node from another shard.

## 5. Rollback (any time)

```bash
pkill -f "driver_shard[.]sh"
scancel -u wsx                     # cancel driver-submitted rule jobs
cd $HOME/Projects/TCGA-mixcr/shards
setsid nohup ./supervisor.sh >/dev/null 2>&1 < /dev/null &   # old path
squeue -u wsx -o "%.10i %j %T %R"   # expect tcga-shardN jobs again
```

The old checkpoint is untouched by the new path (the cluster run keeps its own
run directories), and the new workflow's extra sentinel outputs are additive,
so the old path still sees its declared outputs present.

## 6. Post-cutover monitoring

- `logs/driver_shardN.log` — one line per start/exit, plus the
  `unified=<actual>/<expected>` accounting line after every run.
- `shards/shardN/.oxo-flow/runs/latest/events.jsonl` — SUBMITTED/COMPLETED/
  FAILED stream.
- Driver restarts: only triggered by a non-zero exit or an output gap; each
  retry cancels leftover jobs and re-plans.
- Disk: `du -sh shards/shardN/tmp` should stay near zero between rules now
  that the `cleanup` rule is in the workflow.
