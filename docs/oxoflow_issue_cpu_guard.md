# Add a CPU-oversubscription guard: warn when planned concurrency × rule threads saturates detected CPUs; let `[resource_budget] max_threads` apply to local runs

## Context

Real incident, oxo-flow **v0.30.0** (`384db709`), used as the orchestrator for an
11k-sample bioinformatics cohort on a 4-node Slurm cluster. Each node runs one
shard workflow **locally** (inside a `sbatch` allocation), 64 CPUs per node.

Workflow rules declare `resources = { threads = 8, memory = "32G" }`.
We launched with `oxo-flow run workflow.oxoflow -j 8 --provenance`.

## What happened

`-j 8` scheduled up to **8 concurrent rules × 8 threads = 64 threads on a
64-CPU node** — full core saturation for hours. Nodes sat at load ~50–52,
package temperatures reached 70 °C, and the machine room raised a **fan-noise
complaint from operations**. Nothing in the oxo-flow output warned that the
requested concurrency would consume every core of the host.

## Evidence that the engine already has all the information

1. The engine tracks per-rule thread demand precisely (it pools threads):

```
WARN waiting for resources: memory: need 32768MB, have 28672MB;
top thread holders: [mixcr_all_004064de-...(8t), mixcr_all_004589d9-...(8t),
mixcr_all_00564877-...(8t)]; top memory holders: [...]
rule=mixcr_all_00f54652-... waited_secs=245
```

2. It detects host capacity (`--max-threads 0`/`--max-memory 0` auto-detect,
   "Auto-detected capacity is soft", and the cgroup memory warning shows host
   detection working):

```
WARN cgroup memory limit is below the detected host memory — clamping the
auto-detected ceiling (issue #855: scheduling past the cgroup cap OOM-kills
the run; override with --max-memory if this is wrong) host_mb=1417875 ...
```

3. For **memory** the engine warns and clamps. For **CPU** there is no
   equivalent signal at 100 % demand — we grepped the entire run log for
   CPU/thread warnings: none besides the resource-wait line above.

So the engine knew "8 rules × 8 threads" and knew the host has 64 CPUs, but
never said "this plan will saturate all cores."

## Requests

### R1 — plan-time oversubscription warning

When planning/starting a run, compute the achievable thread demand
(e.g. job concurrency × max declared rule threads over the ready set, or the
sum of currently runnable rules' threads) and compare it with the detected
CPU count (cgroup-aware like the memory path). If demand reaches, say, ≥ 80 %
of detected CPUs, emit a WARN with the numbers and a one-line remedy:

```
WARN planned thread demand 64 = concurrency 8 × max rule threads 8, vs 64
detected CPUs (100%) — the run will occupy all cores; pass --max-threads N
or lower -j to leave headroom for the OS and co-tenants
```

This is the exact analogue of the existing cgroup-memory WARN, and it matters
most on shared HPC nodes / desktops where "how much of the box do I take" is
an operational (thermal, noise, fairness) decision, not just a performance one.
Suggested pairing: when auto-detection is in effect, this warning should be on
by default; a `--max-cpu-percent`-style percentage cap would be a natural
companion to `--max-threads` (absolute) and `--max-memory`.

### R2 — `[resource_budget] max_threads` (and `max_memory`) should also cap **local** runs

The workflow-format reference states:

> `[resource_budget]` — max_jobs clamps concurrency on both paths;
> max_threads/max_memory enforced on cluster submissions — local runs use
> `--max-threads`/`--max-memory`

Because our shards run locally inside `sbatch` (not through oxo-flow's cluster
backend), the 20-CPU cap we were forced to adopt **cannot live in the versioned
workflow** — it has to sit in an external `run_shard.sh` wrapper. That breaks
the "configuration lives in the workflow, so it is reproducible and reviewable"
property. If `[resource_budget] max_threads = 20` clamped local runs too, the
cap would travel with the `.oxoflow` file (and `validate`/`dry-run` could
surface it).

## What we do meanwhile (not a blocker)

`--max-threads 20` works as documented (pre-flight fast-fail + pool clamping);
we now set `-c 20` in `sbatch`, `-j 2`, and `--max-threads 20`, and verified
temps fall back to 43–59 °C. The problem is **discoverability**: we learned
about it from an operations complaint, not from the engine.

## Environment

- oxo-flow 0.30.0 (`384db709`), local execution inside Slurm allocations
- 4 × 64-CPU nodes, 257 GB–1.4 TB RAM each; pixi backend
- Rules: `resources = { threads = 8, memory = "32G" }`, `-j 8` (before),
  `-j 2 --max-threads 20` (after)

Related but distinct: #1078 (placeholder-validation and lifecycle gaps).
