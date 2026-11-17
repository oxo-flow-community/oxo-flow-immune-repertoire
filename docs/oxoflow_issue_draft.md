# Usability report from a real-world pipeline rollout (v0.30.0): expansion semantics, placeholder validation, and lifecycle gaps

## Context

We deployed oxo-flow v0.30.0 (commit `384db709`, built from source) as the orchestrator for a TCGA immune-repertoire pipeline (11,123 samples; this report covers the 8-sample validation phase): per-sample rules for MiXCR 4.7 / TRUST4 1.1.5 over GDC RNA-seq BAMs, pixi backend, sample domain via `[workflow] metadata_file` + `[[sample_groups]]`, runs with `-j 2 --provenance`.

Everything below was hit in real runs (log excerpts are verbatim from `oxo-flow run` output). Overall the core worked well — see **Positives** at the end. These are the friction points, ordered by impact.

---

## F1 [HIGH] Sample-domain expansion only triggers on `{sample}` in `input`/`output`; failure surfaces at run time with literal placeholders

**Symptom.** A rule whose `output` uses only `{meta.*}` (and `{sample}` only in `log`/`shell`) did not expand per sample. It ran once with literal `{meta.*}` tokens:

```
WARN unresolved placeholder remained in the rendered command — no value bound to it in this run
     (the rule will execute with a literal brace token) ... placeholder="{sample}" rule=trust4
...
✗ rule 'trust4' failed
stderr: ... Could not find file {meta.bam_genomic}
```

The per-sample instance (e.g. `trust4_test_<uuid>`) never entered the DAG — the rule produced exactly **one unexpanded instance**. `dry-run` faithfully showed that single instance, which is correct output, but easy to misread as expanded. (Why the WARN names `placeholder="{sample}"` while the failing token was `{meta.bam_genomic}`: with no expansion there is no per-sample binding context, so neither token resolves.)

**Suggested fix.** Catch this at the static gates, without breaking legitimate singletons:
1. At `lint`: warn when a rule references `{meta.*}`/`{sample}`/`{group}` anywhere but has no sample-domain binding in `input`/`output` — producing one unexpanded instance in which those tokens stay unbound. Single-instance rules reading constant metadata columns are a valid pattern, so this should be a warning by default (with a per-rule opt-out, e.g. `singleton = true`, to silence intentional singletons).
2. At render time (see F2): any token still unbound after rendering is an error, which is what actually saves the run.

## F2 [HIGH] Unresolved config/meta/input/output placeholders are WARN-only, then the command executes and fails downstream

Related to F1 (think of F1 as the diagnosis, F2 as the gate) but broader: tokens matching the engine's declared placeholder grammar left unbound after rendering — `{sample}`, `{group}`, `{meta.*}`, `{config.*}`, `{input.*}`, `{output.*}` — only produce a warning today, then execute (and the tool fails on a file literally named `{meta.bam_genomic}`). These should be a hard error, or gated behind `--strict` for one release cycle since it is a behavior change. Unrecognized single-brace content (e.g. bash `{ ...; }` groups, see F7) should keep the current informational path. This check is independent of whether the rule was supposed to expand (F1); singleton rules with no placeholders at all remain perfectly legal.

## F3 [MED] Failed-rule output rollback is invisible in the error message

When a rule fails, its already-created declared outputs are removed. That behavior is good (rerun cleanliness), but the failure message only lists *missing* outputs, e.g.:

```
WARN Shell exited 0 but declared outputs are missing — marking as failed
     rule=mixcr_test_<uuid> missing=.../TCGA-....mixcr.report.json
```

In our case a merge step had **successfully created** the other declared output (exit code 0, 1 of 2 declared outputs present), which then vanished — costing a debugging round to understand the merge had actually run. One clause in this message — "(declared outputs of failed rules are rolled back)" — would make the behavior self-documenting; naming the rolled-back files would be even better.

## F4 [MED] No native wait/chaining primitive for scripted supervision

For "run workflow A, wait for completion, then run validation", the available surface (`on_complete`/`on_error` hooks) covers in-workflow chaining, but there is no `oxo-flow run --wait/--detach` + `oxo-flow wait <run-id>` pair for external scripts/agents to block on completion and survive interrupted connections. We ended up writing a shell supervisor that polls for output files and exit markers — exactly the kind of boilerplate an agent-facing runner should own. A blocking `wait` with structured exit status (success / failed-count / timeout) would make lifecycle automation native. Documenting the `on_complete` cross-run chaining pattern would also help.

## F5 [LOW] cgroup memory-clamp warning is emitted per rule on large-RAM hosts

On a 1.4 TB RAM host with a 512 GB cgroup ceiling, every rule logs:

```
WARN cgroup memory limit is below the detected host memory — clamping the auto-detected ceiling ...
     host_mb=1417875 cgroup_limit_mb=524288 ram_only=true
```

The signal is useful once; per-rule repetition drowns real warnings. Emit once per run, or demote to INFO.

## F6 [LOW] Checkpoint file is per-directory, shared across workflow file variants

Running a `*.smoke.oxoflow` variant next to the main workflow in the same directory shares `./.oxo-flow/checkpoint.json`; the engine handles it correctly ("checkpoint belongs to a different workflow — ignoring its completed/failed records"), but a per-workflow checkpoint identity (workflow name + content hash, so file renames don't break resume; or alternatively a `--checkpoint-dir` flag) would make variant runs (smoke/test/full) first-class instead of silently non-reusable.

## F7 [DOCS] Shell `${var}` syntax coexists safely with `{placeholder}` — worth documenting

Our rule shells contain heavy bash (`${#parts[@]}`, `${ch##*.clones_}`, `{ ...; }` groups). The renderer left all `${...}` constructs intact and only flagged genuine single-brace idioms, with a helpful explanation:

```
... 'this shell idiom is allowed and the rule will run; it is flagged because
such content can also ride in on rendered wildcard values — verify it is intentional
```

This is the right behavior; an explicit docs note ("`${...}` is safe in shell blocks; only `{token}`-shaped single braces are placeholder candidates") would preempt user anxiety.

---

## Positives (for balance)

- **Pixi backend**: lockfile-aware env cache keying (`env=pixi:pixi.toml:<hash>:<hash>`), one-shot bootstrap, in-place verification of prebuilt envs — all worked exactly as designed.
- **Static gates**: `validate`/`dry-run`/`lint` caught real issues pre-run (manifest path resolution, deprecated `threads=`/`memory=` placement via W025).
- **Structured events**: the JSON event stream (`rule_started`/`rule_completed`/`rule_failed` with durations and exit codes) made log-based supervision reliable.
- **Provenance + report snapshots** (`--provenance`, `.oxo-flow/reports/report-*.json`) worked out of the box.
- **Graceful cross-workflow checkpoint handling** (F6) — correct behavior, just worth surfacing.

## Environment

- oxo-flow 0.30.0, built from source at `384db709` (== origin/main at the time)
- Cluster: 4-node Slurm cluster; runs executed on one node (64 cores, 1.4 TB RAM)
- Environment backend: pixi (mixcr 4.7.0, trust4 1.1.5, samtools 1.25, python 3.12)
- Workflow: 3 rules × 8 samples, `[workflow] metadata_file` + `[[sample_groups]]`
