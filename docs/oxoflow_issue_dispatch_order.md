# Dispatch ties are broken by instance *name*, which turns per-sample workflows into phase-ordered runs

**Environment**: oxo-flow v0.30.0 (`384db709`). Reproduced on the Slurm cluster
path and on the local path; this is core scheduler behaviour, not a
cluster-backend quirk.

## Observed

A per-sample pipeline with two independent producers and a join:

```
mixcr_S   (reads the sample's BAM)  ─┐
                                     ├─► standardize_S ─► cleanup_S
trust4_S  (reads the same BAM)      ─┘
```

2,780 samples → 11,120 instances. With equal priorities the run behaves as a
sequence of phases: **all** `mixcr` instances are dispatched before any
`trust4`, and a sample's `standardize` only starts after its own two
producers finish. Nothing is wrong numerically, but:

- the first `unified` output appears only after the slowest producer phase has
  drained (weeks at production scale),
- the two producer branches that are *fully independent* (they only read the
  same input) never overlap across the cohort,
- per-sample QC/progress is unavailable until the end of the run.

## Root cause (source)

`crates/oxo-flow-core/src/scheduler.rs:146-172`:

```rust
    // Sort by name for deterministic ordering
    ready.sort();
```

`crates/oxo-flow-core/src/backend/driver.rs:197-210` (cluster driver):

```rust
    /// local priority aging — ready rules are sorted by effective priority
    /// (declared + rounds waited, ties by name), so fresh high-priority …
    ready.sort_by(|a, b| {
        let pa = plan.rules.get(a).map(|sr| sr.rule.priority).unwrap_or(0) …
```

Instance names are `<rule>_<group>_<sample>`, so an alphabetical tie-break
groups every instance of one rule together: `cleanup_*` < `mixcr_*` <
`standardize_*` < `trust4_*`. The rule grouping is therefore an accident of
naming, not of dependencies.

The name sort itself is clearly deliberate (determinism); the problem is that
it is used as the *only* tie-break, and that it ignores the expansion key that
the engine already knows (the sample).

## Why this matters generically

Fan-out/fan-in per sample is the standard shape of bioinformatics workflows —
`fastqc ∥ trim`, `align ∥ qc`, two aligners feeding a merge, RNA + DNA
branches per sample. In all of them the current tie-break serialises the
branches across the cohort and defers every per-sample result to the end.

Users work around it by chunking the cohort (`--samples @sheet` per batch) or
by hand-tuning priorities, both of which trade one problem for another
(priority is global: raising a downstream rule's priority starves the
producer — observed: 20/20 slots on `trust4`, no `mixcr` for two hours).

## Proposal

Add a dispatch tie-break mode, off by default so existing runs are unchanged:

```toml
[engine]
# "rule" (default, current behaviour) | "sample"
dispatch_order = "sample"
```

With `"sample"`, ready rules sort by **(expansion key, rule order in the
workflow)** instead of by name — i.e. all ready stages of sample A are
dispatched before sample B's, and inside a sample the rule order from the
file decides. The engine already tracks the expansion key per instance, so no
name parsing is needed.

Effect on the pipeline above: `mixcr_S`/`trust4_S` start together, and as soon
as both finish, `standardize_S` and `cleanup_S` run — per-sample results from
the first completed sample onward, with the same total work and the same
`max_submitted` cap.

Optional refinement: keep priority as the primary key (unchanged), and apply
the sample-aware ordering only to the tie-break, so users who set priorities
keep full control.

## Current workaround

We drive each shard in chunks of 100 samples (`--samples @chunk.tsv`), so each
run is small enough that phase ordering costs hours rather than weeks. It
works, but the chunk size is now the knob that trades plan overhead against
time-to-first-result — a knob that would not be needed with a sample-aware
tie-break.
