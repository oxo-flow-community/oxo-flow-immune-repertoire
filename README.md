# flow-mixcr — TCGA RNA-seq immune repertoire pipeline (MiXCR + TRUST4)

> Part of the [oxo-flow-community](https://github.com/oxo-flow-community) pipeline catalogue.

A reproducible, oxo-flow–driven pipeline that extracts TCR/BCR clonotypes from
bulk RNA-seq BAMs (TCGA GDC harmonized) and produces per-sample, cross-tool
normalized clonotype tables plus repertoire metrics.

Built and operated by the Wang Shixiang Lab on a 4-node Slurm cluster; designed
to be reusable: any site with Slurm + a reference genome can run it.

## What it does

```
BAM (genomic [+ chimeric])
   │
   ├── mixcr        MiXCR rna-seq preset → per-chain clones merged into one table
   ├── trust4       TRUST4 de novo assembly → AIRR table
   ├── standardize  both tools → one 14-column unified schema (per-chain fractions)
   └── cleanup      integrity check + scratch reclamation
```

Per-sample outputs (`results/{cancer}/{barcode_12}/{sample_uuid}/`):

| file | content |
|---|---|
| `<barcode>.mixcr.clones.tsv` | MiXCR clones, all chains, `chain` column appended |
| `<barcode>.mixcr.qc.json` | MiXCR QC report (`[]` for zero-clonotype samples) |
| `<barcode>.trust4_airr.tsv` / `_cdr3.out` / `_report.tsv` | TRUST4 AIRR + native outputs |
| `<barcode>.unified.tsv` | 14-column unified clonotype table (MiXCR + TRUST4) |
| `<barcode>.*.done` | per-rule completion sentinels (see "Crash safety") |

Aggregate metrics: `oxo-flow run qc_metrics.oxoflow` → `qc/sample_metrics.tsv`
(~40 columns per sample; column dictionary in `docs/metrics_reference.md`).

## Requirements

- oxo-flow ≥ 0.30 (`oxo-flow --version`), Rust binary on `PATH`
- Slurm (`sbatch`/`squeue`/`scancel`), pixi for environment provisioning
- NCBI/IMGT references for TRUST4 (`reference/trust4/hg38_bcrtcr.fa`,
  `human_IMGT+C.fa` — fetch from https://github.com/liulab-dfci/TRUST4; not
  committed here for licence/size reasons)
- A MiXCR academic licence (`MI_LICENSE` in `~/.config/tcga-mixcr/env`)

## Repository layout

```
main.oxoflow          runnable pipeline (inline sample group + example metadata)
shard.oxoflow         production shape: manifests/metadata.tsv + groups.tsv,
                      one working directory per cluster node
qc_metrics.oxoflow    aggregate per-sample repertoire metrics
setup.oxoflow         environment provisioning + tool self-checks
pixi.toml/lock        pinned tool versions (mixcr 4.7.0, trust4 1.1.5, samtools)
scripts/              standardizer, comparators, metrics, Slurm driver,
                      cutover tooling, dashboard, node health probes
docs/                 methods, metrics reference, validation records, runbooks
manifests/            example manifests (test set + 2-sample shard example)
test/run.sh           acceptance gate: validate + lint + dry-run every workflow
metadata.json         catalogue entry (oxo-flow-community registry)
```

## Quick start

```bash
./test/run.sh                                  # validate + lint + dry-run
oxo-flow dry-run main.oxoflow                  # catalogue quickstart
oxo-flow run main.oxoflow -j 2                 # local mode (needs BAM inputs)
oxo-flow run shard.oxoflow --profile slurm     # cluster mode (per-rule jobs)
```

Production (4 shards, ≤20 CPUs per node):

```bash
bash scripts/cutover_cluster_backend.sh --apply   # deploy workflows + profiles
setsid nohup scripts/driver_shard.sh 1 >/dev/null 2>&1 &   # one driver per shard
python3 scripts/shards_status.py                  # dashboard
bash scripts/cluster_health_probe.sh              # temperature / CPU policy probe
```

## Crash safety (why the `.done` sentinels exist)

The pipeline runs for weeks. oxo-flow's cluster backend treats *output
existence* as completion, so a killed job that left a half-written file could
be skipped forever. Every rule therefore declares a **sentinel written as its
last action**; "all outputs exist" then means "the rule really finished".
The driver additionally cancels leftover jobs from a previous generation and
an output-count guard turns `rc=0` with missing samples into a retry. Details
and the upstream analysis: `docs/oxoflow_issue_cluster_backend.md`,
`docs/oxoflow_issue_cluster_status.md`.

## Validation record

- `docs/phase4_validation.md` — per-sample cross-method agreement (Jaccard /
  Spearman / Pearson) and the literature conformance audit.
- `docs/t1_cluster_validation.md` — cluster-backend end-to-end validation,
  including an independent read-only re-verification of every output.

## Licence / data

TCGA data are controlled-access; do not redistribute BAMs. Code is released
for academic reuse; MiXCR requires its own academic licence.
