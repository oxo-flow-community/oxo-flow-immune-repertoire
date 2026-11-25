# Repertoire Metrics Reference

> Purpose: define every column of `qc/sample_metrics.tsv`, how it is computed,
> why it is low or high, and how to use it downstream (covariates / filters).
> Applies to per-sample output of `scripts/sample_metrics.py` over the unified
> clonotype tables (MiXCR 4.7 + TRUST4 1.1.5 on GDC RNA-seq BAMs).
> Written 2026-11-15.

## 1. The two "lowness" questions (read this first)

**`mixcr_aligned_frac` ~ 0.01%–0.26% — why so low?**
This is NOT the BAM's overall alignment rate (that is STAR's 80–90%). It is
the fraction of reads that MiXCR could align to V/D/J/C **germline immune
templates**. In bulk tumor RNA-seq, receptor transcripts are inherently rare:
the transcriptome is dominated by housekeeping/ribosomal/tumor genes, and
lymphocyte receptor mRNA is ~1e-4 to 1e-3 of total. This matches the yields
reported in Bolotin 2017 (MiXCR RNA-seq) and Song 2021 (TRUST4). Low absolute
value is expected biology, not pipeline failure.

**`*_cdr3_frac` — even lower, why?**
It is the end of a three-stage funnel: total reads → immune-template-aligned
(aligned_frac) → reads yielding a complete CDR3 clonotype. Further losses
come from: alignments outside CDR3 (framework/UTR/intron transcription),
incomplete D-region or hypermutated templates that fail anchoring, and the
rna-seq preset's requirement that both CDR3 anchors (conserved C / F·W) be
present. BCR (IGH) suffers extra loss from somatic hypermutation.

**Interpretation rule**: absolute values are meaningless across datasets;
*relative* values across samples are the signal. `cdr3_frac` spanning
2e-6 … 1.2e-3 across our 8-sample test set tracks immune infiltration, not
technical quality.

## 2. Column dictionary

### Library size (covariates)

| column | source | meaning / use |
|---|---|---|
| `total_reads` | MiXCR `align.report` ("Total sequencing reads") | sequenced pairs entering align; primary sequencing-depth covariate |
| `mixcr_aligned_reads`, `_frac` | MiXCR align report | immune-template-aligned reads (frac stored as fraction of total) |
| `mixcr_cdr3_reads`, `trust4_cdr3_reads` | sum of clonotype read counts | repertoire captured by each tool (library size #2) |
| `bam_size_gb` | manifest | raw file size; coarse depth proxy |

### Immune content (the recommended filter)

| column | definition | use |
|---|---|---|
| `mixcr_cdr3_frac` | `mixcr_cdr3_reads / total_reads` | **immune-content covariate; primary filter**. Suggested floor: `> 1e-5` (≈ clonotypes ≥ 50 in our data) to call a sample "has repertoire". Cross-sample range in test set: 2e-6 (HNSC, "cold") → 1.2e-3 (OV, infiltrated) |
| `trust4_cdr3_frac` | same, TRUST4 AIRR `consensus_count` sum | TRUST4's de novo capture is template-independent and typically several-fold higher than MiXCR's — the difference itself quantifies template-dependence vs assembly sensitivity |

### Diversity / clonality (per tool)

| column | definition | notes |
|---|---|---|
| `clonotypes` | productive rows in unified table | richness (uncorrected; not Chao1 — rarefaction-level estimators deliberately deferred) |
| `clonality` | `1 − Shannon H / ln(S)` on clonotype read counts | 0 = polyclonal, 1 = oligoclonal; standard immuno-oncology metric |
| `gini` | Gini coefficient of clonotype read distribution | skew, complementary to clonality |
| `top_clone_frac` | largest clonotype's read share | AGENTS.md §7.1 flags > 0.5 as possible contamination/oligoclonality |

### Chain composition (per tool, read-share)

`frac_IGH / frac_IGK / frac_IGL / frac_TRA / frac_TRB / frac_TRD / frac_TRG`
— B vs T lineage infiltration and heavy/light balance. Examples from the
test set: OV 33% IGH (B-rich), HNSC 10% TRB (T-skewed, IGH low).

### Other QC / covariates

| column | meaning |
|---|---|
| `mean_cdr3_nt` | mean CDR3 nucleotide length (QC: framework contamination shortens it) |
| `sample_type` | barcode field 4 code: 01 primary, 02 recurrent, 03 blood-derived normal, 06 metastatic, 11 solid normal (see `manifests/sample_type_distribution.txt`) |
| `has_chimeric` | 1 = genomic+chimeric analyzed; 0 = genomic-only (642 samples: UCEC 380, COAD 190, READ 72) — detection-sensitivity covariate |
| `mixcr_generation` | provenance of the MiXCR table: `current` = produced by the rule definition now in the repository; `pre-cutover` = produced by the earlier local-mode workflow and carried across the cluster-backend cutover with a backfilled sentinel; `unknown` = sentinel absent. Both generations ran the identical MiXCR preset and filters (`rna-seq --species hsa --assemble-longest-contigs`, primary-alignment-only `-F 0x900`); only threading, memory, logging and scratch layout differ. Retain it as a stratification column in mixed analyses |

### 原始表的结构事实（2026-11-22 独立复核确认）

- **读数加权 ≠ 克隆计数加权**：按读数份额，**IGK 在 20/20 个样本中是第一链**；
  按克隆（行）数，IGH 在 13/20 中第一。链读数是表达量与组装效率的混合，不是
  细胞计数——跨样本比较链组成时须固定同一归一方式，并把 `total_reads` 作为
  协变量，不可把 `frac_*` 直接当作细胞比例。
- **out-of-frame 尾部**：MiXCR 合并表中典型 0.2–5% 的行（读加权通常 <1.2%）
  `aaSeqCDR3` 含 `_`/`*` 或 `cdr3_nt` 长度不是 3 的倍数（部分组装或移码）。
  2026-11-23 扩大到 38 个样本复核：全批 2.11%，范围 0.37–**9.14%**（离群样本
  `6efd0395`，IGK 单链 13.24%）。严格分析应按 `len(cdr3_nt) % 3 == 0` 过滤；
  读加权时影响可忽略，按行计数的丰富度指标会放大。
- **非生产性行会保留在统一表中**（复核 3 个 TRUST4 样本）：TRUST4 中
  `productive=F` 的行（`junction_aa` 为空、`junction` 非空）仍写入统一表，
  `cdr3_aa` 为空字符串（HU-A4GJ 4,852/120,923 行，BR-A4J5 405/8,935，
  24-1422 62/4,308）。这是有意的：核苷酸分析仍可用。`sample_metrics.py` 的
  `is_productive()`（非空、无 `*`/`_`、以 C 开头、以 F/W 结尾、且 `cdr3_nt`
  长度为 3 的倍数）已只统计生产性行，并保留 `*_raw` 计数；**任何下游氨基酸
  分析必须应用同一过滤**，不能直接按行计数。（2026-11-23 修正：原实现只排除
  `*`，会把 `CALPAL_SKLTF` 这类移码克隆算作生产性；回归测试
  `test/sample_metrics_selftest.py`。）
- **TRUST4 的 in-frame 判据**：`junction` 是核苷酸且密码子完整（含 IMGT 锚点），
  直接用 `len(junction) % 3 == 0`；**不要按字面剥离首 C / 末 W**——绝大多数行以
  TGT/TGC 密码子开头（BR-8369 中 3,842/3,940 行以 Cys 密码子开头，只有 8 行
  首字母恰为 `C`），字面剥离反而制造伪 out-of-frame。
- **MiXCR 的 `readCount` 是浮点字符串**（如 `1356.0`），不是整数；消费原始表时用
  `float()`/`int(float())`（`standardize_clonotypes.py::safe_int` 已处理），直接
  `int()` 会抛异常。已两次让独立复核者的自写解析器踩坑。
- **chimeric 克隆的链归属**：MiXCR 的 `clonalChainUsage` 会把跨链克隆记到某一链
  （本批 3 个：2 个 LUSC、1 个 STAD 样本，均为 TRBV/IGHV 混合），而合并规则把
  空链后缀记为 `ambiguous`（不猜测）。因此按统一表统计的链计数可能与 MiXCR
  报告差这几行，总数守恒；**链组成以统一表为准**。
- 另注：`readFraction` 是**每链内**归一（每链和 = 1），跨链相加会得到 2–8；
  `standardize_clonotypes.py` 已按 (sample, tool, chain) 重新计算 `clone_fraction`，
  但任何直接消费原始表的脚本都必须按链分层。`clone_fraction` 用 `%.6g` 格式化，
  格式化后的链内和会偏离 1：实测最大 **2.2e-6**（60 行的 TRB 链），而用精确分数
  重算的和**恰好 = 1**——纯属 6 位有效数字舍入。因此：严格校验应**从 `read_count`
  重算** `read_count/链内总和`（精确），不要对格式化字符串求和；若坚持对字符串
  求和，容差取 **5e-6**；项目自带 `verify_batch.py` 用 0.999–1.001（1e-3），不会
  被这种舍入误判。

## 2a. Authoritative thresholds (single source of truth, 2026-11-17)

> **This table is the only definition of "has repertoire" in the project.**
> AGENTS.md §7, METHODS.md and any downstream script must reference it rather
> than restating numbers — the old documents disagreed (1e-5 vs 0.1 %, which
> differ by 100x). `cdr3_frac` is **pooled across chains**, so it can pass on
> Ig reads while the sample has only a handful of TRB reads; T-cell statements
> need the chain-specific row.

| Purpose | Metric | Provisional threshold | Basis |
|---|---|---|---|
| "has BCR repertoire" | `mixcr_IGH_cdr3_reads` | ≥ 50 reads | **provisional** — calibrate against this cohort's detection-saturation curve |
| "has TCR repertoire" | `mixcr_TRB_cdr3_reads` | ≥ 20 reads | **provisional** — same |
| Effective depth for diversity | `mixcr_cdr3_reads` | ≥ 1,000 (rarefaction depth) | VDJtools requires equal depth before comparing diversity |
| Any cross-sample diversity comparison | `mixcr_clonality_rare` / `mixcr_gini_rare` | use the rarefied columns | raw richness/clonality track sequencing depth (4 orders of magnitude here) |
| Lineage-specific statements | `mixcr_tcr_*` / `mixcr_bcr_*` | never the pooled `mixcr_clonality` | Ig dominates the pool (>95 % of reads in most samples) |
| Covariates in every model | `total_reads`, `cdr3_frac`, `has_chimeric`, `read_layout`, `read_length`, `sample_type` | always include | Bolotin 2017 (read length), Thorsson 2018 (layout), review H1/H2 |

No literature defines a universal per-sample read floor for unselected bulk
RNA-seq; the rows marked *provisional* must be replaced by values derived from
this cohort's own saturation analysis before any publication claim.

## 3. Downstream usage guide

- **Filter "has repertoire"**: use §2a (chain-specific reads + rarefaction depth). The old pooled rule `mixcr_cdr3_frac > 1e-5 AND mixcr_clonotypes >= 50` is superseded: it passes samples whose TRB reads are in the single digits.
- **Regression covariates**: `total_reads` + `cdr3_frac` together (depth vs
  content are orthogonal — 221M reads with 2e-6 content is still "cold").
- **Do NOT** use raw `total_reads` alone as a size factor: repertoire yield
  is driven by infiltration, not depth.
- **Sensitivity flag**: analyses mixing `has_chimeric` 0/1 samples should
  include the flag as a covariate (or exclude the 642).
- **Comparing clonality across tools**: MiXCR and TRUST4 clonality agree
  directionally but differ in granularity (de novo assembly splits clones);
  use one tool consistently for cross-sample comparisons.

### 3a. TCR and BCR are not the same measurement (2026-11-19)

Measured across the first 10 end-to-end samples, `mixcr_TRB` read totals are
**22–322** while `mixcr_IGH` totals are **13–62,753** — TCRβ is 200–1000× rarer
because it is a signalling receptor at ~1–2 mRNA copies per T cell, whereas
Ig is the secreted effector of plasma cells, and tumour B cells are often
clonally expanded. Consequences, to be respected in any downstream design:

| Analysis | TCR (TRA/TRB/TRG/TRD) | BCR (IGH/IGK/IGL) |
|---|---|---|
| Per-sample diversity (clonality, Gini, top_clone_frac) | **No** — 17–227 distinct CDR3 with 1–2 reads each is not a repertoire estimate, and rarefaction cannot repair it | Yes, when the §2a floor is met |
| Presence / absence of repertoire | Yes | Yes |
| V–J usage frequencies | Yes, only **aggregated across samples** | Yes, per sample and aggregated |
| CDR3 length distribution | Yes, aggregated | Yes |
| Cross-tool Jaccard as agreement | **No** — the TRB universe (3–227) is so small that two random sets already reach J≈0.65; see `verify_batch`'s `chance_*` and `*_degenerate` columns | Yes, Ig J is 10–90× the chance baseline |
| Isotype / class-switch (`c_gene`) | n/a | Yes |

Calibration against the literature: Li, Li, Liu et al. (Nat Genet 2016,
doi:10.1038/ng.3581 — `sources.tsv:li2016ng`) recovered >600,000 CDR3 from
9,142 TCGA samples, i.e. **~66 per sample on average**. Our 3–227 per sample,
median in the low hundreds, sits in the same order of magnitude — the sparsity
is a property of the assay, not a defect of this pipeline.

**Rule of thumb**: report TCR results as population-level summaries with
per-sample depth shown; report BCR results per sample. Never put a per-sample
TCR clonality into a clinical regression.

## 4. Provenance & regeneration

```bash
cd $HOME/Projects/TCGA-mixcr
oxo-flow run qc_metrics.oxoflow -j 4 --provenance   # after tcga-mixcr.oxoflow
# → qc/sample_metrics.tsv (full cohort: 11,123 rows)
```
