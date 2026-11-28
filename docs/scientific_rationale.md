# 方法学依据汇总：每个选择背后的文献与工具文档

> 目的：把流程里每一处会影响科学结论的选择，对应到**官方文档或同行评议文献**，
> 并标明我们自己是否做过本地验证、以及有没有被本地实验推翻的说法。
> 状态含义：**文献支持** = 有文档/论文明确依据；**本地验证** = 我们自己的实测；
> **经验选择** = 无公认标准，属本项目的工程/统计决策（下游必须当协变量或做敏感性分析）。
> 维护规则：任何方法学改动都必须同步更新本表，否则视为文档缺陷。

---

## 1. 工具与版本

| 组件 | 版本 | 依据 |
|---|---|---|
| MiXCR | 4.7.0 | [MiXCR 官方文档](https://mixcr.com/mixcr/)；软件论文 Bolotin 2015 (Nat Methods 12:380–381)；RNA-seq 方法论文 Bolotin 2017 (Nat Biotechnol 35:908–911) |
| TRUST4 | 1.1.5 | Song 2021 (Nat Methods 18:627–630)；[liulab-dfci/TRUST4 README](https://github.com/liulab-dfci/TRUST4) |
| 参考 FASTA | `hg38_bcrtcr.fa` (md5 `236a217d…`), `human_IMGT+C.fa` (md5 `5511c50f…`) | TRUST4 README 指定；与上游 master commit `286335cc…` 一致（记录见 `reference/trust4/download.log` 与 `reference/trust4/PINNED.txt`）。**转载/复现时必须用这两个 md5**，上游 master 变动会导致结果不可复现 |
| oxo-flow | v0.31.0 | 项目引擎，工作流格式见 `docs/guide`（oxo-flow-src） |

## 2. 输入与预处理

| 选择 | 依据 | 本地验证 | 状态 |
|---|---|---|---|
| 直接读 GDC harmonized BAM（不转 FASTQ） | GDC BAM 不带 READ1/READ2 SAM flag（0x40/0x80）；`samtools fastq -1/-2` 会把所有 reads 丢进 single/other 流 | 2026-11-14 实测：转换后 R1/R2 为空 | 文献支持 + 本地验证 |
| **只用一个 BAM：genomic** | GDC 把 chimeric alignments 同时保留在 genomic BAM 中（chimeric BAM 是其子集）；两文件都喂给 MiXCR 会重复计数 | 2026-11-17 实测样本 `44f9e45d`：chimeric 前 5,000 条读段的 4,839 个唯一 read name 在 genomic 命中 9,741 条记录（≈2/name，双端 mate 都在）；MiXCR 报告 `totalReadsProcessed` 精确等于 `(genomic主比对 + chimeric主比对)/2` | **本地验证（推翻了早期"chimeric 含 genomic 没有的 V-J 读段"的假设）** |
| `samtools view -F 0x900`（去 secondary/supplementary） | SAM/BAM 规范：0x100=secondary、0x800=supplementary 是同一条物理读段的额外比对记录，保留会虚增克隆计数；Bolotin 2017 的 RNA-seq 流程同样只用主比对 | MiXCR 内部是否也过滤未验证，因此显式过滤 | 文献支持 + 保守做法 |

## 3. MiXCR 分支

| 选择 | 依据 | 本地验证 | 状态 |
|---|---|---|---|
| `mixcr analyze rna-seq --species hsa` | MiXCR 官方 presets：`rna-seq` 专为“Non-enriched fragmented (shotgun) RNA-Seq data”设计，正是 TCGA bulk RNA-seq 的场景 | validate/lint + 实跑 | 文献支持 |
| `--assemble-longest-contigs` | **MiXCR 文档明确写：`rna-seq` 预设必须提供 `--species` 以及 `--assemble-contigs-by <feature>` 或 `--assemble-longest-contigs`** —— 这不是可选项 | A/B 实测：对已有 `.clna` 分别用默认与 `--assemble-longest-contigs` 跑 assembleContigs，逐链 `exportClones` 的 md5 **完全相同**；`mixcr analyze` 生成的 assembleContigs 命令行本身也不含该 flag | 文档要求 + 本地验证（**对导出克隆表无影响；不要用它宣称"最大化 contig 信息"**） |
| `-Xmx40g` 作为**参数**而非 `JAVA_TOOL_OPTIONS` | MiXCR 启动脚本在参数里没有 `-Xmx` 时会自行追加 `min(12g, 内存余量)`，命令行 `-Xmx` 会覆盖环境变量 | 2026-11-17 实测：`JAVA_TOOL_OPTIONS=-Xmx16m` + 参数 `-Xmx2g` 可以正常运行 | 本地验证（源码 + 实验） |
| 每链导出再合并（加 `chain` 列） | `rna-seq` 预设按链分别导出 `<prefix>.clones_<CHAIN>.tsv`，链集合随样本变化，无法静态声明为输出 | 实跑 39 个样本 | 工具行为 |
| 空 chain 后缀 → `ambiguous` | 导出文件名以 `.clones_.tsv` 结尾表示克隆无法归入单一 locus（嵌合/跨位点） | 实测 `TCGA-63-7020`：V=TRAV22、J=TRBJ2-6，2 reads，早期被误归 TRA | 本地验证 |

## 4. TRUST4 分支

| 选择 | 依据 | 本地验证 | 状态 |
|---|---|---|---|
| 不用 `--repseq` | README：`--repseq` 仅用于 "bulk, **non-UMI-based TCR-seq or BCR-seq**"，即**靶向**组库数据；TCGA 是未富集 bulk RNA-seq | — | 文献支持（早期文档曾写反，已更正） |
| 只传一个 `-b` | README 的 `-b STRING` 是单值；源码 `run-trust4` 只取 `$bamFiles[0]` | 实测：传两个 `-b` 时 `bam-extractor` 命令只含第一个 | 文档 + 本地验证（**第二个 `-b` 一直被静默忽略**） |
| `--abnormalUnmapFlag` | README：处理 "the flag in the BAM for the unmapped read-pair is nonconcordant"；GDC BAM 的 flag 约定不标准 | 2026-11-14 实测：不加该参数，大样本必现 "Two reads from the unaligned fragment are not showing up together" | 文档 + 本地验证 |
| `--clean 1` | README：清理中间文件 | — | 文档 |
| 解析 `*_airr.tsv` 而非 `trust_cdr3.out` | AIRR Community schema（`locus`/`junction`/`junction_aa`/`consensus_count` 有明确定义）；`trust_cdr3.out` 的列数随版本变化 | 实测 v1.1.5 的 `cdr3.out` 是 13 列且带前导 contig id，按位置解析会错位 | 标准 + 本地验证 |

## 5. 统一表与指标

| 定义 | 依据 | 状态 |
|---|---|---|
| `chain` 由 V 基因前缀推断（IMGT 命名 `TRAV/TRBV/TRDV/IGHV/IGKV/IGLV`），工具自带 locus 优先 | IMGT 基因命名；TRUST4 的 `locus` 字段遵循 AIRR | 标准 |
| `clone_fraction` 按 **(tool, chain)** 归一，不跨链合并 | BCR 与 TCR 来自不同细胞群，跨链分母会混淆两个受体系统；MiXCR 自身的 readFraction 也是同分母 | 经验选择（有生物学依据） |
| 两个工具的 CDR3 都含 IMGT Cys/Phe-Trp 锚点 | 实测：MiXCR 的 `CATRGGRSSEGHANW` 对应 `TGTGCCACT…AACTGG`（起 Cys、止 Trp）；TRUST4 的 `junction` 同样含锚点。对 **productive** 行（当前快照 4,564,170 行）做 nt→aa 翻译校验，**100% 一致**；341,407 个非 productive 行含 `_`/`*`，约 92% 可译（预期，不参与生产性分析） | 标准（**原始 aa 可直接匹配**，`verify_batch` 的跨工具 Jaccard 正是如此）。若下游要剥离锚点，必须两侧同时剥离，不能只截一侧；旧版本文档写的“MiXCR 不含锚点、不可直接匹配”是错误陈述，已更正 |
| clonality = 1 − Shannon/ln S、Gini、top clone fraction | 免疫肿瘤学常用指标；Thorsson 2018 用同类多样性指标描述 TCGA 免疫景观；VDJtools 提供标准实现 | 文献支持（**但必须做深度校正，见下**） |
| 多样性/丰富度比较前需 downsampling/rarefaction | VDJtools 论文与官方教程：比较 richness/Shannon/clonality 前先 `DownSample` 到共同深度；Jaccard 是 presence/absence 且被浅样本主导 | 文献支持（**当前实现未做 —— 列为待修项 H6**） |
| TCR 与 BCR 指标必须分开算 | Thorsson 2018 的 TCR（MiTCR）与 BCR（VDJer）多样性是分开计算的；混合池化会让 TCR 结论被 B 细胞稀释 | 文献支持（**当前实现池化 —— 待修项 H7**） |
| read length / layout 必须作为协变量 | Bolotin 2017 指出 RNA-seq 组库提取强烈依赖读长与单双端；Thorsson 2018 方法学明确“优先 76bp，否则 50bp” | 文献支持（**当前未记录 —— 待修项 H2**） |

## 6. QC 阈值的依据（目前**没有**公认标准）

| 阈值 | 现状 | 依据 | 结论 |
|---|---|---|---|
| `mixcr_cdr3_frac > 1e-5` | 项目自定 | **无文献来源**；Bolotin 2017 报告的 RNA-seq 组库读段比例量级为 TCR ~1e-5–1e-7、Ig ~1e-4–1e-7，只能说明"低是正常的" | **经验选择**：必须按本队列分布校准，不能称为"TCGA 标准" |
| `clonotypes ≥ 50` | 项目自定 | 无文献 | 经验选择 |
| `cdr3 reads > 0.1%`（AGENTS §7.1） | 与上面冲突 | 来自靶向测序直觉，不适用于 unselected bulk RNA-seq | **应删除** |
| 链特异阈值（IGH reads ≥ N / TRB reads ≥ N） | 待定 | 近期 TCGA-ACC TRUST4 再分析（Cancer Med 2025）以"能否检出"为筛选，未给统一阈值 | 待用本队列分布 + 重复样本验证后确定 |

## 7. 被本地实验修正/推翻的文献或早期假设

1. **"chimeric BAM 含 genomic 没有的 V-J 读段，必须合并"** → 错。实测为子集，合并导致重复计数（C1）。
2. **"TRUST4 多个 `-b` 会合并输入"**（AGENTS.md §3.3 旧表述）→ 错。只使用第一个，第二个静默忽略。
3. **"`--assemble-longest-contigs` 能最大化 contig 信息"** → 对**导出的克隆表**无任何可观测影响（md5 相同）；它是 preset 的必需参数，不是增益开关。
4. **"`JAVA_TOOL_OPTIONS=-Xmx40g` 生效"** → 错。启动器会覆盖，必须写成参数。

## 7a. X. Shirley Liu 实验室的工作作为直接方法学参照

TRUST4 本身就是 **X. Shirley Liu 实验室**的软件（Song L, Cohen D, Ouyang Z, Cao Y,
Hu X, **X. Shirley Liu**，Nat Methods 2021；Crossref 已核对作者顺序）——所以
本流程的 TRUST4 分支不是"参考了某篇文献"，而是**直接使用该实验室的工具并遵循其
README 的参数语义**（单 `-b`、`--abnormalUnmapFlag`、bulk RNA-seq 不用 `--repseq`）。

同一实验室更早的 TCGA 组库工作提供了**方法学与可比性基准**：

| 工作 | 与本项目的关系 |
|---|---|
| **Li B, Li T, Liu JS, Liu XS. Landscape of tumor-infiltrating T cell repertoire of human cancers. *Nat Genet* 2016;48:725–732** | 首个大规模 TCGA RNA-seq TCR 组库全景（**9,142 样本 / 29 癌种 / >60 万 CDR3**）。他们报告 CDR3 检测率、长度分布、V 基因使用与健康外周血的异同，以及多样性–突变负荷关联——这些正是我们要复现/对照的量。**注意**：他们用的是自研 CDR3 推断方法，不是 TRUST4/MiXCR，因此只能做**趋势/聚合层面**对照，不能做样本级数值比较。 |
| **TRUST4 (Song 2021)** | 我们 TRUST4 分支的实现在此；其 README 对 `-b`/`--repseq` 的语义即我们的依据。 |
| **Zenodo record 14251977（TRUST4 pan-TCGA 去标识化报告，1.2 GB）** | 项目开始时由用户提供（AGENTS.md §4.4）。用于**聚合层面对照**：链比例 Pearson **r = 0.898**（30 癌种 × 7 链）、CDR3 aa 长度分布几乎重合（我们 11aa 36.8% vs 32.3%，12aa 15.8% vs 16.5%）。因原始 ID 已去标识化（`ACC_0` 这类），**不能做样本级匹配**，只能做汇总分布检验。 |

**从 Liu 实验室工作直接吸收的注意点**

1. **检测率与深度是首要混杂**：Li 2016 在 9,142 个样本上报告的是"能检出 CDR3 的
   样本比例"与长度/基因使用分布，而不是单样本绝对克隆数——与我们的 H2/H6 结论
   （读长/深度必须做协变量、多样性需 rarefaction）一致。
2. **脑与肾肿瘤的组库特征明显不同**（Li 2016 明确指出）→ 我们做癌种比较时，
   LGG/GBM/KIRC/KIRP 需单独看待，不能仅按"泛癌"解释。
3. **"TCR 多样性 vs 突变负荷"**是他们的核心发现之一 → 我们下游若做类似关联，
   必须先按 H6 做深度校正，否则多样性差异可能只是测序深度差异。

**明确不做的事**：Li 2016 的新抗原预测（SPAG5/TSSK6、PRAMEF4 突变）与 HLA 结合
预测不在本项目范围内。

## 7b. 参数级溯源：我们传的每个参数与"没传"的默认值

完整的展开式参数已固化在仓库里，**不是靠记忆或转述**：

- `docs/parameters/mixcr_rna-seq_v4.7.0.yaml` —— `mixcr exportPreset --preset-name
  rna-seq --species hsa --assemble-longest-contigs` 的**完整展开**（9.5 KB）：
  pipeline = align → assemblePartial → extend → assemble → assembleContigs →
  exportClones，以及每一步的全部打分/阈值参数。MiXCR 4.7.0 下导出，可逐字段复核。
- `docs/parameters/trust4_v1.1.5_usage.txt` —— `run-trust4` 的原样 usage 文本
  （v1.1.5-r573），列出每个选项的语义与**默认值**。

| 我们传的参数 | 文档依据 | 我们**未传**、因此取默认值的关键项 |
|---|---|---|
| `mixcr analyze rna-seq` | MiXCR presets 文档：`rna-seq` 面向 "Non-enriched fragmented (shotgun) RNA-Seq data" | 预设内部全部打分/阈值（见 YAML） |
| `--species hsa` | 同上，preset 必填 | — |
| `--assemble-longest-contigs` | 同上，**必填**（二选一，另一为 `--assemble-contigs-by`） | — |
| `--force-overwrite` | MiXCR 4.x 默认拒绝覆盖已有输出；批量重跑必须显式声明 | — |
| `-t 4` | `mixcr analyze` 通用参数 | — |
| `-Xmx40g`（**作为参数**） | MiXCR 启动脚本行为（源码 + 本地实验） | — |
| `run-trust4 -b` | README/usage：`-b STRING`（单值） | — |
| `-f hg38_bcrtcr.fa` / `--ref human_IMGT+C.fa` | usage：`-f` 坐标+序列；`--ref` IMGT 注释（recommended） | — |
| `-o` / `-t 4` | usage | — |
| `--clean 1` | usage：0 不清理 / **1 清中间文件** / 2 只留 AIRR | — |
| `--abnormalUnmapFlag` | usage 原文；GDC BAM flag 约定不标准（本地实测） | — |
| **不传** `--repseq` | usage：仅用于 "bulk, non-UMI-based TCR-seq or BCR-seq"（靶向），我们是未富集 RNA-seq | 默认关闭 ✅ 正确 |
| **不传** `-k` | — | k-mer = 9（默认） |
| **不传** `--contigMinCov` / `--minHitLen` / `--skipMateExtension` 等 | — | 全部取默认（见 usage 文本） |

**结论**：**每一条显式传入的参数都能在官方文档或工具 usage 中找到对应条目**；
我们**没有**自定义任何算法内部参数——全部沿用作者默认（展开式 YAML 与 usage
文本已入库，可随时比对）。若日后要改内部参数，必须先在此表登记并说明依据。

## 8. 可直接引用的参考文献

- Bolotin DA, et al. **Antigen receptor repertoire profiling from RNA-seq data.** *Nat Biotechnol.* 2017;35(10):908–911. doi:[10.1038/nbt.3979](https://doi.org/10.1038/nbt.3979)
- Li B, Li T, Liu JS, **Liu XS**. **Landscape of tumor-infiltrating T cell repertoire of human cancers.** *Nat Genet.* 2016;48(7):725–732. doi:[10.1038/ng.3581](https://doi.org/10.1038/ng.3581)
- Bolotin DA, et al. **MiXCR: software for comprehensive adaptive immunity profiling.** *Nat Methods.* 2015;12(5):380–381. doi:[10.1038/nmeth.3364](https://doi.org/10.1038/nmeth.3364)
- Song L, et al. **TRUST4: immune repertoire reconstruction from bulk and single-cell RNA-seq data.** *Nat Methods.* 2021;18(6):627–630. doi:[10.1038/s41592-021-01142-2](https://doi.org/10.1038/s41592-021-01142-2)
- Thorsson V, et al. **The Immune Landscape of Cancer.** *Immunity.* 2018;48(4):812–830.e14. doi:[10.1016/j.immuni.2018.03.023](https://doi.org/10.1016/j.immuni.2018.03.023)
- Shugay M, et al. **VDJtools: Unifying Post-analysis of T Cell Receptor Repertoires.** *PLoS Comput Biol.* 2015;11(11):e1004503. doi:[10.1371/journal.pcbi.1004503](https://doi.org/10.1371/journal.pcbi.1004503)
- MiXCR 官方 presets 文档：<https://mixcr.com/mixcr/reference/overview-built-in-presets/>
- TRUST4 README：<https://github.com/liulab-dfci/TRUST4>
- AIRR Community schema：<https://docs.airr.io/>
- GDC 参考文件与 BAM 处理说明：<https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files>
- Zenodo TRUST4 TCGA（聚合对照）：<https://zenodo.org/records/14251977>

## 9. 仍需补文献支撑的项（诚实列出）

1. **每链最低读段数阈值**：没有公认值，必须用本队列的检测饱和度曲线/重复样本确定。
2. **rarefaction 的目标深度**：取决于下游问题，需在分析计划中写明。
3. **肿瘤纯度 / 白细胞分数校正**：TCGA 免疫景观提供了 leukocyte fraction，但与本数据的整合方式尚无标准做法。
4. **MiXCR vs TRUST4 的绝对定量一致性**：文献一致认为两者粒度不同（de novo 组装 vs 比对组装），只能做相对/秩相关比较，不能假设绝对一致。
