# Methodological Basis for TCGA RNA-seq Immune Repertoire Analysis

> **Purpose**: Document the scientific rationale for tool selection, parameters, and validation strategies.
> **Date**: 2026-11-13

---

## 1. Background

### 1.1 Immune Repertoire from RNA-seq

Bulk tumor RNA-seq contains a mixture of:
- Tumor cell transcripts (~90-99%)
- Stromal cell transcripts
- **Infiltrating immune cell transcripts** (1-10%, including TCR/BCR mRNAs)

TCR/BCR mRNAs encode the hypervariable CDR3 region that defines clonotype specificity. Extracting these from RNA-seq requires:
1. Identifying reads containing V(D)J segments
2. Assembling full or partial CDR3 sequences
3. Correcting for sequencing/PCR errors

### 1.2 Two Complementary Approaches

| Aspect | MiXCR | TRUST4 |
|--------|-------|--------|
| **Core algorithm** | Alignment to reference V/D/J/C database → clonotype assembly | De novo assembly → realignment to IMGT reference |
| **Input** | FASTQ or BAM | BAM (preferred) or FASTQ |
| **Error correction** | Quality-guided mapping + fuzzy clustering | K-mer based, coverage threshold |
| **CDR3 extension** | Yes (extend step, TCR only) | Yes (assembly-based) |
| **Full-length** | assembleContigs (longest contig per clone) | Consensus contigs |
| **Best for** | High-confidence clonotypes, V/J assignment | De novo detection, novel alleles |
| **Citation** | Bolotin et al. 2017 Nat Biotechnol | Song et al. 2021 Nat Methods |

**Why run both?**
- MiXCR is more sensitive for low-abundance clonotypes (better error correction)
- TRUST4 can detect novel V/J alleles not in reference (de novo advantage)
- Cross-method agreement increases confidence in called clonotypes
- Discrepancies highlight challenging regions (e.g., hypermutated IGH)

---

## 2. MiXCR Methodology

### 2.1 Key Steps (from official documentation)

1. **Alignment** (`align`)
   - k-mer seed-and-vote algorithm
   - Two aligners: linear (TCR) and affine (BCR, handles indels)
   - Optional: paired-end merge, quality trimming

2. **Partial Assembly** (`assemblePartial`)
   - Rescues reads partially covering CDR3
   - Merges overlapping partial alignments from same molecule
   - Critical for short-read RNA-seq (~50bp fragments)

3. **CDR3 Extension** (`extend`)
   - Fills missing germline nucleotides at CDR3 edges
   - TCR only (BCR hypermutations make extension unreliable)
   - Uses reference germline sequences

4. **Clonotype Assembly** (`assemble`)
   - Groups alignments by CDR3 similarity
   - Two-layer error correction:
     - Quality-guided: rescue low-Q reads
     - Fuzzy clustering: correct PCR errors
   - Assembling feature: CDR3 (default for RNA-seq)

5. **Contig Assembly** (`assembleContigs`)
   - Reconstructs longest possible receptor sequence per clone
   - Alignment-guided consensus
   - Detects hypermutations (BCR)

### 2.2 RNA-seq Preset Parameters

From built-in `rna-seq` preset (v4.7.0):
```
inheritFrom: shotgun-base
pipeline:
  - align
  - assemblePartial
  - extend
  - assemble
  - assembleContigs
  - exportClones

align:
  trimmingQualityThreshold: 10
  trimmingWindowSize: 6

assemblePartial:
  repeats: 2

assemble:
  assemblingFeatures: [CDR3]
  clnaOutput: true  # preserve alignment-to-clone mapping
```

### 2.3 Required Mix-ins

- `--species hsa`: Homo sapiens reference
- `--assemble-longest-contigs` OR `--assemble-contigs-by VDJRegion`:
  - Longest: assemble maximum length contig per clone
  - VDJRegion: only clones with complete VDJ region

**Our choice**: `--assemble-longest-contigs` for maximum information recovery.

---

## 3. TRUST4 Methodology

### 3.1 Key Steps (from Song et al. 2021)

1. **Candidate Read Extraction**
   - From BAM: extract reads partially aligned to V/J/C genes
   - Uses BAM flags and alignment positions

2. **De Novo Assembly**
   - K-mer based assembly (default k=9)
   - Builds contigs from candidate reads
   - No reference required (novel allele detection)

3. **Realignment and Annotation**
   - Realign contigs to IMGT reference
   - Identify V, D, J, C genes and CDR3 boundaries

4. **Filtering**
   - Coverage threshold (default: contigMinCov=0)
   - CDR3 must be in-frame and contain conserved residues

### 3.2 Reference Files

- `hg38_bcrtcr.fa`: Genomic coordinates of V/D/J/C genes (for BAM extraction)
- `human_IMGT+C.fa`: IMGT reference + constant region sequences (for annotation)

Both included in TRUST4 distribution.

---

## 4. Input Data Specifics

### 4.1 GDC RNA-seq BAMs

TCGA RNA-seq BAMs from GDC are:
- **Aligned to**: GRCh38.d1.vd1 (GDC reference)
- **Sort order**: Coordinate-sorted
- **Read groups**: `@RG SM:` contains full 28-character TCGA barcode
- **Types**:
  - `*.rna_seq.genomic.gdc_realn.bam`: All reads (including unmapped)
  - `*.rna_seq.chimeric.gdc_realn.bam`: Chimeric reads (spanning splice junctions)

### 4.2 BAM → FASTQ for MiXCR

MiXCR can read BAM directly, but we convert to FASTQ for:
1. Merging genomic + chimeric (TRUST4 also benefits)
2. Consistent input format for both tools
3. Easier parallelization (FASTQ is splittable)

**Command**:
```bash
samtools fastq -1 R1.fq.gz -2 R2.fq.gz -0 /dev/null -s /dev/null -n \
  genomic.bam chimeric.bam
```

**Read naming**: `-n` keeps read names (no /1 /2 suffix) for MiXCR paired-end handling.

### 4.3 TRUST4 BAM Input

TRUST4 is designed for BAM input:
```bash
run-trust4 -b genomic.bam -f hg38_bcrtcr.fa --ref human_IMGT+C.fa
```

**Note**: TRUST4 internally extracts candidate reads. We do NOT need to merge chimeric BAM — TRUST4 handles both if given genomic BAM (chimeric reads are marked in BAM flags).

---

## 5. Validation Strategy

### 5.1 Internal Validation (Pipeline Correctness)

**Test samples**: 8 samples, one per disk, varying sizes.

**Checks**:
1. Barcode extraction matches expected format
2. FASTQ read counts match BAM read counts
3. MiXCR and TRUST4 both complete without errors
4. Output files non-empty and well-formed

### 5.2 Cross-Method Validation (MiXCR vs TRUST4)

**Per-sample comparison**:
- **CDR3 overlap**: Jaccard index of top-100 CDR3 amino acid sequences
  - Expected: > 0.5 (tools agree on abundant clonotypes)
  - Low overlap may indicate: low immune infiltration, hypermutated BCRs, or technical artifacts

- **V gene usage correlation**: Spearman ρ
  - Expected: > 0.7 (tools assign same V genes at similar frequencies)

- **Read count correlation**: Pearson r on log-transformed per-clone counts
  - Expected: > 0.7 (tools quantify clonotypes similarly)

**Interpretation**:
- High agreement: high-confidence result
- Low agreement: flag for manual review, check QC metrics

### 5.3 External Validation (vs Published TRUST4 Results)

**Zenodo dataset**: https://zenodo.org/records/14251977

- Same tool (TRUST4), same cancer types, but:
  - IDs are de-identified (e.g., `COAD_0`)
  - Amino acid CDR3 only (no nucleotide)
  - No abundance in reads (only relative abundance within sample)

**Comparison method**:
1. Compute cancer-type aggregate statistics from our TRUST4 results:
   - CDR3 length distribution
   - V gene usage frequencies
   - J gene usage frequencies
   - Chain type frequencies (TRA:TRB:IGH:IGK:IGL ratio)

2. Compute same statistics from Zenodo data

3. Compare:
   - Distributions: KS test or chi-square
   - Frequencies: Spearman correlation
   - Ratios: log-fold difference

**Expected**: High correlation (> 0.8) for V/J usage and chain type ratios. CDR3 length distributions should match closely.

### 5.4 Aggregate QC

**Per-sample metrics** (flag outliers):
- Total reads
- CDR3-containing reads (% of total)
- Clonotype count
- Top clone fraction (clonality)
- Tool agreement score (Jaccard)

**Per-cancer-type metrics**:
- Sample count
- Median clonotype count
- Median CDR3 read fraction
- Outlier samples (potential batch effects)

---

## 6. Error Sources and Mitigations

| Error Source | Impact | Mitigation |
|--------------|--------|------------|
| Low immune infiltration | Few CDR3 reads, unreliable clonotypes | Flag samples with < 0.1% CDR3 reads |
| RNA degradation | Short fragments, truncated CDR3 | assemblePartial rescues partial alignments |
| PCR errors | Artificial clonotype diversity | MiXCR fuzzy clustering; TRUST4 coverage threshold |
| Sequencing errors | Wrong bases in CDR3 | Quality-guided correction (both tools) |
| BCR hypermutations | V/J misassignment | TRUST4 de novo advantage; MiXCR affine aligner |
| Reference bias | Missed novel alleles | TRUST4 de novo assembly |
| Batch effects | Technical variation across TCGA batches | Include plate/batch in metadata; normalize in downstream analysis |

---

## 7. Output Standardization

### 7.1 Unified Schema

Both tools output different formats. Standardize to:

| Column | Type | Source |
|--------|------|--------|
| `barcode_28` | string | BAM @RG SM: |
| `barcode_15` | string | Derived (first 4 fields) |
| `barcode_12` | string | Derived (first 3 fields) |
| `cancer_type` | string | Directory name / barcode prefix |
| `tool` | string | "mixcr" or "trust4" |
| `chain` | string | TRA, TRB, IGH, IGK, IGL |
| `cdr3_nt` | string | Nucleotide sequence |
| `cdr3_aa` | string | Amino acid sequence |
| `v_gene` | string | IMGT gene name |
| `d_gene` | string | IMGT gene name (or empty) |
| `j_gene` | string | IMGT gene name |
| `c_gene` | string | Constant region (or empty) |
| `read_count` | int | Number of reads supporting this clone |
| `clone_fraction` | float | Fraction of total CDR3 reads |
| `clone_count` | int | Absolute clonotype count (if applicable) |

### 7.2 MiXCR to Unified

MiXCR `.clones.tsv` columns (selected):
- `cloneCount` → `clone_count`
- `cloneFraction` → `clone_fraction`
- `nSeqCDR3` → `cdr3_nt`
- `aaSeqCDR3` → `cdr3_aa`
- `vGene` → `v_gene`
- `jGene` → `j_gene`
- `cGene` → `c_gene` (if present)

### 7.3 TRUST4 to Unified

TRUST4 `trust_cdr3.out` columns:
- Column 1: consensus_id (ignore)
- Column 2: chain_type → `chain` (map: `TRA`, `TRB`, `IGH`, `IGK`, `IGL`)
- Column 3: CDR3_AminoAcid → `cdr3_aa`
- Column 4: VGene → `v_gene`
- Column 5: JGene → `j_gene`
- Column 6: ConstantGene → `c_gene`
- Column 7: Abundance → `read_count` (note: TRUST4 abundance is relative weight, not raw reads)

**Note**: TRUST4 does not output nucleotide CDR3 directly in `trust_cdr3.out`. Extract from `trust_annot.fa` or re-translate from `trust_final.out`.

---

## 8. References

1. **MiXCR**: Bolotin DA, et al. "Antigen receptor repertoire profiling from RNA-seq data." *Nat Biotechnol* 35, 908–911 (2017). https://doi.org/10.1038/nbt.3979

2. **MiXCR methods**: Bolotin DA, et al. "MiXCR: software for comprehensive adaptive immunity profiling." *Nat Methods* 12, 380–381 (2015). https://doi.org/10.1038/nmeth.3364

3. **TRUST4**: Song L, et al. "TRUST4: immune repertoire reconstruction from bulk and single-cell RNA-seq data." *Nat Methods* 18, 627–630 (2021). https://doi.org/10.1038/s41592-021-01142-2

4. **TRUST4 TCGA**: Zenodo record 14251977, https://doi.org/10.5281/zenodo.14251977

5. **GDC**: "GDC Reference Files," https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files

6. **TCGA Barcode**: "TCGA Barcode," https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/

---

*End of METHODS.md*
