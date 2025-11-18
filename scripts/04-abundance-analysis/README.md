# Abundance Analysis Scripts

Scripts for mapping reads to assembly, differential abundance analysis, and taxonomic annotation.

## Prerequisites

- minimap2 (v2.30+)
- samtools (v1.19+)
- R (v4.0+) with packages: edgeR, tidyverse, ggplot2
- Python 3.7+ with pandas

## Workflow Overview
```
Assembled contigs + Filtered reads
    ↓
[1] map_reads_to_assembly.sh → BAM files (per sample)
    ↓
[2] build_count_matrix.sh → Count matrix (samples × contigs)
    ↓
[3] run_edger_analysis.R → Differential abundance results
    ↓
[4] empirical_fdr_calibration.R → Calibrated p-value thresholds
    ↓
[5] annotate_contigs.py → Taxonomic annotation
```

---

## Configuration

Replace placeholders in scripts:

- `<PROJECT>`: NCI project code
- `<MINIMAP2_PATH>`: minimap2 binary directory
- `<INPUT_CONTIGS_FASTA>`: Assembled/polished contigs
- `<INPUT_FILTERED_FASTQ_DIR>`: Quality-filtered reads from step 01
- `<OUTPUT_BAM_DIR>`: Directory for BAM files
- `<INPUT_METADATA_CSV>`: Sample metadata (barcode, treatment, etc.)
- `<SMAG_TAXONOMY_CSV>`: SMAG taxonomy database

---

## Step 1: Map Reads to Assembly

Map filtered reads back to assembled contigs using minimap2.
```bash
# Edit script
nano map_reads_to_assembly.sh

# Key parameters:
# CONTIGS_FILE: Assembled contigs (from step 03)
# READS_DIR: Filtered reads (from step 01)
# OUTPUT_DIR: Output directory for BAM files

# Submit job
qsub map_reads_to_assembly.sh
```

**Output**: Sorted and indexed BAM files (one per sample)

**Resource requirements**:
- Queue: hugemem
- Memory: 500 GB
- Walltime: 24 hours
- Uses minimap2 preset `-ax map-ont` for Nanopore reads

---

## Step 2: Build Count Matrix

Extract read mappings from BAM files and generate count matrix.
```bash
# Edit script
nano build_count_matrix.sh

# Set paths:
# BAM_DIR: Output from step 1
# OUTPUT_CSV: Mapped reads table
# MATRIX_OUTPUT: Count matrix (samples × contigs)

# Key parameter:
# MIN_MAPQ: 60 (high-quality mappings only)

# Submit job
qsub build_count_matrix.sh
```

**Outputs**:
- `mapped_reads_table.csv`: Individual read mappings
- `count_matrix.csv`: Aggregated counts (samples × contigs)

**Filtering**:
- MAPQ ≥ 60 (high-quality alignments)
- Excludes unmapped, secondary, and supplementary alignments

---

## Step 3: Differential Abundance Analysis

Identify differentially abundant contigs using edgeR.
```bash
# Edit script
nano run_edger_analysis.R

# Set paths:
# count_matrix_file: Output from step 2
# metadata_file: Sample metadata
# output_dir: Results directory

# Run
Rscript run_edger_analysis.R
```

**Outputs**:
- `edger_results_categorical.csv`: Results for categorical comparison
- `edger_results_continuous.csv`: Results for continuous covariate
- `edger_*_fit.rds`: Saved edgeR objects for downstream analysis

**Analysis modes**:
1. **Categorical**: Compare groups (e.g., small vs. large biomass)
2. **Continuous**: Test correlation with continuous variable (e.g., shoot mass)

---

## Step 4: Empirical FDR Calibration

Calibrate FDR thresholds using permutation testing and generate comprehensive visualizations.
```bash
# Edit script
nano empirical_fdr_calibration.R

# IMPORTANT: Set paths to plotting functions
# These should point to the 05-visualization scripts
# PLOT_VOLCANO_R: ../05-visualization/plot_volcano.R
# PLOT_HEATMAP_R: ../05-visualization/plot_heatmap.R

# Set paths (same as step 3)

# Key parameter:
# B: 500 (number of permutations; increase to 1000 for final analysis)

# Run
Rscript empirical_fdr_calibration.R
```

**Outputs**:
1. **FDR Calibration**:
   - `empirical_fdr_curve.png`: FDR vs. p-value threshold
   - `observed_vs_null_discoveries.png`: Observed vs. null discoveries
   - `empirical_tau.csv`: Calibrated p-value threshold

2. **Volcano Plots**:
   - `volcano_empirical.png`: Using empirical τ threshold
   - `volcano_bh.png`: Using Benjamini-Hochberg FDR

3. **Heatmaps**:
   - `heatmap_empirical.png`: Top 50 contigs (empirical threshold)
   - `heatmap_bh.png`: Top 50 contigs (BH FDR)

4. **Results Tables**:
   - `DA_contigs_empirical_fdr.csv`: Significant contigs
   - `DA_contigs_list.txt`: List of significant contig IDs
   - `all_contigs_results.csv`: Complete results table

**Dependencies**:
- Requires `plot_volcano.R` and `plot_heatmap.R` from `05-visualization/`
- Make sure these are sourced correctly at the top of the script

**Why multiple visualization modes?**
- **Empirical (τ)**: Conservative, permutation-based, controls false discoveries empirically
- **BH FDR**: Standard, assumes theoretical null distribution, commonly used in publications
- **Comparison**: Shows difference between theoretical and empirical approaches

---

## Step 5: Annotate Contigs

Assign taxonomy to contigs based on MAG mapping.
```bash
# Edit script
nano annotate_contigs.py

# Set paths:
# MAPPED_CONTIGS_CSV: Contig-to-MAG mapping
# SMAG_TAXONOMY_CSV: SMAG taxonomy database
# PLANT_MAG_FILES: Optional plant genome MAG lists

# Run
python annotate_contigs.py
```

**Outputs**:
- `annotated_contigs.csv`: Contigs with full taxonomy
- `taxonomy_summary_phylum.csv`: Counts by phylum
- `taxonomy_summary_family.csv`: Counts by family

**Taxonomy assignment**:
1. Match contigs to SMAGs → assign microbial taxonomy
2. Match contigs to plant genomes → assign host taxonomy
3. Clean and standardize taxonomy labels

---

## Expected Results

### Mapping Statistics
- **Mapping rate**: 70-95% (depends on assembly completeness)
- **High-quality mappings (MAPQ ≥ 60)**: 50-80%

### Differential Abundance
- **Significant contigs (FDR < 0.05)**: 10-100 (depends on treatment effect)
- **Typical log fold changes**: 0.5-3.0

### Taxonomy
- **Phyla**: Typically 5-15 major phyla
- **Most abundant**: Proteobacteria, Bacteroidota, Actinobacteriota
- **Plant-derived contigs**: 1-10% (if plant genomes included)

---

## Troubleshooting

### Low Mapping Rate

**Symptoms**: <50% of reads map to assembly

**Possible causes**:
- Assembly incomplete
- Different samples used for assembly vs. mapping
- Contamination

**Solutions**:
1. Check assembly completeness (N50, total length)
2. Verify sample metadata
3. Run QC on both assembly and reads

---

### No Significant Contigs

**Symptoms**: All FDR > 0.05

**Possible causes**:
- Small treatment effect
- High biological variability
- Insufficient sample size
- Low coverage

**Solutions**:
1. Check power analysis (need n ≥ 5 per group)
2. Filter for high-coverage contigs (≥10× coverage)
3. Consider alternative analysis methods
4. Increase sample size if possible

---

### Permutations Too Slow

**Symptoms**: Empirical FDR calibration takes >24 hours

**Solutions**:
1. Reduce number of permutations (B = 100 for testing)
2. Filter contigs before analysis (keep high coverage only)
3. Use parallel processing (if available)

---

## Integration with Downstream Analysis

After abundance analysis:
```bash
# Abundance complete
qsub scripts/04-abundance-analysis/map_reads_to_assembly.sh
qsub scripts/04-abundance-analysis/build_count_matrix.sh
Rscript scripts/04-abundance-analysis/run_edger_analysis.R
Rscript scripts/04-abundance-analysis/empirical_fdr_calibration.R
python scripts/04-abundance-analysis/annotate_contigs.py

# Next: Binning and functional annotation
cd ../05-binning
# Or: Functional analysis of DA contigs
cd ../06-functional-annotation
```

---

## Parameter Tuning

### Mapping Quality Threshold (MAPQ)

| MAPQ | Interpretation | Use Case |
|------|---------------|----------|
| 0    | Any alignment | Exploratory |
| 20   | ≥99% correct  | Standard |
| 60   | ≥99.9999% correct | Conservative (recommended) |

### edgeR Filtering
```R
# Default: filterByExpr (adaptive)
# Manual:
keep <- rowSums(cpm(y) > 1) >= 3  # CPM > 1 in ≥3 samples
```

### Log Fold Change Threshold
```R
# Conservative (recommended)
abs(logFC) >= 0.5  # 1.4-fold change

# Moderate
abs(logFC) >= 1.0  # 2-fold change

# Liberal
abs(logFC) >= 0    # Any change
```

---

## References

- edgeR: https://bioconductor.org/packages/edgeR/
- minimap2: https://github.com/lh3/minimap2
- Permutation FDR: Storey & Tibshirani (2003) PNAS