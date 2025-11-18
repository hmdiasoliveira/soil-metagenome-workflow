# Visualization Scripts

Scripts for generating publication-quality plots of assembly statistics, differential abundance results, and taxonomic composition.

## Prerequisites

- R (v4.0+) with packages:
  - tidyverse, ggplot2, plotly
  - vegan, compositions, zCompositions
  - pheatmap, EnhancedVolcano, patchwork
  - RColorBrewer, dendextend

## Scripts Overview

| Script | Input | Output | Description |
|--------|-------|--------|-------------|
| `plot_assembly_stats.R` | Contig statistics CSV | Length/coverage distributions | Assembly QC plots |
| `plot_pca.R` | Count matrix + metadata | PCA ordination + PERMANOVA | Beta-diversity analysis |
| `plot_volcano.R` | DA results | Volcano plot | Differential abundance visualization |
| `plot_heatmap.R` | DA results + DGE object | Heatmap | Expression patterns |
| `plot_taxonomy.R` | Annotated contigs | Circular barplot | Taxonomic composition |

---

## Usage

### 1. Assembly Statistics
```bash
# Edit script
nano plot_assembly_stats.R

# Set paths:
# CONTIG_STATS_CSV: Output from extract_contig_headers.py
# OUTPUT_DIR: Output directory for plots

# Run
Rscript plot_assembly_stats.R
```

**Outputs**:
- Length distributions (linear and log scale)
- Coverage distributions
- Length vs. coverage hexbin plots
- Circular vs. non-circular comparisons
- Summary statistics CSV

---

### 2. PCA Ordination
```bash
# Edit script
nano plot_pca.R

# Set paths:
# COUNT_MATRIX_CSV: Output from build_count_matrix.sh
# METADATA_CSV: Sample metadata
# OUTPUT_DIR: Output directory

# Run
Rscript plot_pca.R
```

**Outputs**:
- PCA plot with marginal correlations
- PERMANOVA statistics box
- Dendrogram colored by shoot mass
- Compositional data transformation artifacts

**Key features**:
- Uses Aitchison distance (appropriate for compositional data)
- Zero imputation with Bayesian-Multiplicative replacement
- CLR transformation for Euclidean space
- PERMANOVA for statistical testing

---

### 3. Volcano Plot
```bash
# Source functions in R
source("plot_volcano.R")

# Create volcano plot
p <- plot_volcano(
  df = da_results,
  title = "Differential Abundance",
  cutoff_mode = "empirical",  # or "BH"
  tau = 0.001,                 # empirical threshold
  logfc_cutoff = 0.5,
  label_n_genes = 20
)

# Save
ggsave("volcano_plot.png", p, width = 8, height = 6, dpi = 300)
```

**Options**:
- `cutoff_mode`: "BH" (Benjamini-Hochberg) or "empirical" (permutation-based)
- `fdr_cutoff`: FDR threshold for BH mode (default: 0.05)
- `tau`: P-value threshold for empirical mode
- `logfc_cutoff`: Log fold change threshold
- `label_n_genes`: Number of top genes to label

---

### 4. Heatmap
```bash
# Source functions in R
source("plot_heatmap.R")

# Create heatmap
plot_heatmap_significant(
  sig_df = da_results,
  dge_object = dge_fit,
  metadata_df = metadata,
  col_annotation_cols = "shoot_mass",
  title = "Significant Contigs",
  cutoff_mode = "empirical",
  tau = 0.001,
  logfc_cutoff = 0.5,
  filename = "heatmap.png"
)
```

**Options**:
- `col_annotation_cols`: Metadata columns to annotate (e.g., treatment, biomass)
- `clustering_distance_rows`: Distance metric for rows ("correlation", "euclidean")
- `clustering_method`: Hierarchical clustering method ("ward.D2", "average")
- `show_col_names`: Show sample names (default: FALSE)

---

### 5. Taxonomy Circular Barplot
```bash
# Edit script
nano plot_taxonomy.R

# Set paths:
# ANNOTATED_CONTIGS_CSV: Output from annotate_contigs.R
# OUTPUT_PREFIX: "circular_contigs" or "DA_contigs"

# Run
Rscript plot_taxonomy.R
```

**Outputs**:
- Circular barplot showing phylum → genus hierarchy
- Taxonomy summary table

**Customization**:
- `MIN_COUNT`: Filter taxa with <N contigs (default: 2)
- `LABEL_QUANTILE`: Label top N% by count (default: 0.5)

---

## Integration Example
```bash
# Complete visualization workflow
cd scripts/05-visualization

# 1. Assembly QC
Rscript plot_assembly_stats.R

# 2. Beta-diversity
Rscript plot_pca.R

# 3. Differential abundance (in R)
source("plot_volcano.R")
source("plot_heatmap.R")

# Load DA results
da_results <- read.csv("../04-abundance-analysis/DA_contigs.csv")
dge_fit <- readRDS("../04-abundance-analysis/edger_fit.rds")

# Create plots
p_volcano <- plot_volcano_no_title(
  df = da_results,
  cutoff_mode = "empirical",
  tau = 0.001,
  logfc_cutoff = 0.5,
  label_n_genes = 20
)
ggsave("volcano_empirical.png", p_volcano, width = 8, height = 6, dpi = 300)

plot_heatmap_significant_no_title(
  sig_df = da_results,
  dge_object = dge_fit$y,
  metadata_df = metadata,
  col_annotation_cols = "shoot_mass",
  cutoff_mode = "empirical",
  tau = 0.001,
  filename = "heatmap_empirical.png"
)

# 4. Taxonomy
Rscript plot_taxonomy.R
```

---

## Plot Customization

### Color Schemes
```R
# Assembly plots
fill_color <- "#A5C9CA"
line_color <- "#395B64"
highlight_color <- "#D1495B"

# PCA
shoot_mass_gradient <- c("lightgrey", "darkgreen")

# Volcano
colors <- c("grey30", "forestgreen", "lightblue", "pink")

# Taxonomy
# Automatic color selection based on number of phyla
# ≤12 phyla: RColorBrewer Set3
# >12 phyla: rainbow palette
```

### Font Sizes
```R
# For presentations
base_size <- 14
label_size <- 4

# For publications
base_size <- 11
label_size <- 3
```

---

## Troubleshooting

### PCA: Zero imputation too slow

**Solution**: Save `imputed_data.rds` and reuse:
```R
if (file.exists("imputed_data.rds")) {
  imputed_data <- readRDS("imputed_data.rds")
} else {
  # Compute and save
}
```

---

### Heatmap: Too many contigs

**Solution**: Filter to top N:
```R
top_contigs <- da_results %>%
  arrange(PValue) %>%
  head(50)
```

---

### Taxonomy plot: Overlapping labels

**Solution**: Adjust `LABEL_QUANTILE`:
```R
LABEL_QUANTILE <- 0.7  # Label only top 30%
```

---

## References

- EnhancedVolcano: https://github.com/kevinblighe/EnhancedVolcano
- pheatmap: https://cran.r-project.org/package=pheatmap
- vegan: https://cran.r-project.org/package=vegan
- compositions: https://cran.r-project.org/package=compositions