#!/usr/bin/env Rscript
# ============================================================================
# Differential Abundance Analysis with edgeR
# ============================================================================
# Description: Identify differentially abundant contigs using edgeR
# Input: Count matrix (samples × contigs) and sample metadata
# Output: edgeR results with p-values, log fold changes, and FDR
# ============================================================================

library(edgeR)
library(tidyverse)

# ---- I/O Paths ----
count_matrix_file <- "<INPUT_COUNT_MATRIX_CSV>"
metadata_file <- "<INPUT_METADATA_CSV>"
output_dir <- "<OUTPUT_EDGER_DIR>"
dir.create(output_dir, recursive = TRUE)

# ---- Load Data ----
cat("Loading count matrix and metadata...\n")

# Load count matrix (contigs × samples)
count_matrix <- read.csv(count_matrix_file, row.names = 1)

# Load metadata
metadata <- read.csv(metadata_file)

# ---- Validate and Align Data ----
# Keep only shared samples
barcodes_shared <- intersect(colnames(count_matrix), metadata$barcode)
stopifnot(length(barcodes_shared) > 2)

# Subset and order by shared barcodes
count_matrix <- count_matrix[, barcodes_shared]
metadata <- metadata %>%
  filter(barcode %in% barcodes_shared) %>%
  arrange(barcode)

cat(paste("Analyzing", ncol(count_matrix), "samples and", 
          nrow(count_matrix), "contigs\n"))

# ---- Add Categorical Classification ----
# Example: Classify samples by biomass (adjust as needed)
metadata <- metadata %>%
  mutate(biomass_class = cut(shoot_mass,
                             breaks = c(0, 0.5, 1.0),
                             labels = c("small", "large"),
                             right = FALSE))
metadata$biomass_class <- factor(metadata$biomass_class)

# ---- edgeR Helper Function ----
run_edgeR <- function(count_matrix, metadata, model_formula) {
  # Create design matrix
  design <- model.matrix(model_formula, data = metadata)
  
  # Create DGEList object
  y <- DGEList(counts = count_matrix)
  
  # Filter low-expression contigs
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  cat(paste("Retained", sum(keep), "contigs after filtering\n"))
  
  # Normalize
  y <- calcNormFactors(y, method = "TMMwsp")
  
  # Estimate dispersion
  y <- estimateDisp(y, design)
  
  # Fit model
  fit <- glmQLFit(y, design)
  
  # Test last coefficient (treatment effect)
  k <- ncol(design)
  qlf <- glmQLFTest(fit, coef = k)
  
  # Extract results
  tab <- topTags(qlf, n = Inf)$table
  tab$FDR <- p.adjust(tab$PValue, method = "BH")
  tab$direction <- ifelse(tab$logFC > 0, "enriched_in_large", "enriched_in_small")
  tab$contig <- rownames(tab)
  
  list(table = tab, y = y, design = design, fit = fit)
}

# ---- Run Differential Abundance Analysis ----

# Categorical analysis (small vs. large)
cat("\nRunning categorical analysis (biomass class)...\n")
results_categorical <- run_edgeR(count_matrix, metadata, ~ biomass_class)

# Continuous analysis (shoot mass as continuous variable)
cat("\nRunning continuous analysis (shoot mass)...\n")
results_continuous <- run_edgeR(count_matrix, metadata, ~ shoot_mass)

# ---- Save Results ----
cat("\nSaving results...\n")

# Categorical results
write.csv(results_categorical$table, 
          file.path(output_dir, "edger_results_categorical.csv"),
          row.names = FALSE)

# Continuous results
write.csv(results_continuous$table,
          file.path(output_dir, "edger_results_continuous.csv"),
          row.names = FALSE)

# Save DGE objects for downstream analysis
saveRDS(results_categorical, file.path(output_dir, "edger_categorical_fit.rds"))
saveRDS(results_continuous, file.path(output_dir, "edger_continuous_fit.rds"))

# ---- Summary Statistics ----
cat("\n=== Categorical Analysis Summary ===\n")
cat(paste("Total contigs tested:", nrow(results_categorical$table), "\n"))
cat(paste("Significant at FDR < 0.05:", 
          sum(results_categorical$table$FDR < 0.05), "\n"))
cat(paste("  Enriched in large:", 
          sum(results_categorical$table$FDR < 0.05 & 
              results_categorical$table$direction == "enriched_in_large"), "\n"))
cat(paste("  Enriched in small:", 
          sum(results_categorical$table$FDR < 0.05 & 
              results_categorical$table$direction == "enriched_in_small"), "\n"))

cat("\n=== Continuous Analysis Summary ===\n")
cat(paste("Total contigs tested:", nrow(results_continuous$table), "\n"))
cat(paste("Significant at FDR < 0.05:", 
          sum(results_continuous$table$FDR < 0.05), "\n"))
cat(paste("  Positive correlation:", 
          sum(results_continuous$table$FDR < 0.05 & 
              results_continuous$table$logFC > 0), "\n"))
cat(paste("  Negative correlation:", 
          sum(results_continuous$table$FDR < 0.05 & 
              results_continuous$table$logFC < 0), "\n"))

cat("\nAnalysis complete!\n")