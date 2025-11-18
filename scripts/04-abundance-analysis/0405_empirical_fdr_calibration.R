#!/usr/bin/env Rscript
# ============================================================================
# Empirical FDR Calibration via Permutation Testing
# ============================================================================
# Description: Perform permutation-based FDR control for differential abundance
# Input: Count matrix, metadata, and edgeR results
# Output: Empirical FDR curves, calibrated p-value thresholds, plots
# ============================================================================

library(edgeR)
library(tidyverse)
library(ggplot2)

# Source plotting functions
source("<PLOT_VOLCANO_R>")  # plot_volcano.R from 05-visualization
source("<PLOT_HEATMAP_R>")  # plot_heatmap.R from 05-visualization

# ---- I/O Paths ----
count_matrix_file <- "<INPUT_COUNT_MATRIX_CSV>"
metadata_file <- "<INPUT_METADATA_CSV>"
output_dir <- "<OUTPUT_FDR_DIR>"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load Data ----
cat("Loading count matrix and metadata...\n")

count_matrix <- read.csv(count_matrix_file, row.names = 1)
metadata <- read.csv(metadata_file)

# Align data
barcodes_shared <- intersect(colnames(count_matrix), metadata$barcode)
stopifnot(length(barcodes_shared) > 2)

count_matrix <- count_matrix[, barcodes_shared]
metadata <- metadata %>%
  filter(barcode %in% barcodes_shared) %>%
  arrange(barcode)

# Add biomass classification
metadata <- metadata %>%
  mutate(biomass_class = cut(shoot_mass,
                             breaks = c(0, 0.5, 1.0),
                             labels = c("small", "large"),
                             right = FALSE))
metadata$biomass_class <- factor(metadata$biomass_class)

# ---- edgeR Helper Function ----
run_edgeR <- function(count_matrix, metadata, model_formula) {
  design <- model.matrix(model_formula, data = metadata)
  y <- DGEList(counts = count_matrix)
  
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  y <- calcNormFactors(y, method = "TMMwsp")
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  
  k <- ncol(design)
  qlf <- glmQLFTest(fit, coef = k)
  
  tab <- topTags(qlf, n = Inf)$table
  tab$FDR <- p.adjust(tab$PValue, method = "BH")
  tab$direction <- ifelse(tab$logFC > 0, "enriched_in_large", "enriched_in_small")
  tab$MAG <- rownames(tab)
  
  list(table = tab, y = y, design = design)
}

# ---- Permutation Testing ----
set.seed(202509)
B <- 500  # Number of permutations

cat(paste("Running", B, "permutations for empirical FDR calibration...\n"))

# Observed analysis
cat("Running observed analysis...\n")
fit_obs_cont <- run_edgeR(count_matrix, metadata, ~ shoot_mass)
obs_tab_cont <- fit_obs_cont$table
obs_p_cont <- obs_tab_cont$PValue

# Permutation analyses
perm_p_list <- vector("list", B)

for (b in seq_len(B)) {
  if (b %% 50 == 0) {
    cat(paste("  Permutation", b, "/", B, "\n"))
  }
  
  # Permute covariate
  meta_perm <- metadata
  meta_perm$shoot_mass <- sample(meta_perm$shoot_mass)
  
  # Run edgeR on permuted data
  fit_perm <- run_edgeR(count_matrix, meta_perm, ~ shoot_mass)
  perm_p_list[[b]] <- fit_perm$table$PValue
}

# Save permutation results
saveRDS(perm_p_list, file.path(output_dir, "permutation_p_values.rds"))
saveRDS(obs_tab_cont, file.path(output_dir, "observed_results.rds"))

# ---- Calculate Empirical FDR ----
cat("\nCalculating empirical FDR...\n")

calculate_empirical_fdr <- function(obs_p, perm_p_list) {
  # Grid of p-value thresholds
  grid <- sort(unique(c(seq(1e-6, 0.2, by = 1e-3), 0.01, 0.05, 0.1)))
  
  # Count discoveries at each threshold
  obs_hits <- sapply(grid, function(tau) sum(obs_p <= tau, na.rm = TRUE))
  
  # Count null discoveries
  null_hits_mat <- sapply(perm_p_list, function(pv) {
    sapply(grid, function(tau) sum(pv <= tau, na.rm = TRUE))
  })
  null_mean_hits <- rowMeans(null_hits_mat, na.rm = TRUE)
  
  # Empirical FDR
  empirical_fdr <- null_mean_hits / pmax(obs_hits, 1)
  
  tibble(
    tau = grid,
    obs_hits = obs_hits,
    null_mean_hits = null_mean_hits,
    empirical_fdr = pmin(empirical_fdr, 1)
  )
}

fdr_curve <- calculate_empirical_fdr(obs_p_cont, perm_p_list)

# ---- Find Tau at FDR = 0.05 ----
target_fdr <- 0.05

df_uni <- fdr_curve %>%
  arrange(empirical_fdr, tau) %>%
  distinct(empirical_fdr, .keep_all = TRUE)

approx_tau <- approx(
  x = df_uni$empirical_fdr,
  y = df_uni$tau,
  xout = target_fdr,
  ties = "ordered"
)$y

cat(paste("\nEmpirical p-value threshold (τ) at FDR =", target_fdr, ":", 
          signif(approx_tau, 4), "\n"))

# Save tau value
write.csv(data.frame(target_fdr = target_fdr, tau = approx_tau),
          file.path(output_dir, "empirical_tau.csv"),
          row.names = FALSE)

# ---- Generate Plots ----
cat("\nGenerating plots...\n")

# Plot 1: Empirical FDR Curve
p_curve <- ggplot(fdr_curve, aes(x = tau, y = empirical_fdr)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  geom_vline(xintercept = approx_tau, linetype = "dotted", color = "darkgreen") +
  annotate("text", x = approx_tau, y = 0.9,
           label = paste0("τ ≈ ", signif(approx_tau, 3)),
           angle = 90, vjust = -0.5, color = "darkgreen") +
  annotate("text", x = max(fdr_curve$tau), y = 0.05,
           label = "FDR = 0.05", hjust = 1.1, vjust = -0.5,
           color = "red", size = 3.5) +
  labs(title = "Empirical FDR Curve",
       x = expression(paste("p-value threshold ", tau)),
       y = "Empirical FDR") +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1))

ggsave(file.path(output_dir, "empirical_fdr_curve.png"), p_curve, 
       width = 7, height = 5, dpi = 300)

# Plot 2: Discovery Counts (Observed vs. Null)
p_counts <- ggplot(fdr_curve, aes(x = tau)) +
  geom_line(aes(y = obs_hits, color = "Observed"), linewidth = 1) +
  geom_line(aes(y = null_mean_hits, color = "Null mean"), linewidth = 1) +
  scale_color_manual(values = c("Observed" = "steelblue", "Null mean" = "grey40")) +
  labs(title = "Discoveries: Observed vs. Permutation Null",
       x = expression(paste("p-value threshold ", tau)),
       y = "# of discoveries", color = NULL) +
  theme_bw()

ggsave(file.path(output_dir, "observed_vs_null_discoveries.png"), p_counts, 
       width = 7, height = 5, dpi = 300)

# Plot 3: Volcano Plot with Empirical Tau Cutoff
volcano_plot_emp <- plot_volcano_no_title(
  obs_tab_cont, 
  title = NULL, 
  cutoff_mode = "empirical", 
  tau = approx_tau,
  logfc_cutoff = 0.5, 
  label_n_genes = 20, 
  show_titles = FALSE
)

ggsave(file.path(output_dir, "volcano_empirical.png"), volcano_plot_emp, 
       width = 8, height = 6, dpi = 1000)

# Plot 4: Volcano Plot with BH FDR Cutoff
volcano_plot_bh <- plot_volcano_no_title(
  obs_tab_cont, 
  title = NULL, 
  cutoff_mode = "BH", 
  fdr_cutoff = 0.05,
  logfc_cutoff = 0.5, 
  label_n_genes = 20, 
  show_titles = FALSE
)

ggsave(file.path(output_dir, "volcano_bh.png"), volcano_plot_bh, 
       width = 8, height = 6, dpi = 1000)

# Plot 5 & 6: Heatmaps (Empirical and BH)
# Prepare metadata for heatmap
metadata_aligned <- metadata
rownames(metadata_aligned) <- metadata_aligned$barcode

# Keep only samples present in DGE object
sample_ids <- colnames(edgeR::cpm(fit_obs_cont$y))
metadata_aligned <- metadata_aligned[sample_ids, , drop = FALSE]

# Sanity check
stopifnot(identical(rownames(metadata_aligned), sample_ids))

# Get top 50 contigs for heatmap
obs_tab_ordered <- obs_tab_cont[order(obs_tab_cont$FDR), ]
obs_tab_top <- obs_tab_ordered[1:min(50, nrow(obs_tab_ordered)), ]

# Heatmap with empirical tau cutoff
plot_heatmap_significant_no_title(
  sig_df = obs_tab_top, 
  dge_object = fit_obs_cont$y, 
  metadata_df = metadata_aligned,
  col_annotation_cols = "shoot_mass",
  title = NULL, 
  cutoff_mode = "empirical", 
  tau = approx_tau, 
  logfc_cutoff = 0.5,
  filename = file.path(output_dir, "heatmap_empirical.png"),
  show_title = FALSE
)

# Heatmap with BH FDR cutoff
plot_heatmap_significant_no_title(
  sig_df = obs_tab_top, 
  dge_object = fit_obs_cont$y, 
  metadata_df = metadata_aligned,
  col_annotation_cols = "shoot_mass",
  title = NULL, 
  cutoff_mode = "BH", 
  fdr_cutoff = 0.05, 
  logfc_cutoff = 0.5,
  filename = file.path(output_dir, "heatmap_bh.png"),
  show_title = FALSE
)

# ---- Save Significant Contigs ----
da_contigs <- obs_tab_cont %>%
  filter(PValue <= approx_tau, abs(logFC) >= 0.5)

sig_contigs <- obs_tab_cont %>%
  filter(PValue <= approx_tau)

effect_contigs <- sig_contigs %>%
  filter(abs(logFC) >= 0.5)

write.csv(da_contigs, 
          file.path(output_dir, "DA_contigs_empirical_fdr.csv"),
          row.names = FALSE)

write_lines(da_contigs$MAG, 
            file.path(output_dir, "DA_contigs_list.txt"))

write.csv(obs_tab_cont,
          file.path(output_dir, "all_contigs_results.csv"),
          row.names = FALSE)

# ---- Summary Statistics ----
cat("\n=== Summary ===\n")
cat(paste("Total contigs tested:", nrow(obs_tab_cont), "\n"))
cat(paste("Significant contigs (p ≤ τ):", nrow(sig_contigs), "\n"))
cat(paste("DA contigs (p ≤ τ & |logFC| ≥ 0.5):", nrow(da_contigs), "\n"))
cat(paste("  Enriched in large:", 
          sum(da_contigs$direction == "enriched_in_large"), "\n"))
cat(paste("  Enriched in small:", 
          sum(da_contigs$direction == "enriched_in_small"), "\n"))

cat("\nEmpirical FDR calibration complete!\n")
cat(paste("Output directory:", output_dir, "\n"))