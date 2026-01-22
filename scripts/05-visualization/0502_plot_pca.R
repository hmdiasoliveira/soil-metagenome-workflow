#!/usr/bin/env Rscript
# ============================================================================
# PCA Ordination with PERMANOVA
# ============================================================================
# Description: Perform PCA ordination and PERMANOVA on compositional data
# Input: Count matrix (samples × contigs) and metadata
# Output: PCA plots with marginal correlations and PERMANOVA statistics
# ============================================================================

library(vegan)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(dplyr)
library(compositions)
library(zCompositions)
library(dendextend)
library(colorspace)

set.seed(123)
theme_set(theme_bw(base_size = 12))

# ---- I/O Paths ----
COUNT_MATRIX_CSV <- "<INPUT_COUNT_MATRIX_CSV>"
METADATA_CSV <- "<INPUT_METADATA_CSV>"
OUTPUT_DIR <- "<OUTPUT_PCA_DIR>"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Load Data ----
cat("Loading count matrix and metadata...\n")

matrix <- read.csv(COUNT_MATRIX_CSV, row.names = 1)
metadata <- read.csv(METADATA_CSV)

# ---- Filter Low-Abundance Contigs ----
min_samples <- 3
min_total_reads <- 10

matrix_filtered <- matrix[rowSums(matrix > 0) >= min_samples | 
                           rowSums(matrix) >= min_total_reads, ]

cat("Original contigs:", nrow(matrix), "\n")
cat("Filtered contigs:", nrow(matrix_filtered), "\n")
cat("Samples:", ncol(matrix_filtered), "\n\n")

# Transpose for vegan (samples as rows)
matrix_t <- t(matrix_filtered)

# Ensure metadata order matches matrix
metadata <- metadata[match(colnames(matrix_filtered), metadata$barcode), ]
rownames(metadata) <- metadata$barcode

# Convert variables to appropriate types
metadata$host_plant <- as.factor(metadata$host_plant)
metadata$shoot_mass_scaled <- as.numeric(scale(metadata$shoot_mass))

# ---- Compositional Data Transformation ----
cat("Performing compositional data transformation...\n")

# Zero imputation using Bayesian-Multiplicative replacement
imputed_data_file <- file.path(OUTPUT_DIR, "imputed_data.rds")

# Note: GBM (Geometric Bayesian Multiplicative) is very computationally expensive and memory-hungry. 
#       It may be killed by OOM after running for a long time.
if (file.exists(imputed_data_file)) {
  cat("Loading pre-computed imputed data...\n")
  imputed_data <- readRDS(imputed_data_file)
} else {
  cat("Computing zero imputation (this may take a while)...\n")
  imputed_data <- cmultRepl(
    X = matrix_t,
    label = 0,
    method = "GBM",  # Geometric Bayesian Multiplicative
    output = "p-counts"
  )
  saveRDS(imputed_data, imputed_data_file)
}

# Centered Log Ratio (CLR) transformation
# The CLR transformation moves data from the compositional simplex to Euclidean space,
# enabling the use of standard multivariate statistics.
clr_data <- cenLR(imputed_data)$x

# If cmultRepl is crashing, try to add pseudo-count (usually 1 or 0.5) to handle zeros
#matrix_pseudo <- matrix + 0.5
#res_cenLR <- cenLR(matrix_pseudo)
#clr_data <- res_cenLR$x.clr
#or, clr_data <- clr(acomp(matrix_pseudo))

# Calculate Aitchison distance
# The Aitchison distance is the Euclidean distance calculated on CLR-transformed data.
# It is the statistically appropriate beta-diversity metric for compositional data.
# aDist() {composisions} expects raw positive compositions. CLR coordinates (clr_data), which contain negative values by design.
#dist_aitchison <- aDist(clr_data)
dist_aitchison <- dist(clr_data, method = "euclidean")
# or use vegan::vegdist(clr_data, method = "euclidean")

# ---- PERMANOVA Analysis ----
cat("\nRunning PERMANOVA...\n")

# Test multiple models
permanova_shoot_mass <- adonis2(dist_aitchison ~ shoot_mass, 
                                data = metadata, 
                                permutations = 999)

permanova_host_plant <- adonis2(dist_aitchison ~ host_plant, 
                                data = metadata, 
                                permutations = 999)

permanova_full <- adonis2(dist_aitchison ~ host_plant + shoot_mass, 
                          data = metadata, 
                          permutations = 999,
                          by = "margin")

cat("\n=== PERMANOVA Results ===\n")
cat("\nShoot Mass:\n")
print(permanova_shoot_mass)
cat("\nHost Plant:\n")
print(permanova_host_plant)
cat("\nFull Model (marginal effects):\n")
print(permanova_full)

# Extract statistics for plotting
r2_shoot <- permanova_shoot_mass$R2[1]
p_shoot <- permanova_shoot_mass$`Pr(>F)`[1]

# ---- PCA Ordination ----
cat("\nPerforming PCA...\n")

pca_result <- prcomp(clr_data, scale. = FALSE)
var_exp <- round((pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100, 1)

# Create PCA data frame
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  host_plant = metadata$host_plant,
  shoot_mass = metadata$shoot_mass,
  sample = metadata$barcode
)

# ---- Define Plot Limits ----
pc1_range <- range(pca_df$PC1, na.rm = TRUE)
pc2_range <- range(pca_df$PC2, na.rm = TRUE)

pc1_padding <- diff(pc1_range) * 0.05
pc2_padding <- diff(pc2_range) * 0.05

pc1_limits <- c(pc1_range[1] - pc1_padding, pc1_range[2] + pc1_padding)
pc2_limits <- c(pc2_range[1] - pc2_padding, pc2_range[2] + pc2_padding)

# ---- Main PCA Plot ----
main_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, 
                                 color = shoot_mass, 
                                 shape = host_plant)) +
  geom_point(size = 4, alpha = 0.8, stroke = 1) +
  scale_color_gradient(
    low = "lightgrey",
    high = "darkgreen",
    name = "Shoot Mass (g)"
  ) +
  scale_shape_manual(
    values = c(15, 16, 17, 18, 19, 25, 21),
    name = "Host Plant"
  ) +
  scale_x_continuous(limits = pc1_limits, expand = c(0, 0)) +
  scale_y_continuous(limits = pc2_limits, expand = c(0, 0)) +
  labs(
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)")
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.margin = margin(0, 0, 4, 4)
  )

# ---- Marginal Scatter Plots ----

# Calculate correlations
lm_pc1 <- lm(shoot_mass ~ PC1, data = pca_df)
lm_pc2 <- lm(PC2 ~ shoot_mass, data = pca_df)

r2_pc1 <- summary(lm_pc1)$r.squared
r2_pc2 <- summary(lm_pc2)$r.squared
p_lm_pc1 <- summary(lm_pc1)$coefficients[2, 4]
p_lm_pc2 <- summary(lm_pc2)$coefficients[2, 4]

# Format p-values
format_pval <- function(p) {
  if (p < 0.001) return("< 0.001")
  else if (p < 0.01) return("< 0.01")
  else if (p < 0.05) return(paste0("= ", sprintf("%.3f", p)))
  else return(paste0("= ", sprintf("%.2f", p)))
}

pc1_annotation <- paste0("R² = ", sprintf("%.3f", r2_pc1), ", p ", format_pval(p_lm_pc1))
pc2_annotation <- paste0("R² = ", sprintf("%.3f", r2_pc2), ", p ", format_pval(p_lm_pc2))

# PC1 vs shoot_mass scatter
pc1_scatter <- ggplot(pca_df, aes(x = PC1, y = shoot_mass)) +
  geom_point(aes(color = shoot_mass, shape = host_plant), size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", linewidth = 0.8) +
  annotate("text", 
           x = pc1_limits[1] + diff(pc1_limits) * 0.05,
           y = max(pca_df$shoot_mass, na.rm = TRUE),
           label = pc1_annotation,
           hjust = 0.25, vjust = -1, size = 3.5, fontface = "bold") +
  scale_color_gradient(low = "lightgrey", high = "darkgreen", guide = "none") +
  scale_shape_manual(values = c(15, 16, 17, 18, 19, 25, 21), guide = "none") +
  scale_x_continuous(limits = pc1_limits, expand = c(0, 0)) +
  labs(y = "Shoot Mass (g)", x = NULL) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(4, 0, 0, 4)
  )

# PC2 vs shoot_mass scatter
pc2_scatter <- ggplot(pca_df, aes(x = shoot_mass, y = PC2)) +
  geom_point(aes(color = shoot_mass, shape = host_plant), size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", linewidth = 0.8) +
  annotate("text", 
           x = min(pca_df$shoot_mass, na.rm = TRUE),
           y = pc2_limits[2],
           label = pc2_annotation,
           hjust = 0, vjust = 1.5, size = 3.5, fontface = "bold") +
  scale_color_gradient(low = "lightgrey", high = "darkgreen", guide = "none") +
  scale_shape_manual(values = c(15, 16, 17, 18, 19, 25, 21), guide = "none") +
  scale_y_continuous(limits = pc2_limits, expand = c(0, 0)) +
  labs(x = "Shoot Mass (g)", y = NULL) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0, 4, 4, 0)
  )

# ---- PERMANOVA Statistics Box ----
permanova_text <- paste0(
  "PERMANOVA\n",
  "R² = ", sprintf("%.2f", r2_shoot), "\n",
  "p ", ifelse(p_shoot < 0.001, "< 0.001", 
               ifelse(p_shoot < 0.01, "< 0.01", 
                      paste0("= ", sprintf("%.3f", p_shoot))))
)

permanova_box <- ggplot() +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
           fill = "white", color = "black", linewidth = 1.5) +
  annotate("text", x = 0.5, y = 0.5, 
           label = permanova_text,
           size = 4.5, fontface = "bold",
           hjust = 0.5, vjust = 0.5) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(4, 4, 4, 4))

# ---- Combine Plots ----
layout <- "
AAB
CCD
CCD
"

combined_plot <- pc1_scatter + permanova_box + 
  main_plot + pc2_scatter +
  plot_layout(
    design = layout,
    widths = c(4, 4, 1.5),
    heights = c(1, 1.5, 1.5, 1.5)
  )

# Save
ggsave(file.path(OUTPUT_DIR, "PCA_permanova_shootmass.pdf"), combined_plot, 
       width = 10, height = 8, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "PCA_permanova_shootmass.png"), combined_plot, 
       width = 10, height = 8, dpi = 300)

cat("\nPCA plot saved successfully!\n")

# ---- Hierarchical Clustering ----
cat("\nPerforming hierarchical clustering...\n")

dist_matrix <- as.dist(dist_aitchison)

# Test different linkage methods
hc_ward <- hclust(dist_matrix, method = "ward.D2")
hc_average <- hclust(dist_matrix, method = "average")
hc_complete <- hclust(dist_matrix, method = "complete")

# Calculate cophenetic correlation
coph_ward <- cor(dist_matrix, cophenetic(hc_ward))
coph_average <- cor(dist_matrix, cophenetic(hc_average))
coph_complete <- cor(dist_matrix, cophenetic(hc_complete))

cat("Cophenetic correlations:\n")
cat("  Ward.D2:  ", round(coph_ward, 3), "\n")
cat("  Average:  ", round(coph_average, 3), "\n")
cat("  Complete: ", round(coph_complete, 3), "\n")

# Use best method
hc_best <- hc_average

# Create dendrogram colored by shoot mass
shoot_mass_colors <- sequential_hcl(100, palette = "Greens 3")
shoot_mass_breaks <- cut(metadata$shoot_mass, breaks = 100)
shoot_mass_cols <- shoot_mass_colors[as.numeric(shoot_mass_breaks)]

dend <- as.dendrogram(hc_best)
labels_colors(dend) <- shoot_mass_cols[order.dendrogram(dend)]

pdf(file.path(OUTPUT_DIR, "dendrogram_by_shootmass.pdf"), width = 12, height = 6)
par(mar = c(8, 4, 4, 6))
plot(dend, main = "Hierarchical Clustering - Colored by Shoot Mass",
     ylab = "Aitchison Distance")

# Add color legend
legend_image <- as.raster(matrix(rev(shoot_mass_colors), ncol = 1))
rasterImage(legend_image, 
            xlim = par("usr")[2] * 1.02, 
            ylim = par("usr")[3], 
            xlim + (par("usr")[2] - par("usr")[1]) * 0.03,
            ylim + (par("usr")[4] - par("usr")[3]))
text(x = par("usr")[2] * 1.08, 
     y = c(par("usr")[3], par("usr")[4]), 
     labels = c(round(min(metadata$shoot_mass), 2), 
                round(max(metadata$shoot_mass), 2)),
     cex = 0.8)
dev.off()

cat("\nDendrogram saved successfully!\n")
cat(paste("Output directory:", OUTPUT_DIR, "\n"))