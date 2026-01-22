# ============================================================================
# Metagenome Count Matrix QC + Filtering
# ============================================================================
# Purpose:
#   - Clean and standardize sample names
#   - Generate QC plots for count matrix
#   - Filter low-quality samples and contigs
#   - Export filtered count matrix and metadata for downstream analysis
#
# Input:
#   - count_matrix.csv (rows = samples, cols = contigs)
#   - metadata.csv (sample metadata with sample_name column)
#
# Output:
#   - count_matrix_filtered.csv (filtered count matrix)
#   - metadata_filtered.csv (metadata for retained samples)
#   - QC plots in <OUTPUT_DIR>/plots/
# ============================================================================

library(tidyverse)
library(scales)

# Input files
count_matrix_file <- "<INPUT_COUNT_MATRIX_CSV>"  # e.g., "count_matrix.csv"
metadata_file <- "<INPUT_METADATA_CSV>"          # e.g., "metadata.csv"

# Output directory
output_dir <- "<OUTPUT_DIR>"                     # e.g., "output/qc_filtering"

# Filtering parameters
cutoff <- <MIN_LIBRARY_SIZE>                     # e.g., 10000 (minimum total counts per sample)
min_samples <- <MIN_PREVALENCE>                  # e.g., 3 (minimum samples a contig must appear in), 0.05 * nrow(count_matrix)
min_total_reads <- <MIN_TOTAL_READS>             # e.g., 100 (minimum total reads for a contig)

# Create output directories
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)

# ---- LOAD DATA ----

cat("Loading data...\n")
metadata <- read.csv(metadata_file)
count_matrix <- read.csv(count_matrix_file, row.names = 1)

cat("Initial dimensions:\n")
cat("  Samples (rows):", nrow(count_matrix), "\n")
cat("  Contigs (cols):", ncol(count_matrix), "\n\n")

#  ---- CLEAN SAMPLE NAMES ---- 

# Standardize sample names to format: {ID}_{condition}_{rep1}_{rep2}
# Adjust regex patterns based on your naming convention

cat("Cleaning sample names...\n")

# Clean metadata sample names
# Pattern: extract components and reformat
# Adjust this regex pattern to match your naming convention
m <- str_match(metadata$sample_name, "^(\\d+_[^_]+)_(\\d+)_(\\d+)$")
metadata$sample_name <- sprintf("%s_%02d_%02d", m[,2], as.integer(m[,3]), as.integer(m[,4]))

# Clean count matrix row names
# Adjust this regex pattern to match your naming convention
rn <- rownames(count_matrix)
m <- regexec("^[^_]+_(\\d+)_([^_]+)_(\\d+)_(\\d+)_.*$", rn)
mm <- regmatches(rn, m)

new_rn <- vapply(mm, function(x) {
  # x = c(fullmatch, group1, group2, group3, group4)
  sprintf("%s_%s_%02d_%02d", x[2], x[3], as.integer(x[4]), as.integer(x[5]))
}, character(1))

rownames(count_matrix) <- new_rn

# Find common samples between metadata and count matrix
common_samples <- intersect(metadata$sample_name, rownames(count_matrix))

if (length(common_samples) == 0) {
  stop("ERROR: No overlapping sample names between metadata and count_matrix after cleaning.\n",
       "Check that sample naming patterns match in both files.")
}

# Keep only shared samples
metadata <- metadata %>% filter(sample_name %in% common_samples)
count_matrix <- count_matrix[common_samples, , drop = FALSE]

# Reorder metadata to match count matrix row order
metadata <- metadata[match(rownames(count_matrix), metadata$sample_name), ]
stopifnot(all(metadata$sample_name == rownames(count_matrix)))

cat("After name cleaning:\n")
cat("  Samples:", nrow(count_matrix), "\n")
cat("  Contigs:", ncol(count_matrix), "\n\n")


#  ---- COMPUTE QC STATISTICS ---- 


cat("Computing QC statistics...\n")

# Sample statistics
sample_stats <- tibble(
  sample_name  = rownames(count_matrix),
  total_counts = rowSums(count_matrix),
  prevalence   = rowSums(count_matrix > 0)
)

# Contig statistics
contig_stats <- tibble(
  contig_id    = colnames(count_matrix),
  total_counts = colSums(count_matrix),
  prevalence   = colSums(count_matrix > 0)
)


#  ---- GENERATE QC PLOTS ---- 


cat("Generating QC plots...\n")

# 1. Library size distribution
p1 <- ggplot(sample_stats, aes(x = total_counts)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  scale_x_log10(labels = comma) +
  geom_vline(xintercept = cutoff, linetype = "dashed", color = "red", size = 1) +
  labs(x = "Total counts per sample (log10)", 
       y = "Number of samples",
       title = "Library Size Distribution",
       subtitle = paste("Red line: cutoff =", comma(cutoff))) +
  theme_bw()

ggsave(file.path(output_dir, "plots", "library_size_distribution.png"), 
       p1, width = 8, height = 6)

# 2. Contig prevalence distribution
p2 <- ggplot(contig_stats, aes(x = prevalence)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  geom_vline(xintercept = min_samples, linetype = "dashed", color = "red", size = 1) +
  labs(x = "Number of samples where contig > 0", 
       y = "Number of contigs",
       title = "Contig Prevalence Distribution",
       subtitle = paste("Red line: min_samples =", min_samples)) +
  theme_bw()

ggsave(file.path(output_dir, "plots", "contig_prevalence_distribution.png"), 
       p2, width = 8, height = 6)

# 3. Contig total abundance distribution
p3 <- ggplot(contig_stats, aes(x = total_counts)) +
  geom_histogram(bins = 80, fill = "steelblue", color = "black") +
  scale_x_log10(labels = comma) +
  geom_vline(xintercept = min_total_reads, linetype = "dashed", color = "red", size = 1) +
  labs(x = "Total counts per contig (log10)", 
       y = "Number of contigs",
       title = "Contig Total Abundance Distribution",
       subtitle = paste("Red line: min_total_reads =", comma(min_total_reads))) +
  theme_bw()

ggsave(file.path(output_dir, "plots", "contig_abundance_distribution.png"), 
       p3, width = 8, height = 6)


#  ---- FILTER SAMPLES BY LIBRARY SIZE ---- 


cat("Filtering samples...\n")

keep_samples <- sample_stats %>% 
  filter(total_counts >= cutoff) %>% 
  pull(sample_name)

count_matrix_f <- count_matrix[keep_samples, , drop = FALSE]
metadata_f <- metadata %>% filter(sample_name %in% keep_samples)

# Reorder metadata to match filtered count matrix
metadata_f <- metadata_f[match(rownames(count_matrix_f), metadata_f$sample_name), ]
stopifnot(all(metadata_f$sample_name == rownames(count_matrix_f)))

cat("After sample filtering (total_counts >= ", comma(cutoff), "):\n", sep = "")
cat("  Samples:", nrow(count_matrix_f), "\n")
cat("  Contigs:", ncol(count_matrix_f), "\n")
cat("  Removed:", nrow(count_matrix) - nrow(count_matrix_f), "samples\n\n")


#  ---- FILTER CONTIGS BY PREVALENCE AND ABUNDANCE ---- 


cat("Filtering contigs...\n")

# Keep contigs that meet EITHER criterion:
# - Present in >= min_samples samples, OR
# - Total reads >= min_total_reads across all samples
keep_contigs <- (colSums(count_matrix_f > 0) >= min_samples) | 
                (colSums(count_matrix_f) >= min_total_reads)

count_matrix_ff <- count_matrix_f[, keep_contigs, drop = FALSE]

cat("After contig filtering (prevalence >= ", min_samples, 
    " OR total_reads >= ", comma(min_total_reads), "):\n", sep = "")
cat("  Samples:", nrow(count_matrix_ff), "\n")
cat("  Contigs:", ncol(count_matrix_ff), "\n")
cat("  Removed:", ncol(count_matrix_f) - ncol(count_matrix_ff), "contigs\n\n")


#  ---- SAVE FILTERED DATA ---- 


cat("Saving filtered data...\n")

# Save filtered count matrix (rows = samples, cols = contigs)
write.csv(count_matrix_ff, 
          file.path(output_dir, "count_matrix_filtered.csv"),
          row.names = TRUE)

# Save filtered metadata
write.csv(metadata_f, 
          file.path(output_dir, "metadata_filtered.csv"),
          row.names = FALSE)

cat("\nFinal summary:\n")
cat("  Retained samples:", nrow(count_matrix_ff), "/", nrow(count_matrix), 
    sprintf("(%.1f%%)\n", 100 * nrow(count_matrix_ff) / nrow(count_matrix)))
cat("  Retained contigs:", ncol(count_matrix_ff), "/", ncol(count_matrix),
    sprintf("(%.1f%%)\n", 100 * ncol(count_matrix_ff) / ncol(count_matrix)))

cat("\nOutputs saved to:\n")
cat("  -", file.path(output_dir, "count_matrix_filtered.csv"), "\n")
cat("  -", file.path(output_dir, "metadata_filtered.csv"), "\n")
cat("  - QC plots in", file.path(output_dir, "plots/"), "\n")

cat("\nQC and filtering complete!\n")