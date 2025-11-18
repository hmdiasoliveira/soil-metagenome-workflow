#!/usr/bin/env Rscript
# ============================================================================
# Contig Analysis and Filtering
# ============================================================================
# Description: Analyze assembly statistics and filter contigs by criteria
# Input: CSV file with contig metadata (from extract_contig_headers.py)
# Output: Summary statistics, filtered contig lists, and plots
# ============================================================================

#---- SETUP ----#
library(tidyverse)
library(plotly)
library(htmlwidgets)

# I/O paths (replace with your paths)
stats_file_path <- "<INPUT_CONTIG_HEADERS_CSV>"
out_dir <- "<OUTPUT_ANALYSIS_DIR>"
dir.create(out_dir, recursive = TRUE)

# Load data
contig_data <- read_csv(stats_file_path, show_col_types = FALSE)

# Human-readable labels for plots
fancy_labels <- function(x) {
  ifelse(x >= 1e6, paste0(x / 1e6, "M"),
         ifelse(x >= 1e3, paste0(x / 1e3, "k"),
                as.character(x)))
}

#---- CALCULATE N50 ----#
# Sort contig lengths in descending order
sorted_lengths <- sort(contig_data$length, decreasing = TRUE)

# Compute cumulative sum
cumulative_lengths <- cumsum(sorted_lengths)

# Total assembly length
total_length <- sum(sorted_lengths)

# N50 is the contig length at which cumulative length >= 50% of total
n50_index <- which(cumulative_lengths >= total_length * 0.5)[1]
n50 <- sorted_lengths[n50_index]

# Output N50 result
cat("Total assembly length:", total_length, "\n")
cat("N50 contig length:", n50, "\n")

#---- FILTER CONTIGS BY CRITERIA ----#

# Filter: ≥10kb OR ≥10x coverage
contigs_10kb_or_10x <- contig_data %>%
  filter(length >= 1e4 | coverage >= 10) %>%
  distinct(contig_id, .keep_all = TRUE)

# Filter: ≥20kb OR ≥20x coverage
contigs_20kb_or_20x <- contig_data %>%
  filter(length >= 2e4 | coverage >= 20) %>%
  distinct(contig_id, .keep_all = TRUE)

# Filter: ≥1Mb
contigs_1Mb <- contig_data %>% 
  filter(length >= 1e6) %>%
  distinct(contig_id, .keep_all = TRUE)

# Filter: Circular contigs
contigs_circular <- contig_data %>% 
  filter(circular == "yes") %>% 
  distinct(contig_id, .keep_all = TRUE)

# Filter: Circular AND (≥10kb OR ≥10x)
contigs_circular_10kb_or_10x <- contig_data %>%
  filter(length >= 1e4 | coverage >= 10,
         circular == "yes") %>%
  distinct(contig_id, .keep_all = TRUE)

# Filter: Circular AND (≥20kb OR ≥20x)
contigs_circular_20kb_or_20x <- contig_data %>%
  filter(length >= 2e4 | coverage >= 20,
         circular == "yes") %>%
  distinct(contig_id, .keep_all = TRUE)

# Filter: Circular AND ≥1Mb
contigs_circular_1Mb <- contig_data %>%
  filter(length >= 1e6,
         circular == "yes") %>%
  distinct(contig_id, .keep_all = TRUE)

#---- SUMMARY STATISTICS ----#
# Create comprehensive QC summary table
qc_counts <- t(tibble(
  n_contigs = nrow(contig_data),
  tot = sum(contig_data$length),
  n_ge_10kb = sum(contig_data$length >= 1e4, na.rm = TRUE),
  tot_10kb =  sum(contig_data$length[contig_data$length >= 1e4], na.rm = TRUE),
  n_ge_20kb = sum(contig_data$length >= 2e4, na.rm = TRUE),
  tot_20kb =  sum(contig_data$length[contig_data$length >= 2e4], na.rm = TRUE),
  n_ge_10x = sum(contig_data$coverage >= 10, na.rm = TRUE),
  tot_10x =  sum(contig_data$length[contig_data$coverage >= 10], na.rm = TRUE),
  n_ge_20x = sum(contig_data$coverage >= 20, na.rm = TRUE),
  tot_20x =  sum(contig_data$length[contig_data$coverage >= 20], na.rm = TRUE),
  n_ge_10kb_or_10x = nrow(contigs_10kb_or_10x),
  tot_10kb_or_10x = sum(contigs_10kb_or_10x$length, na.rm = TRUE),
  n_ge_20kb_or_20x = nrow(contigs_20kb_or_20x),
  tot_20kb_or_20x = sum(contigs_20kb_or_20x$length, na.rm = TRUE),
  n_ge_1Mb = nrow(contigs_1Mb),
  n_circular = nrow(contigs_circular),
  tot_circular = sum(contigs_circular$length, na.rm = TRUE),
  n_circular_ge_10kb_or_10x = sum(contigs_10kb_or_10x$circular == "yes", na.rm = TRUE),
  tot_circular_10kb_or_10x = sum(contigs_10kb_or_10x$length[contigs_10kb_or_10x$circular == "yes"], na.rm = TRUE),
  n_circular_ge_20kb_or_20x = sum(contigs_20kb_or_20x$circular == "yes", na.rm = TRUE),
  tot_circular_20kb_or_20x = sum(contigs_20kb_or_20x$length[contigs_20kb_or_20x$circular == "yes"], na.rm = TRUE),
  n_circular_ge_1Mb = sum(contigs_1Mb$circular == "yes", na.rm = TRUE),
  tot_circular_1Mb = sum(contigs_1Mb$length[contigs_1Mb$circular == "yes"], na.rm = TRUE)
))

print(qc_counts)

#---- SAVE OUTPUTS ----#

# Save summary statistics
write.csv(qc_counts, file.path(out_dir, "postassembly_qc_counts.csv"), row.names = TRUE)

# Save filtered contig lists
write_lines(contigs_1Mb$contig_id, file.path(out_dir, "contigs_list_1Mb.txt"))
write_lines(contigs_10kb_or_10x$contig_id, file.path(out_dir, "contigs_list_10kb_or_10x.txt"))
write_lines(contigs_20kb_or_20x$contig_id, file.path(out_dir, "contigs_list_20kb_or_20x.txt"))
write_lines(contigs_circular$contig_id, file.path(out_dir, "circular_contigs_list.txt"))
write_lines(contigs_circular_10kb_or_10x$contig_id, file.path(out_dir, "circular_contigs_list_10kb_or_10x.txt"))
write_lines(contigs_circular_20kb_or_20x$contig_id, file.path(out_dir, "circular_contigs_list_20kb_or_20x.txt"))

cat("\nAnalysis complete! Results saved to:", out_dir, "\n")