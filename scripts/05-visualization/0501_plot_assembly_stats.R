#!/usr/bin/env Rscript
# ============================================================================
# Assembly Statistics Visualization
# ============================================================================
# Description: Generate comprehensive plots of assembly quality metrics
# Input: Contig statistics CSV (from extract_contig_headers.py)
# Output: Multiple plots showing length, coverage, and circularity distributions
# ============================================================================

library(tidyverse)
library(plotly)
library(htmlwidgets)

# ---- I/O Paths ----
CONTIG_STATS_CSV <- "<INPUT_CONTIG_STATS_CSV>"
OUTPUT_DIR <- "<OUTPUT_PLOTS_DIR>"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Load Data ----
cat("Loading contig statistics...\n")
contig_data <- read_csv(CONTIG_STATS_CSV, show_col_types = FALSE)
cat(paste("Loaded", nrow(contig_data), "contigs\n"))

# ---- Helper Functions ----

# Human-readable labels for axes
fancy_labels <- function(x) {
  ifelse(x >= 1e6, paste0(x / 1e6, "M"),
         ifelse(x >= 1e3, paste0(x / 1e3, "k"),
                as.character(x)))
}

# ---- Calculate Summary Statistics ----

# Sort contig lengths in descending order
sorted_lengths <- sort(contig_data$length, decreasing = TRUE)
cumulative_lengths <- cumsum(sorted_lengths)
total_length <- sum(sorted_lengths)

# Calculate N50
n50_index <- which(cumulative_lengths >= total_length * 0.5)[1]
n50 <- sorted_lengths[n50_index]

cat("\n=== Assembly Statistics ===\n")
cat("Total assembly length:", fancy_labels(total_length), "bp\n")
cat("N50 contig length:", fancy_labels(n50), "bp\n")
cat("Number of contigs:", nrow(contig_data), "\n")

# Filter contigs by quality criteria
contigs_10kb_or_10x <- contig_data %>%
  filter(length >= 1e4 | coverage >= 10) %>%
  distinct(contig_id, .keep_all = TRUE)

contigs_20kb_or_20x <- contig_data %>%
  filter(length >= 2e4 | coverage >= 20) %>%
  distinct(contig_id, .keep_all = TRUE)

contigs_1Mb <- contig_data %>% 
  filter(length >= 1e6) %>%
  distinct(contig_id, .keep_all = TRUE)

contigs_circular <- contig_data %>% 
  filter(circular == "yes") %>% 
  distinct(contig_id, .keep_all = TRUE)

contigs_circular_10kb_or_10x <- contig_data %>%
  filter(length >= 1e4 | coverage >= 10, circular == "yes") %>%
  distinct(contig_id, .keep_all = TRUE)

contigs_circular_1Mb <- contig_data %>%
  filter(length >= 1e6, circular == "yes") %>%
  distinct(contig_id, .keep_all = TRUE)

# Create comprehensive QC summary
qc_counts <- t(tibble(
  n_contigs = nrow(contig_data),
  tot = sum(contig_data$length),
  n_ge_10kb = sum(contig_data$length >= 1e4, na.rm = TRUE),
  tot_10kb = sum(contig_data$length[contig_data$length >= 1e4], na.rm = TRUE),
  n_ge_20kb = sum(contig_data$length >= 2e4, na.rm = TRUE),
  tot_20kb = sum(contig_data$length[contig_data$length >= 2e4], na.rm = TRUE),
  n_ge_10x = sum(contig_data$coverage >= 10, na.rm = TRUE),
  tot_10x = sum(contig_data$length[contig_data$coverage >= 10], na.rm = TRUE),
  n_ge_20x = sum(contig_data$coverage >= 20, na.rm = TRUE),
  tot_20x = sum(contig_data$length[contig_data$coverage >= 20], na.rm = TRUE),
  n_ge_10kb_or_10x = nrow(contigs_10kb_or_10x),
  tot_10kb_or_10x = sum(contigs_10kb_or_10x$length, na.rm = TRUE),
  n_ge_20kb_or_20x = nrow(contigs_20kb_or_20x),
  tot_20kb_or_20x = sum(contigs_20kb_or_20x$length, na.rm = TRUE),
  n_ge_1Mb = nrow(contigs_1Mb),
  n_circular = nrow(contigs_circular),
  tot_circular = sum(contigs_circular$length, na.rm = TRUE),
  n_circular_ge_10kb_or_10x = sum(contigs_10kb_or_10x$circular == "yes", na.rm = TRUE),
  tot_circular_10kb_or_10x = sum(contigs_10kb_or_10x$length[contigs_10kb_or_10x$circular == "yes"], na.rm = TRUE),
  n_circular_ge_1Mb = sum(contigs_1Mb$circular == "yes", na.rm = TRUE),
  tot_circular_1Mb = sum(contigs_1Mb$length[contigs_1Mb$circular == "yes"], na.rm = TRUE)
))

print(qc_counts)
write.csv(qc_counts, file.path(OUTPUT_DIR, "assembly_qc_counts.csv"), row.names = TRUE)

# Save filtered contig lists
write_lines(contigs_1Mb$contig_id, file.path(OUTPUT_DIR, "contigs_list_1Mb.txt"))
write_lines(contigs_10kb_or_10x$contig_id, file.path(OUTPUT_DIR, "contigs_list_10kb_or_10x.txt"))
write_lines(contigs_20kb_or_20x$contig_id, file.path(OUTPUT_DIR, "contigs_list_20kb_or_20x.txt"))
write_lines(contigs_circular$contig_id, file.path(OUTPUT_DIR, "circular_contigs_list.txt"))
write_lines(contigs_circular_10kb_or_10x$contig_id, file.path(OUTPUT_DIR, "circular_contigs_list_10kb_or_10x.txt"))

# ---- PLOT 1: Contig Length Distribution ----
cat("\nGenerating length distribution plots...\n")

# Convert circular to factor for better plotting
contig_data$circular <- factor(contig_data$circular,
                               levels = c("no", "yes"),
                               labels = c("non-circular", "circular"))

# 1A: Histogram (linear scale)
p1a <- ggplot(contig_data, aes(x = length)) +
  geom_histogram(bins = 100, fill = "#A5C9CA", color = "#395B64") +
  scale_x_continuous(labels = fancy_labels) +
  scale_y_continuous(labels = fancy_labels) +
  labs(x = "Contig Length (bp)", y = "Number of Contigs") +
  theme_minimal()

ggsave(file.path(OUTPUT_DIR, "contig_length_distribution.png"), p1a, width = 7, height = 5, dpi = 300)

# 1B: Histogram (log scale) with N50 annotation by circularity
n50_by_group <- contig_data %>%
  group_by(circular) %>%
  summarise(
    total_len = sum(length, na.rm = TRUE),
    N50 = {
      L <- sort(length, decreasing = TRUE)
      C <- cumsum(L)
      L[which(C >= 0.5 * sum(L))[1]]
    },
    n_ge_1Mb = sum(length >= 1e6, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(n50_by_group, file.path(OUTPUT_DIR, "n50_by_circularity.csv"), row.names = FALSE)

p1b <- ggplot(contig_data, aes(x = length)) +
  geom_histogram(bins = 100, fill = "#A5C9CA", color = "#395B64") +
  facet_wrap(~ circular, ncol = 2, scales = "free_y") +
  scale_x_log10(labels = fancy_labels) +
  geom_vline(data = n50_by_group, aes(xintercept = N50), 
             color = "#D1495B", linetype = "dashed", linewidth = 0.8) +
  geom_text(data = n50_by_group, aes(x = N50, y = Inf,
                                     label = paste0("N50 = ", fancy_labels(N50))),
            vjust = 1.3, hjust = 1.1, angle = 90, color = "#D1495B", size = 3.3) +
  geom_label(data = n50_by_group, aes(x = 1000, y = Inf,
                                      label = paste0(n_ge_1Mb, " contigs ≥ 1 Mb")),
             inherit.aes = FALSE, vjust = 4, hjust = -1,
             fill = "lightyellow", color = "black", size = 3.3) +
  coord_cartesian(clip = "off") +
  labs(x = "Contig Length (bp, log10)", y = "Count") +
  theme_minimal()

ggsave(file.path(OUTPUT_DIR, "contig_length_distribution_log.png"), p1b, width = 7, height = 5, dpi = 300)

# ---- PLOT 2: Coverage vs. Length ----
cat("Generating coverage vs. length plots...\n")

p2 <- ggplot(contig_data, aes(x = length, y = coverage)) +
  geom_hex(bins = 50) +
  scale_fill_viridis_c(trans = "log10") +
  scale_x_log10(labels = fancy_labels) +
  scale_y_log10(labels = fancy_labels) +
  geom_vline(xintercept = 10000, linetype = "dashed", color = "#D1495B", linewidth = 0.8) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "#D1495B", linewidth = 0.8) +
  labs(x = "Contig Length (bp, log10)", y = "Coverage (X, log10)", fill = "Count") +
  theme_minimal()

ggsave(file.path(OUTPUT_DIR, "contig_coverage_length_hexbin.png"), p2, width = 7, height = 5, dpi = 300)

# Coverage vs. Length by circularity
p2b <- ggplot(contig_data, aes(x = length, y = coverage)) +
  geom_hex(bins = 50) +
  scale_fill_viridis_c(trans = "log10") +
  scale_x_log10(labels = fancy_labels) +
  scale_y_log10(labels = fancy_labels) +
  geom_vline(xintercept = 10000, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  facet_wrap(~circular) +
  labs(x = "Contig Length (bp, log10)", y = "Coverage (X, log10)", fill = "Count") +
  theme_minimal()

ggsave(file.path(OUTPUT_DIR, "contigs_length_coverage_by_circularity_hexbin.png"), p2b, width = 7, height = 5, dpi = 300)

# ---- PLOT 3: Coverage Distribution ----
cat("Generating coverage distribution plots...\n")

# 3A: Histogram of coverage
p3a <- ggplot(contig_data, aes(x = coverage)) +
  geom_histogram(bins = 20, fill = "#A5C9CA", color = "white") +
  geom_vline(xintercept = c(10, 20), linetype = "dashed", color = c("red", "blue")) +
  scale_x_log10(labels = fancy_labels) +
  scale_y_continuous(labels = fancy_labels) +
  labs(x = "Coverage (X)", y = "Number of Contigs") +
  theme_minimal()

ggsave(file.path(OUTPUT_DIR, "contig_coverage_distribution.png"), p3a, width = 7, height = 5, dpi = 300)

# 3B: Violin plot by circularity
n_circular <- nrow(contig_data %>% filter(circular == "circular"))

p3b <- ggplot(contig_data, aes(x = circular, y = coverage)) +
  geom_violin(fill = "#A5C9CA", color = "#395B64", alpha = 0.6, width = 0.8, trim = TRUE) +
  geom_boxplot(width = 0.05, alpha = 0.2) +
  scale_y_log10() +
  annotate("label", x = 1.5, y = 200, 
           label = paste("circular contigs:", n_circular), 
           fill = "lightyellow") +
  labs(x = "Circular", y = "Coverage (X, log10)") +
  theme_minimal()

ggsave(file.path(OUTPUT_DIR, "contig_coverage_by_circularity.png"), p3b, width = 7, height = 5, dpi = 300)

# ---- PLOT 4: Circular Status ----
cat("Generating circular status plot...\n")

p4 <- ggplot(contig_data, aes(x = circular)) +
  geom_bar(fill = "#A5C9CA", color = "#395B64") +
  scale_y_log10(labels = fancy_labels) +
  labs(x = "Circular Status", y = "Count (log10)") +
  theme_minimal()

ggsave(file.path(OUTPUT_DIR, "circular_contigs_barplot.png"), p4, width = 7, height = 5, dpi = 300)

# ---- PLOT 5: High-Quality Contigs (≥10kb OR ≥10x) ----
cat("Generating high-quality contig plots...\n")

p5 <- ggplot(contigs_10kb_or_10x, aes(x = length, y = coverage, color = circular)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("yes" = "darkgreen", "no" = "gray80")) +
  scale_x_log10(labels = fancy_labels) +
  scale_y_log10() +
  labs(x = "Contig Length (bp, log10)",
       y = "Coverage (X, log10)",
       color = "Circular?") +
  theme_minimal()

ggsave(file.path(OUTPUT_DIR, "contigs_10kb_or_10x_by_circularity_dot.png"), p5, width = 7, height = 5, dpi = 300)

cat("\nAll assembly statistics plots generated successfully!\n")
cat(paste("Output directory:", OUTPUT_DIR, "\n"))