#!/usr/bin/env Rscript
# ============================================================================
# Taxonomy Circular Barplot
# ============================================================================
# Description: Create circular barplot of taxonomic composition
# Input: Annotated contigs with taxonomy
# Output: Circular barplot showing hierarchical taxonomy
# ============================================================================

library(tidyverse)
library(ggplot2)
library(RColorBrewer)

# ---- I/O Paths ----
ANNOTATED_CONTIGS_CSV <- "<INPUT_ANNOTATED_CONTIGS_CSV>"
OUTPUT_DIR <- "<OUTPUT_TAXONOMY_PLOTS_DIR>"
OUTPUT_PREFIX <- "<OUTPUT_FILE_PREFIX>"  # e.g., "circular_contigs" or "DA_contigs"

# Filtering options
MIN_COUNT <- 2           # Show only taxa with >2 contigs
LABEL_QUANTILE <- 0.5    # Label top 50% by count

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Helper Functions ----

# Clean GTDB taxonomy prefixes
clean_taxonomy <- function(tax_string) {
  if (is.na(tax_string) || tax_string == "" || tax_string %in% c("g__", "s__")) {
    return("Unclassified")
  }
  gsub("^[a-z]__", "", tax_string)
}

# Prepare data for each taxonomic level
prepare_taxonomy_data <- function(df, tax_level) {
  df %>%
    mutate(
      tax_clean = sapply(get(tax_level), clean_taxonomy),
      phylum_clean = sapply(phylum, clean_taxonomy)
    ) %>%
    group_by(tax_clean, phylum_clean) %>%
    summarise(count = n(), .groups = 'drop') %>%
    arrange(desc(count)) %>%
    mutate(
      taxonomy_level = tax_level,
      percentage = round(count / sum(count) * 100, 1)
    )
}

# ---- Load Data ----
cat("Loading annotated contigs...\n")
contigs <- read_csv(ANNOTATED_CONTIGS_CSV, show_col_types = FALSE)
cat(paste("Loaded", nrow(contigs), "contigs\n"))

# ---- Prepare Taxonomy Data ----
cat("\nPreparing taxonomy data for all levels...\n")

taxonomy_levels <- c("phylum", "class", "order", "family", "genus")

all_tax_data <- lapply(taxonomy_levels, function(level) {
  cat(paste("  Processing", level, "...\n"))
  prepare_taxonomy_data(contigs, level)
}) %>% bind_rows()

# For phylum level, phylum_clean should match tax_clean
all_tax_data <- all_tax_data %>%
  mutate(phylum_clean = ifelse(taxonomy_level == "phylum", tax_clean, phylum_clean))

# Filter by count
all_tax_data <- all_tax_data %>%
  filter(count > MIN_COUNT)

cat(paste("Filtered to", nrow(all_tax_data), "taxonomy entries\n"))

# ---- Add Empty Bars Between Levels ----
empty_bar <- 2
combined_data <- all_tax_data %>%
  group_by(taxonomy_level) %>%
  group_split() %>%
  map_df(function(group_df) {
    to_add <- data.frame(
      tax_clean = rep(NA, empty_bar),
      phylum_clean = rep(NA, empty_bar),
      count = rep(0, empty_bar),
      taxonomy_level = rep(unique(group_df$taxonomy_level), empty_bar),
      percentage = rep(NA, empty_bar)
    )
    bind_rows(group_df, to_add)
  })

# ---- Arrange Data ----
level_order <- c("phylum", "class", "order", "family", "genus")
combined_data$taxonomy_level <- factor(combined_data$taxonomy_level, levels = level_order)
combined_data <- combined_data %>% arrange(taxonomy_level, desc(count))
combined_data$id <- seq(1, nrow(combined_data))

# ---- Prepare Label Data ----
label_data <- combined_data %>%
  filter(!is.na(tax_clean) & count > 0)

number_of_bar <- nrow(combined_data)
label_data$angle <- 90 - 360 * (label_data$id - 0.5) / number_of_bar
label_data$hjust <- ifelse(label_data$angle < -90, 1, 0)
label_data$angle <- ifelse(label_data$angle < -90, 
                           label_data$angle + 180, 
                           label_data$angle)

# ---- Prepare Base Lines ----
base_data <- combined_data %>%
  group_by(taxonomy_level) %>%
  summarise(
    start = min(id),
    end = max(id) - empty_bar,
    .groups = 'drop'
  ) %>%
  rowwise() %>%
  mutate(title = mean(c(start, end)))

# ---- Prepare Grid ----
grid_data <- base_data
grid_data$end <- grid_data$end[c(nrow(grid_data), 1:(nrow(grid_data)-1))] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1, ]

# ---- Calculate Scale ----
max_count <- max(combined_data$count, na.rm = TRUE)
scale_values <- pretty(c(0, max_count), n = 4)
scale_values <- scale_values[scale_values > 0 & scale_values <= max_count]

# ---- Create Color Palettes ----
unique_phyla <- unique(combined_data$phylum_clean[!is.na(combined_data$phylum_clean)])
n_phyla <- length(unique_phyla)

if (n_phyla <= 12) {
  phylum_colors <- setNames(brewer.pal(max(3, n_phyla), "Set3"), unique_phyla)
} else {
  phylum_colors <- setNames(rainbow(n_phyla), unique_phyla)
}

tax_colors <- c(
  "phylum" = "#E41A1C",
  "class" = "#377EB8", 
  "order" = "#4DAF4A",
  "family" = "#984EA3",
  "genus" = "#FF7F00"
)

# ---- Create Circular Barplot ----
cat("\nGenerating circular barplot...\n")

p <- ggplot(combined_data, aes(x = id, y = count, fill = phylum_clean)) +
  
  # Grid lines
  {if(length(scale_values) > 0) {
    lapply(scale_values, function(val) {
      geom_segment(data = grid_data, 
                   aes(x = end, y = val, xend = start, yend = val),
                   colour = "grey80", alpha = 1, linewidth = 0.3, 
                   inherit.aes = FALSE)
    })
  }} +
  
  # Scale labels
  {if(length(scale_values) > 0) {
    annotate("text", 
             x = rep(max(combined_data$id), length(scale_values)), 
             y = scale_values,
             label = as.character(scale_values),
             color = "grey40", size = 3, angle = 0, 
             fontface = "bold", hjust = 1)
  }} +
  
  # Bars
  geom_bar(stat = "identity", alpha = 0.7) +
  
  # Scale
  ylim(-max_count * 1.5, max_count * 1.5) +
  
  # Colors
  scale_fill_manual(values = phylum_colors,
                    name = "Phylum",
                    na.value = "grey90") +
  
  # Theme
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  
  # Polar coordinates
  coord_polar(start = 0) +
  
  # Labels
  geom_text(
    data = filter(label_data, count >= quantile(count, LABEL_QUANTILE, na.rm = TRUE)),
    aes(x = id, y = count + max_count * 0.05, 
        label = paste0(tax_clean, " (", count, ")"),
        hjust = hjust),
    color = "black", fontface = "bold", alpha = 0.7, 
    size = 2.9, 
    angle = filter(label_data, count >= quantile(count, LABEL_QUANTILE, na.rm = TRUE))$angle,
    inherit.aes = FALSE
  ) +
  
  # Base lines
  geom_segment(
    data = base_data,
    aes(x = start, y = -max_count * 0.1, xend = end, yend = -max_count * 0.1),
    colour = "black", alpha = 0.8, linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  
  # Level labels
  geom_text(
    data = base_data,
    aes(x = title, y = -max_count * 0.35, label = toupper(taxonomy_level)),
    colour = tax_colors[base_data$taxonomy_level], 
    alpha = 0.9, size = 3, fontface = "bold", 
    inherit.aes = FALSE
  )

# Save plot
output_file <- file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_taxonomy_circular.png"))
ggsave(output_file, p, width = 10, height = 9, dpi = 1000)
cat(paste("\nSaved plot to:", output_file, "\n"))

# ---- Create Summary Table ----
taxonomy_summary <- all_tax_data %>%
  group_by(taxonomy_level) %>%
  summarise(
    unique_taxa = sum(!is.na(tax_clean) & tax_clean != "Unclassified"),
    total_contigs = sum(count, na.rm = TRUE),
    .groups = 'drop'
  )

cat("\nTaxonomy Summary:\n")
print(taxonomy_summary)

summary_file <- file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_taxonomy_summary.csv"))
write_csv(taxonomy_summary, summary_file)
cat(paste("Saved summary to:", summary_file, "\n"))

cat("\nTaxonomy visualization complete!\n")