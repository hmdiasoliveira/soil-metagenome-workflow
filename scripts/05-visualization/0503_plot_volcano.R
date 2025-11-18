#!/usr/bin/env Rscript
# ============================================================================
# Volcano Plot Functions
# ============================================================================
# Description: Functions for creating volcano plots from differential abundance results
# Usage: Source this file and call plot_volcano() or plot_volcano_no_title()
# ============================================================================

library(EnhancedVolcano)
library(ggplot2)
library(dplyr)

# ---- Helper Function ----

#' Pick the first available column from a list of options
#' @param df Data frame
#' @param cols Character vector of column names to try
#' @return Column name if found, NA otherwise
.pick_col <- function(df, cols) {
  for (col in cols) {
    if (col %in% names(df)) return(col)
  }
  return(NA_character_)
}

# ---- Main Volcano Plot Function ----

#' Create volcano plot with customizable FDR control
#'
#' @param df Data frame with differential abundance results
#' @param title Plot title
#' @param cutoff_mode Either "BH" for Benjamini-Hochberg or "empirical" for empirical FDR
#' @param fdr_cutoff FDR cutoff for BH mode (default 0.05)
#' @param tau P-value threshold for empirical mode
#' @param logfc_cutoff Log fold change threshold (default 1)
#' @param label_n_genes Number of top genes to label (default 0)
#' @return ggplot object
plot_volcano <- function(df,
                         title = NULL,
                         cutoff_mode = c("BH", "empirical"),
                         fdr_cutoff = 0.05,
                         tau = NULL,
                         logfc_cutoff = 1,
                         label_n_genes = 0) {
  
  cutoff_mode <- match.arg(cutoff_mode)
  df <- as.data.frame(df)
  
  # Identify p-value and FDR columns
  p_col <- .pick_col(df, c("PValue", "P.Value", "p", "p_value"))
  fdr_col <- .pick_col(df, c("FDR", "adj.P.Val", "qvalue"))
  
  if (cutoff_mode == "BH" && is.na(fdr_col)) {
    stop("BH mode requires an FDR/adj.P.Val column.")
  }
  if (cutoff_mode == "empirical" && is.na(p_col)) {
    stop("Empirical mode requires a raw p-value column (e.g., PValue).")
  }
  
  # Ensure MAG column exists for labeling
  if (!"MAG" %in% names(df)) {
    df$MAG <- rownames(df)
  }
  
  # Determine significance cutoff
  if (cutoff_mode == "BH") {
    cutoff_p_value <- fdr_cutoff
    cutoff_col <- fdr_col
    subtitle_text <- "BH/FDR control"
  } else {
    if (is.null(tau)) stop("Provide 'tau' when cutoff_mode='empirical'.")
    cutoff_p_value <- tau
    cutoff_col <- p_col
    subtitle_text <- bquote(paste("Empirical ", tau, " control (", tau, " ≈ ", .(signif(tau, 3)), ")"))
  }
  
  # Mark significant contigs
  df$is_significant <- df[[cutoff_col]] <= cutoff_p_value & abs(df$logFC) >= logfc_cutoff
  
  # Select top labels
  selectLab <- NULL
  if (label_n_genes > 0 && any(df$is_significant, na.rm = TRUE)) {
    sig_df <- df %>% filter(is_significant) %>% arrange(!!sym(cutoff_col))
    selectLab <- head(sig_df$MAG, label_n_genes)
  }
  
  # Create volcano plot
  p <- EnhancedVolcano::EnhancedVolcano(
    toptable = df,
    lab = df$MAG,
    x = "logFC",
    y = p_col,
    selectLab = selectLab,
    title = title,
    subtitle = subtitle_text,
    pCutoff = cutoff_p_value,
    pCutoffCol = cutoff_col,
    FCcutoff = logfc_cutoff,
    cutoffLineType = "longdash",
    cutoffLineCol = "black",
    cutoffLineWidth = 0.4,
    pointSize = 1.8,
    labSize = 3.0,
    col = c("grey30", "forestgreen", "lightblue", "pink"),
    colAlpha = 0.8,
    legendPosition = "top",
    legendLabels = c(
      "Not significant",
      "Large effect size",
      "P-value significant",
      "P-value and effect size significant"
    ),
    legendLabSize = 8,
    legendIconSize = 5,
    drawConnectors = label_n_genes > 0,
    max.overlaps = Inf,
    widthConnectors = 0.5
  )
  
  # Add cutoff label
  cutoff_label <- if (cutoff_mode == "BH") {
    paste0("FDR cutoff = ", fdr_cutoff)
  } else {
    paste0("tau = ", signif(tau, 3))
  }
  
  p <- p + ggplot2::annotate("text",
                             x = Inf, y = -log10(cutoff_p_value),
                             label = cutoff_label,
                             hjust = 1.1, vjust = -0.5,
                             size = 3.5, fontface = "italic")
  
  return(p)
}

# ---- Volcano Plot Without Titles ----

#' Create volcano plot without title/subtitle (for publication)
#'
#' @param df Data frame with differential abundance results
#' @param title Plot title (optional)
#' @param cutoff_mode Either "BH" or "empirical"
#' @param fdr_cutoff FDR cutoff for BH mode
#' @param tau P-value threshold for empirical mode
#' @param logfc_cutoff Log fold change threshold
#' @param label_n_genes Number of top genes to label
#' @param show_titles Whether to show title/subtitle
#' @return ggplot object
plot_volcano_no_title <- function(df,
                                   title = NULL,
                                   cutoff_mode = c("BH", "empirical"),
                                   fdr_cutoff = 0.05,
                                   tau = NULL,
                                   logfc_cutoff = 1,
                                   label_n_genes = 20,
                                   show_titles = FALSE) {
  
  cutoff_mode <- match.arg(cutoff_mode)
  df <- as.data.frame(df)
  
  p_col <- .pick_col(df, c("PValue", "P.Value", "p", "p_value"))
  fdr_col <- .pick_col(df, c("FDR", "adj.P.Val", "qvalue"))
  
  if (cutoff_mode == "BH" && is.na(fdr_col)) stop("BH mode requires an FDR/adj.P.Val column.")
  if (cutoff_mode == "empirical" && is.na(p_col)) stop("Empirical mode requires a raw p-value column.")
  
  if (!"MAG" %in% names(df)) df$MAG <- rownames(df)
  
  # Determine cutoff
  if (cutoff_mode == "BH") {
    cutoff_p_value <- fdr_cutoff
    cutoff_col <- fdr_col
    subtitle_text <- "BH/FDR control"
  } else {
    if (is.null(tau)) stop("Provide 'tau' when cutoff_mode='empirical'.")
    cutoff_p_value <- tau
    cutoff_col <- p_col
    subtitle_text <- bquote(paste("Empirical ", tau, " control"))
  }
  
  df$is_significant <- df[[cutoff_col]] <= cutoff_p_value & abs(df$logFC) >= logfc_cutoff
  
  # Create custom color scheme
  keyvals <- ifelse(
    df$logFC < -logfc_cutoff & df[[cutoff_col]] < cutoff_p_value, 'royalblue3',
    ifelse(df$logFC > logfc_cutoff & df[[cutoff_col]] < cutoff_p_value, 'red2',
           ifelse(abs(df$logFC) > logfc_cutoff & df[[cutoff_col]] >= cutoff_p_value, 'grey50', 'grey80')))
  
  keyvals[is.na(keyvals)] <- 'grey80'
  names(keyvals)[keyvals == 'red2'] <- 'Up'
  names(keyvals)[keyvals == 'royalblue3'] <- 'Down'
  names(keyvals)[keyvals == 'grey50'] <- 'Partial'
  names(keyvals)[keyvals == 'grey80'] <- 'NS'
  
  # Select labels
  selectLab <- NULL
  if (label_n_genes > 0 && any(df$is_significant, na.rm = TRUE)) {
    sig_df <- df %>% filter(is_significant) %>% arrange(!!sym(cutoff_col))
    selectLab <- head(sig_df$MAG, label_n_genes)
  }
  
  # Set titles
  title_arg <- if (show_titles) title else NULL
  subtitle_arg <- if (show_titles) subtitle_text else NULL
  
  # Create plot
  p <- EnhancedVolcano::EnhancedVolcano(
    toptable = df,
    lab = df$MAG,
    x = "logFC",
    y = p_col,
    selectLab = selectLab,
    title = title_arg,
    subtitle = subtitle_arg,
    pCutoff = cutoff_p_value,
    pCutoffCol = cutoff_col,
    FCcutoff = logfc_cutoff,
    cutoffLineType = "longdash",
    cutoffLineCol = "black",
    cutoffLineWidth = 0.4,
    pointSize = 1,
    labSize = 3.0,
    axisLabSize = 10,
    legendLabSize = 8,
    legendIconSize = 3,
    colCustom = keyvals,
    colAlpha = 0.7,
    legendPosition = "top",
    legendLabels = c("Not significant", "Log2FC or p-value only", "Down", "Up"),
    drawConnectors = label_n_genes > 0,
    max.overlaps = Inf,
    widthConnectors = 0.5,
    colConnectors = "grey30"
  )
  
  # Remove titles explicitly
  p <- p + ggplot2::labs(title = NULL, subtitle = NULL) +
    ggplot2::theme(plot.title = ggplot2::element_blank(),
                   plot.subtitle = ggplot2::element_blank())
  
  # Add cutoff label
  cutoff_label <- if (cutoff_mode == "BH") {
    paste0("FDR cutoff = ", fdr_cutoff)
  } else {
    paste0("tau = ", signif(tau, 3))
  }
  
  p <- p + ggplot2::annotate("text",
                             x = Inf, y = -log10(cutoff_p_value),
                             label = cutoff_label,
                             hjust = 1.1, vjust = -0.5,
                             size = 3.5, fontface = "italic")
  
  return(p)
}