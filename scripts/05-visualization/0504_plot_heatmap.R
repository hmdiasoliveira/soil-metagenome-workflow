#!/usr/bin/env Rscript
# ============================================================================
# Heatmap Plot Functions
# ============================================================================
# Description: Functions for creating heatmaps of significant contigs
# Usage: Source this file and call plot_heatmap_significant()
# ============================================================================

library(pheatmap)
library(edgeR)
library(dplyr)
library(rlang)

# ---- Helper Function ----

#' Pick the first available column from a list of options
.pick_col <- function(df, cols) {
  for (col in cols) {
    if (col %in% names(df)) return(col)
  }
  return(NA_character_)
}

# ---- Main Heatmap Function ----

#' Plot heatmap of significant contigs
#'
#' @param sig_df Data frame with differential abundance results
#' @param dge_object DGEList object from edgeR
#' @param metadata_df Metadata data frame (samples as rows)
#' @param col_annotation_cols Columns from metadata to annotate
#' @param title Plot title
#' @param cutoff_mode Either "BH" or "empirical"
#' @param fdr_cutoff FDR cutoff for BH mode
#' @param tau P-value threshold for empirical mode
#' @param logfc_cutoff Log fold change threshold
#' @param show_col_names Show column (sample) names
#' @param clustering_distance_rows Distance metric for rows
#' @param clustering_distance_cols Distance metric for columns
#' @param clustering_method Clustering method
#' @param filename Output filename (PNG)
#' @param matrix_outfile Optional CSV output of matrix
#' @return NULL (saves plot to file)
plot_heatmap_significant <- function(sig_df, dge_object, metadata_df,
                                      col_annotation_cols, title = NULL,
                                      cutoff_mode = c("BH", "empirical"),
                                      fdr_cutoff = 0.05, tau = NULL, logfc_cutoff = 1,
                                      show_col_names = FALSE,
                                      clustering_distance_rows = "correlation",
                                      clustering_distance_cols = "euclidean",
                                      clustering_method = "ward.D2",
                                      filename = NULL,
                                      matrix_outfile = NULL) {
  
  cutoff_mode <- match.arg(cutoff_mode)
  if (!is.data.frame(sig_df)) sig_df <- as.data.frame(sig_df)
  if (!"MAG" %in% names(sig_df)) sig_df$MAG <- rownames(sig_df)
  if (!"logFC" %in% names(sig_df)) stop("Need 'logFC' in sig_df.")
  
  p_col <- .pick_col(sig_df, c("PValue", "P.Value", "p", "p_value"))
  fdr_col <- .pick_col(sig_df, c("FDR", "adj.P.Val", "qvalue"))
  
  if (cutoff_mode == "empirical" && is.na(p_col)) {
    stop("Empirical mode requires a raw p-value column in sig_df.")
  }
  if (cutoff_mode == "BH" && is.na(fdr_col)) {
    stop("BH mode requires an FDR/adj.P.Val column in sig_df.")
  }
  
  # Filter contigs
  if (cutoff_mode == "BH") {
    keep_mag <- sig_df %>%
      filter(.data[[fdr_col]] < fdr_cutoff & abs(logFC) > logfc_cutoff) %>%
      pull(MAG)
  } else {
    if (is.null(tau)) stop("Provide 'tau' when cutoff_mode='empirical'.")
    keep_mag <- sig_df %>%
      filter(.data[[p_col]] <= tau & abs(logFC) > logfc_cutoff) %>%
      pull(MAG)
  }
  
  keep_mag <- intersect(keep_mag, rownames(dge_object$counts))
  
  if (length(keep_mag) == 0) {
    message("No significant MAGs found. No heatmap generated.")
    return(invisible(NULL))
  }
  
  # Get log-CPM matrix
  logcpm_matrix <- cpm(dge_object, log = TRUE, prior.count = 1)
  heatmap_data <- logcpm_matrix[keep_mag, , drop = FALSE]
  
  # Prepare annotations
  cols_available <- intersect(col_annotation_cols, names(metadata_df))
  annotation_col <- metadata_df[colnames(heatmap_data), cols_available, drop = FALSE]
  
  if (!identical(rownames(annotation_col), colnames(heatmap_data))) {
    stop("Metadata and heatmap matrix samples do not align.")
  }
  
  # Save matrix if requested
  if (!is.null(matrix_outfile)) {
    write.csv(heatmap_data, file = matrix_outfile)
  }
  
  # Create title
  if (cutoff_mode == "BH") {
    main_title <- sprintf("%s [BH FDR≤%.02f; |logFC|>%g]", 
                          ifelse(is.null(title), "", title), fdr_cutoff, logfc_cutoff)
  } else {
    main_title <- sprintf("%s [τ≈%g; |logFC|>%g]", 
                          ifelse(is.null(title), "", title), signif(tau, 3), logfc_cutoff)
  }
  
  # Plot heatmap
  pheatmap::pheatmap(
    heatmap_data,
    scale = "row",
    cluster_rows = TRUE, cluster_cols = TRUE,
    clustering_distance_rows = clustering_distance_rows,
    clustering_distance_cols = clustering_distance_cols,
    clustering_method = clustering_method,
    show_rownames = TRUE,
    show_colnames = show_col_names,
    annotation_col = annotation_col,
    main = main_title,
    fontsize_row = 6, fontsize_col = 8,
    filename = filename,
    width = 10, height = max(7, nrow(heatmap_data) * 0.05 + 2)
  )
  
  invisible()
}

# ---- Heatmap Without Title ----

#' Plot heatmap without title (for publication)
#'
#' @param sig_df Data frame with differential abundance results
#' @param dge_object DGEList object from edgeR
#' @param metadata_df Metadata data frame
#' @param col_annotation_cols Columns to annotate
#' @param title Plot title (optional)
#' @param cutoff_mode Either "BH" or "empirical"
#' @param fdr_cutoff FDR cutoff for BH mode
#' @param tau P-value threshold for empirical mode
#' @param logfc_cutoff Log fold change threshold
#' @param show_col_names Show column names
#' @param filename Output filename
#' @param matrix_outfile Optional CSV output
#' @param show_title Whether to show title
#' @return NULL
plot_heatmap_significant_no_title <- function(sig_df, dge_object, metadata_df,
                                               col_annotation_cols, title = NULL,
                                               cutoff_mode = c("BH", "empirical"),
                                               fdr_cutoff = 0.05, tau = NULL, logfc_cutoff = 1,
                                               show_col_names = FALSE,
                                               filename = NULL, matrix_outfile = NULL,
                                               show_title = FALSE) {
  
  cutoff_mode <- match.arg(cutoff_mode)
  if (!is.data.frame(sig_df)) sig_df <- as.data.frame(sig_df)
  if (!"MAG" %in% names(sig_df)) sig_df$MAG <- rownames(sig_df)
  if (!"logFC" %in% names(sig_df)) stop("Need 'logFC' in sig_df.")
  
  p_col <- .pick_col(sig_df, c("PValue", "P.Value", "p", "p_value"))
  fdr_col <- .pick_col(sig_df, c("FDR", "adj.P.Val", "qvalue"))
  
  if (cutoff_mode == "empirical" && is.na(p_col)) stop("Empirical mode requires a raw p-value column.")
  if (cutoff_mode == "BH" && is.na(fdr_col)) stop("BH mode requires an FDR/adj.P.Val column.")
  
  # Filter contigs
  if (cutoff_mode == "BH") {
    keep_mag <- sig_df %>% filter(.data[[fdr_col]] < fdr_cutoff & abs(logFC) > logfc_cutoff) %>% pull(MAG)
  } else {
    if (is.null(tau)) stop("Provide 'tau' when cutoff_mode='empirical'.")
    keep_mag <- sig_df %>% filter(.data[[p_col]] <= tau & abs(logFC) > logfc_cutoff) %>% pull(MAG)
  }
  
  keep_mag <- intersect(keep_mag, rownames(dge_object$counts))
  if (length(keep_mag) == 0) { message("No significant MAGs found."); return(invisible(NULL)) }
  
  logcpm_matrix <- edgeR::cpm(dge_object, log = TRUE, prior.count = 1)
  heatmap_data <- logcpm_matrix[keep_mag, , drop = FALSE]
  
  cols_available <- intersect(col_annotation_cols, names(metadata_df))
  annotation_col <- metadata_df[colnames(heatmap_data), cols_available, drop = FALSE]
  if (!identical(rownames(annotation_col), colnames(heatmap_data))) stop("Metadata and matrix samples do not align.")
  
  if (!is.null(matrix_outfile)) write.csv(heatmap_data, file = matrix_outfile)
  
  # Build or suppress title
  if (show_title) {
    if (cutoff_mode == "BH") {
      main_title <- sprintf("%s [BH FDR≤%.02f; |logFC|>%g]", ifelse(is.null(title), "", title), fdr_cutoff, logfc_cutoff)
    } else {
      main_title <- sprintf("%s [τ≈%g; |logFC|>%g]", ifelse(is.null(title), "", title), signif(tau, 3), logfc_cutoff)
    }
  } else {
    main_title <- ""
  }
  
  # Use Spearman correlation distance
  row_dist <- as.dist(1 - cor(t(heatmap_data), method = "spearman", use = "pairwise.complete.obs"))
  col_dist <- as.dist(1 - cor(heatmap_data, method = "spearman", use = "pairwise.complete.obs"))
  
  pheatmap::pheatmap(
    heatmap_data,
    scale = "row",
    cluster_rows = TRUE, cluster_cols = TRUE,
    clustering_distance_rows = row_dist,
    clustering_distance_cols = col_dist,
    clustering_method = "average",
    show_rownames = TRUE,
    show_colnames = show_col_names,
    annotation_col = annotation_col,
    main = main_title,
    fontsize_row = 6, fontsize_col = 8,
    filename = filename,
    width = 10, height = max(7, nrow(heatmap_data) * 0.05 + 2)
  )
  
  invisible()
}