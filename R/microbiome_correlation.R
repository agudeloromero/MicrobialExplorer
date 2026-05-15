# =============================================================================
# MICROBIOME CORRELATION AND ASSOCIATION ANALYSIS
# =============================================================================
# Description : Comprehensive correlation and association pipeline covering
#               taxa-metadata correlations, taxa-taxa co-occurrence,
#               cross-domain associations, network correlations,
#               mixed-effects models, and mediation analysis
# Input       : Filtered phyloseq object (output from microbiome_qc.R)
# Output      : Correlation matrices, heatmaps, scatter plots,
#               association tables, bubble plots, network overlays
# Author      : Patricia
# Dependencies: phyloseq, ggplot2, dplyr, tidyr, patchwork, vegan,
#               scales, tibble, RColorBrewer, stringr, forcats,
#               Hmisc, corrplot, ggrepel, lme4, lmerTest (optional)
# =============================================================================

# --- 1. LOAD LIBRARIES -------------------------------------------------------

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(vegan)
  library(scales)
  library(tibble)
  library(RColorBrewer)
  library(stringr)
  library(forcats)
  library(Hmisc)       # rcorr for p-values on correlation matrices
  library(ggrepel)
})

pkg_available <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# Global theme
theme_microbiome <- function() {
  theme_bw() +
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = "#2c3e50", colour = NA),
      strip.text        = element_text(colour = "white", face = "bold", size = 10),
      axis.title        = element_text(size = 11),
      axis.text         = element_text(size = 9),
      legend.title      = element_text(size = 10, face = "bold"),
      legend.text       = element_text(size = 9),
      plot.title        = element_text(size = 13, face = "bold", hjust = 0),
      plot.subtitle     = element_text(size = 10, colour = "grey40", hjust = 0),
      plot.caption      = element_text(size = 8, colour = "grey60", hjust = 1)
    )
}


# =============================================================================
# SECTION 1 — TAXA-METADATA CORRELATION
# =============================================================================

#' Correlate individual taxon abundances against continuous metadata variables.
#'
#' Tests Spearman correlation between each taxon and each numeric metadata
#' variable, applies Benjamini-Hochberg correction, and returns a tidy
#' results table with significance annotations.
#'
#' @param ps              A phyloseq object (raw counts or relative abundance).
#' @param rank            Taxonomic rank to agglomerate to. Default = "Genus".
#' @param meta_vars       Character vector of numeric metadata variables.
#'                        If NULL, all numeric metadata columns are used.
#' @param transform       Abundance transformation: "clr", "relative", "log10".
#' @param method          Correlation method: "spearman" or "pearson".
#' @param p_adjust        P-value adjustment method. Default = "BH".
#' @param alpha           Significance threshold. Default = 0.05.
#' @param min_prevalence  Minimum taxon prevalence. Default = 0.10.
#' @return A tidy data frame of all taxon-metadata correlations.

correlate_taxa_metadata <- function(ps,
                                     rank           = "Genus",
                                     meta_vars      = NULL,
                                     transform      = "clr",
                                     method         = "spearman",
                                     p_adjust       = "BH",
                                     alpha          = 0.05,
                                     min_prevalence = 0.10) {

  cat("=== Taxa-metadata correlation ===\n")

  # --- Agglomerate and transform -------------------------------------------
  ps_agg   <- tax_glom(ps, taxrank = rank, NArm = FALSE)
  new_names <- as.character(tax_table(ps_agg)[, rank])
  new_names[is.na(new_names)] <- paste0("Unknown_", seq_len(sum(is.na(new_names))))
  dup_idx  <- duplicated(new_names)
  new_names <- make.unique(new_names, sep = "_dup")
  taxa_names(ps_agg) <- new_names

  otu_mat <- as.matrix(otu_table(ps_agg))
  if (!taxa_are_rows(ps_agg)) otu_mat <- t(otu_mat)

  # Prevalence filter
  prev    <- rowSums(otu_mat > 0) / ncol(otu_mat)
  otu_mat <- otu_mat[prev >= min_prevalence, ]

  # Transform
  if (transform == "clr") {
    feat_mat <- apply(otu_mat + 0.5, 2, function(x) log(x) - mean(log(x)))
  } else if (transform == "relative") {
    feat_mat <- apply(otu_mat, 2, function(x) x / sum(x))
  } else {
    feat_mat <- log10(otu_mat + 1)
  }
  feat_t <- t(feat_mat)  # samples as rows

  # --- Metadata -----------------------------------------------------------
  meta_df  <- data.frame(sample_data(ps_agg))

  if (is.null(meta_vars)) {
    meta_vars <- names(meta_df)[sapply(meta_df, is.numeric)]
    if (length(meta_vars) == 0) {
      stop("No numeric metadata variables found. Specify meta_vars manually.")
    }
    cat("  Auto-detected numeric variables:", paste(meta_vars, collapse = ", "), "\n")
  }

  meta_vars <- intersect(meta_vars, colnames(meta_df))
  meta_sub  <- meta_df[, meta_vars, drop = FALSE]

  # Align
  shared    <- intersect(rownames(feat_t), rownames(meta_sub))
  feat_t    <- feat_t[shared, , drop = FALSE]
  meta_sub  <- meta_sub[shared, , drop = FALSE]

  cat("  Taxa tested:", ncol(feat_t), "\n")
  cat("  Metadata variables:", length(meta_vars), "\n")

  # --- Pairwise correlations -----------------------------------------------
  result_list <- lapply(meta_vars, function(mv) {
    y_vals <- meta_sub[[mv]]
    complete_idx <- !is.na(y_vals)
    y_sub  <- y_vals[complete_idx]
    X_sub  <- feat_t[complete_idx, , drop = FALSE]

    taxon_results <- lapply(colnames(X_sub), function(tx) {
      x_vals <- X_sub[, tx]
      ct <- tryCatch(
        cor.test(x_vals, y_sub, method = method, exact = FALSE),
        error = function(e) NULL
      )
      if (is.null(ct)) return(NULL)
      data.frame(
        taxon       = tx,
        metadata    = mv,
        rho         = round(ct$estimate, 4),
        p_value     = ct$p.value,
        n           = sum(complete_idx),
        stringsAsFactors = FALSE
      )
    })
    bind_rows(taxon_results)
  })

  results_df <- bind_rows(result_list) %>%
    mutate(
      q_value      = p.adjust(p_value, method = p_adjust),
      significant  = q_value < alpha,
      direction    = ifelse(rho > 0, "Positive", "Negative"),
      sig_label    = case_when(
        q_value < 0.001 ~ "***",
        q_value < 0.01  ~ "**",
        q_value < 0.05  ~ "*",
        q_value < 0.1   ~ ".",
        TRUE             ~ ""
      )
    ) %>%
    arrange(q_value)

  n_sig <- sum(results_df$significant)
  cat("  Significant associations (q <", alpha, "):", n_sig, "\n\n")

  return(results_df)
}


# =============================================================================
# SECTION 2 — CORRELATION HEATMAP (TAXA × METADATA)
# =============================================================================

#' Plot a heatmap of taxa-metadata Spearman correlations.
#'
#' @param cor_results   Data frame from correlate_taxa_metadata().
#' @param top_n         Top taxa by maximum absolute correlation. Default = 30.
#' @param alpha         Significance threshold for asterisk overlay. Default = 0.05.
#' @param cluster       Whether to hierarchically cluster rows and cols. Default = TRUE.
#' @return A ggplot heatmap.

plot_taxa_metadata_heatmap <- function(cor_results,
                                        top_n   = 30,
                                        alpha   = 0.05,
                                        cluster = TRUE) {

  cat("=== Taxa-metadata correlation heatmap ===\n")

  # Pivot to wide matrix
  rho_wide <- cor_results %>%
    select(taxon, metadata, rho) %>%
    pivot_wider(names_from = metadata, values_from = rho, values_fill = 0)

  sig_wide <- cor_results %>%
    select(taxon, metadata, sig_label) %>%
    pivot_wider(names_from = metadata, values_from = sig_label, values_fill = "")

  # Select top_n taxa by max absolute correlation
  rho_mat   <- as.matrix(rho_wide[, -1])
  rownames(rho_mat) <- rho_wide$taxon

  max_cor   <- apply(abs(rho_mat), 1, max)
  top_taxa  <- names(sort(max_cor, decreasing = TRUE))[seq_len(min(top_n, nrow(rho_mat)))]
  rho_mat   <- rho_mat[top_taxa, , drop = FALSE]

  sig_mat   <- as.matrix(sig_wide[match(top_taxa, sig_wide$taxon), -1])
  rownames(sig_mat) <- top_taxa

  # Hierarchical clustering
  if (cluster && nrow(rho_mat) > 2) {
    row_order <- hclust(dist(rho_mat))$order
    col_order <- hclust(dist(t(rho_mat)))$order
    rho_mat   <- rho_mat[row_order, col_order]
    sig_mat   <- sig_mat[row_order, col_order]
  }

  # Melt
  heat_df <- as.data.frame(rho_mat) %>%
    rownames_to_column("taxon") %>%
    pivot_longer(-taxon, names_to = "metadata", values_to = "rho") %>%
    mutate(
      taxon    = factor(taxon, levels = rev(rownames(rho_mat))),
      metadata = factor(metadata, levels = colnames(rho_mat))
    )

  sig_df <- as.data.frame(sig_mat) %>%
    rownames_to_column("taxon") %>%
    pivot_longer(-taxon, names_to = "metadata", values_to = "sig") %>%
    mutate(
      taxon    = factor(taxon, levels = rev(rownames(rho_mat))),
      metadata = factor(metadata, levels = colnames(rho_mat))
    )

  heat_df <- left_join(heat_df, sig_df, by = c("taxon", "metadata"))

  p <- ggplot(heat_df, aes(x = metadata, y = taxon, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sig), size = 3.5, vjust = 0.75,
              colour = "black", fontface = "bold") +
    scale_fill_gradient2(
      low      = "#3498db",
      mid      = "white",
      high     = "#e74c3c",
      midpoint = 0,
      limits   = c(-1, 1),
      name     = "Spearman\nrho"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      title    = "Taxa-metadata correlation heatmap",
      subtitle = paste0("Top ", nrow(rho_mat), " taxa by max |rho| | * q<0.05, ** q<0.01, *** q<0.001"),
      x        = "Metadata variable",
      y        = NULL
    ) +
    theme_microbiome() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y  = element_text(face = "italic", size = 8),
      panel.border = element_rect(colour = "grey80", fill = NA)
    )

  cat("  Taxa plotted:", nrow(rho_mat),
      "| Metadata variables:", ncol(rho_mat), "\n\n")
  return(p)
}


# =============================================================================
# SECTION 3 — TAXA-METADATA SCATTER PLOTS
# =============================================================================

#' Generate scatter plots for the top significant taxa-metadata associations.
#'
#' @param ps          A phyloseq object.
#' @param cor_results Data frame from correlate_taxa_metadata().
#' @param rank        Taxonomic rank. Default = "Genus".
#' @param top_n       Number of top associations to plot. Default = 12.
#' @param transform   Abundance transformation. Default = "clr".
#' @param group_var   Optional grouping variable for point colouring.
#' @return A patchwork of scatter plots.

plot_top_associations <- function(ps,
                                   cor_results,
                                   rank      = "Genus",
                                   top_n     = 12,
                                   transform = "clr",
                                   group_var = NULL) {

  cat("=== Top association scatter plots ===\n")

  top_assoc <- cor_results %>%
    filter(significant) %>%
    arrange(q_value) %>%
    slice_head(n = top_n)

  if (nrow(top_assoc) == 0) {
    cat("  No significant associations to plot.\n\n")
    return(NULL)
  }

  # Prepare abundance data
  ps_agg   <- tax_glom(ps, taxrank = rank, NArm = FALSE)
  new_names <- as.character(tax_table(ps_agg)[, rank])
  new_names[is.na(new_names)] <- paste0("Unknown_", seq_len(sum(is.na(new_names))))
  taxa_names(ps_agg) <- new_names

  otu_mat  <- as.matrix(otu_table(ps_agg))
  if (!taxa_are_rows(ps_agg)) otu_mat <- t(otu_mat)

  if (transform == "clr") {
    feat_mat <- apply(otu_mat + 0.5, 2, function(x) log(x) - mean(log(x)))
  } else if (transform == "relative") {
    feat_mat <- apply(otu_mat, 2, function(x) x / sum(x))
  } else {
    feat_mat <- log10(otu_mat + 1)
  }
  feat_t   <- t(feat_mat)

  meta_df  <- data.frame(sample_data(ps_agg)) %>%
    rownames_to_column("sample")

  # Build plots
  plots <- lapply(seq_len(nrow(top_assoc)), function(i) {
    tx <- top_assoc$taxon[i]
    mv <- top_assoc$metadata[i]
    rho <- top_assoc$rho[i]
    qv  <- top_assoc$q_value[i]

    if (!tx %in% colnames(feat_t) || !mv %in% colnames(meta_df)) return(NULL)

    scatter_df <- data.frame(
      sample    = rownames(feat_t),
      abundance = feat_t[, tx],
      metadata  = meta_df[[mv]][match(rownames(feat_t), meta_df$sample)],
      stringsAsFactors = FALSE
    )

    # Add group variable
    if (!is.null(group_var) && group_var %in% colnames(meta_df)) {
      scatter_df$group <- meta_df[[group_var]][
        match(rownames(feat_t), meta_df$sample)
      ]
    }

    scatter_df <- scatter_df %>% filter(!is.na(metadata))

    n_groups  <- if (!is.null(group_var) && "group" %in% colnames(scatter_df))
                   n_distinct(scatter_df$group) else 1
    colours   <- brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)]

    p <- ggplot(scatter_df, aes(x = metadata, y = abundance)) +
      {
        if (!is.null(group_var) && "group" %in% colnames(scatter_df)) {
          list(
            geom_point(aes(colour = group), alpha = 0.7, size = 1.8),
            scale_colour_brewer(palette = "Set2", name = group_var)
          )
        } else {
          list(geom_point(colour = "#3498db", alpha = 0.7, size = 1.8))
        }
      } +
      geom_smooth(method = "lm", se = TRUE, colour = "#e74c3c",
                  fill = "#e74c3c", alpha = 0.12, linewidth = 0.9) +
      annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
               label = paste0("rho=", rho, "\nq=", signif(qv, 2)),
               size = 2.8, colour = "grey30", fontface = "bold") +
      labs(
        title = tx,
        x     = mv,
        y     = paste0(transform, " abundance")
      ) +
      theme_microbiome() +
      theme(
        plot.title  = element_text(size = 8, face = "italic", hjust = 0.5),
        axis.title  = element_text(size = 7),
        axis.text   = element_text(size = 7),
        legend.position = "none"
      )

    p
  })

  plots <- Filter(Negate(is.null), plots)
  n_cols <- min(4, length(plots))

  combined <- wrap_plots(plots, ncol = n_cols) +
    plot_annotation(
      title    = paste0("Top ", length(plots), " taxa-metadata associations"),
      subtitle = paste0("Spearman correlation | BH-corrected q < 0.05"),
      theme    = theme(
        plot.title    = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey40")
      )
    )

  cat("  Scatter plots produced:", length(plots), "\n\n")
  return(combined)
}


# =============================================================================
# SECTION 4 — TAXA-TAXA CORRELATION MATRIX
# =============================================================================

#' Compute and visualise pairwise Spearman correlations between taxa.
#'
#' @param ps            A phyloseq object.
#' @param rank          Taxonomic rank. Default = "Genus".
#' @param top_n         Top taxa to include by mean abundance. Default = 40.
#' @param transform     Transformation: "clr", "relative". Default = "clr".
#' @param alpha         Significance threshold for masking. Default = 0.05.
#' @param cluster       Hierarchically cluster the matrix. Default = TRUE.
#' @param show_sig_only Whether to show only significant correlations. Default = FALSE.
#' @return A list: correlation matrix, p-value matrix, and plot.

compute_taxa_taxa_correlation <- function(ps,
                                           rank          = "Genus",
                                           top_n         = 40,
                                           transform     = "clr",
                                           alpha         = 0.05,
                                           cluster       = TRUE,
                                           show_sig_only = FALSE) {

  cat("=== Taxa-taxa correlation matrix ===\n")

  ps_agg   <- tax_glom(ps, taxrank = rank, NArm = FALSE)
  new_names <- as.character(tax_table(ps_agg)[, rank])
  new_names[is.na(new_names)] <- paste0("Unknown_", seq_len(sum(is.na(new_names))))
  dup_idx  <- duplicated(new_names)
  new_names <- make.unique(new_names, sep = "_dup")
  taxa_names(ps_agg) <- new_names

  otu_mat  <- as.matrix(otu_table(ps_agg))
  if (!taxa_are_rows(ps_agg)) otu_mat <- t(otu_mat)

  # Select top_n taxa by mean relative abundance
  ps_rel   <- transform_sample_counts(ps_agg, function(x) x / sum(x))
  rel_mat  <- as.matrix(otu_table(ps_rel))
  if (!taxa_are_rows(ps_rel)) rel_mat <- t(rel_mat)
  mean_abund <- rowMeans(rel_mat)
  top_taxa <- names(sort(mean_abund, decreasing = TRUE))[seq_len(min(top_n, length(mean_abund)))]
  otu_sub  <- otu_mat[top_taxa, ]

  # Transform
  if (transform == "clr") {
    feat_mat <- t(apply(otu_sub + 0.5, 2, function(x) log(x) - mean(log(x))))
  } else {
    feat_mat <- t(apply(otu_sub, 2, function(x) x / sum(x)))
  }

  cat("  Computing pairwise correlations for", ncol(feat_mat), "taxa...\n")

  # rcorr for correlation + p-values simultaneously
  rc_result <- rcorr(feat_mat, type = "spearman")
  cor_mat   <- rc_result$r
  p_mat     <- rc_result$P
  diag(p_mat) <- 1

  # BH correction
  upper_idx <- upper.tri(p_mat)
  p_adj     <- p_mat
  p_adj[upper_idx] <- p.adjust(p_mat[upper_idx], method = "BH")
  p_adj[lower.tri(p_adj)] <- t(p_adj)[lower.tri(p_adj)]

  # Cluster
  if (cluster) {
    hc        <- hclust(as.dist(1 - cor_mat), method = "ward.D2")
    tax_order <- hc$labels[hc$order]
    cor_mat   <- cor_mat[tax_order, tax_order]
    p_adj     <- p_adj[tax_order, tax_order]
  }

  # Melt
  heat_df <- as.data.frame(cor_mat) %>%
    rownames_to_column("taxon1") %>%
    pivot_longer(-taxon1, names_to = "taxon2", values_to = "rho") %>%
    mutate(
      taxon1 = factor(taxon1, levels = rownames(cor_mat)),
      taxon2 = factor(taxon2, levels = rev(rownames(cor_mat)))
    )

  padj_df <- as.data.frame(p_adj) %>%
    rownames_to_column("taxon1") %>%
    pivot_longer(-taxon1, names_to = "taxon2", values_to = "p_adj")

  heat_df <- left_join(heat_df, padj_df, by = c("taxon1", "taxon2")) %>%
    mutate(
      rho_plot = if (show_sig_only) ifelse(p_adj < alpha, rho, NA) else rho,
      sig      = ifelse(p_adj < alpha & taxon1 != taxon2, "*", "")
    )

  n_sig_pairs <- sum(p_adj[upper.tri(p_adj)] < alpha)
  cat("  Significant pairs (q <", alpha, "):", n_sig_pairs, "\n\n")

  p <- ggplot(heat_df, aes(x = taxon1, y = taxon2, fill = rho_plot)) +
    geom_tile(colour = "white", linewidth = 0.1) +
    scale_fill_gradient2(
      low      = "#3498db",
      mid      = "white",
      high     = "#e74c3c",
      midpoint = 0,
      limits   = c(-1, 1),
      na.value = "grey95",
      name     = "Spearman\nrho"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      title    = paste0("Taxa-taxa correlation matrix — ", rank, " level"),
      subtitle = paste0("Top ", ncol(feat_mat), " taxa by mean abundance | ",
                        n_sig_pairs, " significant pairs (q < ", alpha, ")"),
      x        = NULL, y = NULL
    ) +
    theme_microbiome() +
    theme(
      axis.text.x  = element_text(angle = 90, hjust = 1,
                                   vjust = 0.5, face = "italic", size = 7),
      axis.text.y  = element_text(face = "italic", size = 7),
      panel.border = element_rect(colour = "grey80", fill = NA)
    )

  return(list(
    cor_matrix = cor_mat,
    p_matrix   = p_adj,
    plot       = p,
    n_sig      = n_sig_pairs
  ))
}


# =============================================================================
# SECTION 5 — BUBBLE PLOT: TOP ASSOCIATIONS OVERVIEW
# =============================================================================

#' Bubble plot summarising top taxa-metadata associations.
#'
#' Each bubble represents one taxon-metadata pair. Size = |rho|,
#' colour = direction, opacity = significance.
#'
#' @param cor_results   Data frame from correlate_taxa_metadata().
#' @param top_n         Top associations to show. Default = 50.
#' @param alpha         Significance threshold. Default = 0.05.
#' @return A ggplot bubble chart.

plot_association_bubbles <- function(cor_results,
                                      top_n = 50,
                                      alpha = 0.05) {

  cat("=== Association bubble plot ===\n")

  plot_df <- cor_results %>%
    arrange(q_value) %>%
    slice_head(n = top_n) %>%
    mutate(
      neg_log_q = -log10(q_value + 1e-10),
      direction = ifelse(rho > 0, "Positive", "Negative"),
      taxon     = fct_reorder(taxon, abs(rho), .desc = TRUE),
      alpha_pt  = ifelse(significant, 0.9, 0.4)
    )

  colour_map <- c("Positive" = "#e74c3c", "Negative" = "#3498db")

  p <- ggplot(plot_df,
              aes(x = metadata, y = taxon,
                  size    = abs(rho),
                  colour  = direction,
                  alpha   = alpha_pt)) +
    geom_point() +
    geom_text(
      data = plot_df %>% filter(significant),
      aes(label = sig_label),
      size = 2.8, vjust = -1.2, colour = "black", fontface = "bold"
    ) +
    scale_size(range = c(1, 8), name = "|rho|") +
    scale_colour_manual(values = colour_map, name = "Direction") +
    scale_alpha_identity() +
    labs(
      title    = "Taxa-metadata association overview",
      subtitle = paste0("Top ", nrow(plot_df), " associations | Opaque = q < ", alpha),
      x        = "Metadata variable",
      y        = NULL,
      caption  = "* q<0.05  ** q<0.01  *** q<0.001"
    ) +
    theme_microbiome() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(face = "italic", size = 8),
      panel.grid.major = element_line(colour = "grey92")
    )

  cat("  Associations plotted:", nrow(plot_df), "\n\n")
  return(p)
}


# =============================================================================
# SECTION 6 — METADATA-METADATA CORRELATION
# =============================================================================

#' Compute and visualise correlations between all numeric metadata variables.
#'
#' @param ps          A phyloseq object.
#' @param meta_vars   Character vector of metadata variables. Default = all numeric.
#' @param method      Correlation method. Default = "spearman".
#' @param alpha       Significance threshold. Default = 0.05.
#' @return A list: correlation matrix and plot.

plot_metadata_correlations <- function(ps,
                                        meta_vars = NULL,
                                        method    = "spearman",
                                        alpha     = 0.05) {

  cat("=== Metadata-metadata correlation matrix ===\n")

  meta_df <- data.frame(sample_data(ps))

  if (is.null(meta_vars)) {
    meta_vars <- names(meta_df)[sapply(meta_df, is.numeric)]
  }
  meta_vars <- intersect(meta_vars, colnames(meta_df))

  if (length(meta_vars) < 2) {
    cat("  Need at least 2 numeric metadata variables.\n\n")
    return(NULL)
  }

  meta_sub  <- meta_df[, meta_vars, drop = FALSE]
  meta_sub  <- meta_sub[complete.cases(meta_sub), ]

  rc        <- rcorr(as.matrix(meta_sub), type = method)
  cor_mat   <- rc$r
  p_mat     <- rc$P
  diag(p_mat) <- 1

  # BH correction
  upper_idx <- upper.tri(p_mat)
  p_adj     <- p_mat
  p_adj[upper_idx] <- p.adjust(p_mat[upper_idx], method = "BH")
  p_adj[lower.tri(p_adj)] <- t(p_adj)[lower.tri(p_adj)]

  # Cluster
  if (nrow(cor_mat) > 2) {
    hc_order  <- hclust(as.dist(1 - abs(cor_mat)))$order
    cor_mat   <- cor_mat[hc_order, hc_order]
    p_adj     <- p_adj[hc_order, hc_order]
  }

  # Melt
  heat_df <- as.data.frame(cor_mat) %>%
    rownames_to_column("var1") %>%
    pivot_longer(-var1, names_to = "var2", values_to = "rho") %>%
    left_join(
      as.data.frame(p_adj) %>%
        rownames_to_column("var1") %>%
        pivot_longer(-var1, names_to = "var2", values_to = "p_adj"),
      by = c("var1", "var2")
    ) %>%
    mutate(
      var1     = factor(var1, levels = rownames(cor_mat)),
      var2     = factor(var2, levels = rev(rownames(cor_mat))),
      sig_lab  = case_when(
        p_adj < 0.001 & var1 != var2 ~ "***",
        p_adj < 0.01  & var1 != var2 ~ "**",
        p_adj < 0.05  & var1 != var2 ~ "*",
        TRUE                          ~ ""
      )
    )

  cat("  Variables:", length(meta_vars), "\n")
  cat("  Significant pairs:",
      sum(p_adj[upper.tri(p_adj)] < alpha), "\n\n")

  p <- ggplot(heat_df, aes(x = var1, y = var2, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = paste0(round(rho, 2), sig_lab)),
              size = 3, colour = ifelse(abs(heat_df$rho) > 0.6, "white", "black")) +
    scale_fill_gradient2(
      low = "#3498db", mid = "white", high = "#e74c3c",
      midpoint = 0, limits = c(-1, 1), name = "rho"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      title    = "Metadata variable correlations",
      subtitle = paste0(length(meta_vars), " variables | * q<0.05 ** q<0.01 *** q<0.001"),
      x = NULL, y = NULL
    ) +
    theme_microbiome() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      panel.border = element_rect(colour = "grey80", fill = NA)
    )

  return(list(cor_matrix = cor_mat, p_matrix = p_adj, plot = p))
}


# =============================================================================
# SECTION 7 — CROSS-DOMAIN ASSOCIATION (MICROBIOME × METABOLOMICS)
# =============================================================================

#' Compute cross-domain correlations between microbiome taxa and
#' external data (metabolites, cytokines, clinical markers).
#'
#' @param ps            A phyloseq object.
#' @param external_mat  Matrix of external features (features × samples).
#' @param rank          Taxonomic rank. Default = "Genus".
#' @param transform     Microbiome transformation. Default = "clr".
#' @param method        Correlation method. Default = "spearman".
#' @param p_adjust      P-value adjustment method. Default = "BH".
#' @param alpha         Significance threshold. Default = 0.05.
#' @param top_n_taxa    Top taxa to include. Default = 30.
#' @param top_n_ext     Top external features to include. Default = 30.
#' @return A list: results table and heatmap.

correlate_cross_domain <- function(ps,
                                    external_mat,
                                    rank        = "Genus",
                                    transform   = "clr",
                                    method      = "spearman",
                                    p_adjust    = "BH",
                                    alpha       = 0.05,
                                    top_n_taxa  = 30,
                                    top_n_ext   = 30) {

  cat("=== Cross-domain correlation (microbiome × external) ===\n")

  # Prepare microbiome features
  ps_agg   <- tax_glom(ps, taxrank = rank, NArm = FALSE)
  new_names <- as.character(tax_table(ps_agg)[, rank])
  new_names[is.na(new_names)] <- paste0("Unknown_", seq_len(sum(is.na(new_names))))
  taxa_names(ps_agg) <- new_names

  otu_mat  <- as.matrix(otu_table(ps_agg))
  if (!taxa_are_rows(ps_agg)) otu_mat <- t(otu_mat)

  # Select top taxa by mean abundance
  ps_rel   <- transform_sample_counts(ps_agg, function(x) x / sum(x))
  rel_mat  <- as.matrix(otu_table(ps_rel))
  if (!taxa_are_rows(ps_rel)) rel_mat <- t(rel_mat)
  top_taxa <- names(sort(rowMeans(rel_mat), decreasing = TRUE))[
    seq_len(min(top_n_taxa, nrow(rel_mat)))
  ]
  otu_sub  <- otu_mat[top_taxa, ]

  # Transform microbiome
  if (transform == "clr") {
    mb_mat <- t(apply(otu_sub + 0.5, 2, function(x) log(x) - mean(log(x))))
  } else {
    mb_mat <- t(apply(otu_sub, 2, function(x) x / sum(x)))
  }

  # Align samples
  shared      <- intersect(rownames(mb_mat), colnames(external_mat))
  if (length(shared) == 0) stop("No shared samples between phyloseq and external matrix.")
  mb_sub      <- mb_mat[shared, , drop = FALSE]
  ext_sub     <- t(external_mat[, shared, drop = FALSE])

  # Select top external features by variance
  if (ncol(ext_sub) > top_n_ext) {
    ext_var <- apply(ext_sub, 2, var, na.rm = TRUE)
    ext_sub <- ext_sub[, order(ext_var, decreasing = TRUE)[seq_len(top_n_ext)]]
  }

  cat("  Microbiome features:", ncol(mb_sub), "\n")
  cat("  External features:", ncol(ext_sub), "\n")
  cat("  Shared samples:", length(shared), "\n")

  # Pairwise correlations
  result_list <- lapply(colnames(mb_sub), function(tx) {
    lapply(colnames(ext_sub), function(ef) {
      x_vals <- mb_sub[, tx]
      y_vals <- ext_sub[, ef]
      complete_idx <- complete.cases(x_vals, y_vals)
      if (sum(complete_idx) < 5) return(NULL)

      ct <- tryCatch(
        cor.test(x_vals[complete_idx], y_vals[complete_idx],
                 method = method, exact = FALSE),
        error = function(e) NULL
      )
      if (is.null(ct)) return(NULL)

      data.frame(
        taxon    = tx,
        external = ef,
        rho      = round(ct$estimate, 4),
        p_value  = ct$p.value,
        n        = sum(complete_idx),
        stringsAsFactors = FALSE
      )
    }) %>% bind_rows()
  }) %>% bind_rows()

  results_df <- result_list %>%
    mutate(
      q_value     = p.adjust(p_value, method = p_adjust),
      significant = q_value < alpha,
      sig_label   = case_when(
        q_value < 0.001 ~ "***",
        q_value < 0.01  ~ "**",
        q_value < 0.05  ~ "*",
        TRUE             ~ ""
      )
    ) %>%
    arrange(q_value)

  n_sig <- sum(results_df$significant)
  cat("  Significant cross-domain associations:", n_sig, "\n\n")

  # Heatmap
  rho_wide <- results_df %>%
    select(taxon, external, rho) %>%
    pivot_wider(names_from = external, values_from = rho, values_fill = 0)

  rho_mat <- as.matrix(rho_wide[, -1])
  rownames(rho_mat) <- rho_wide$taxon

  sig_wide <- results_df %>%
    select(taxon, external, sig_label) %>%
    pivot_wider(names_from = external, values_from = sig_label, values_fill = "")
  sig_mat <- as.matrix(sig_wide[match(rownames(rho_mat), sig_wide$taxon), -1])

  # Cluster
  if (nrow(rho_mat) > 2 && ncol(rho_mat) > 2) {
    row_ord <- hclust(dist(rho_mat))$order
    col_ord <- hclust(dist(t(rho_mat)))$order
    rho_mat <- rho_mat[row_ord, col_ord]
    sig_mat <- sig_mat[row_ord, col_ord]
  }

  heat_df <- as.data.frame(rho_mat) %>%
    rownames_to_column("taxon") %>%
    pivot_longer(-taxon, names_to = "external", values_to = "rho") %>%
    mutate(
      taxon    = factor(taxon, levels = rev(rownames(rho_mat))),
      external = factor(external, levels = colnames(rho_mat))
    )

  sig_df <- as.data.frame(sig_mat) %>%
    rownames_to_column("taxon") %>%
    pivot_longer(-taxon, names_to = "external", values_to = "sig")

  heat_df <- left_join(heat_df, sig_df, by = c("taxon", "external"))

  p <- ggplot(heat_df, aes(x = external, y = taxon, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sig), size = 3, vjust = 0.75,
              colour = "black", fontface = "bold") +
    scale_fill_gradient2(
      low = "#3498db", mid = "white", high = "#e74c3c",
      midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      title    = "Cross-domain correlation heatmap",
      subtitle = paste0("Microbiome (", ncol(mb_sub), " taxa) × External (",
                        ncol(ext_sub), " features) | ",
                        n_sig, " significant pairs (q < ", alpha, ")"),
      x = "External features", y = NULL
    ) +
    theme_microbiome() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y  = element_text(face = "italic", size = 7),
      panel.border = element_rect(colour = "grey80", fill = NA)
    )

  return(list(results = results_df, plot = p, n_sig = n_sig))
}


# =============================================================================
# SECTION 8 — MIXED-EFFECTS ASSOCIATION MODEL
# =============================================================================

#' Test taxa-metadata associations using linear mixed-effects models.
#'
#' Accounts for repeated measures or clustered data (e.g. multiple
#' samples per subject) by including a random intercept for subject.
#'
#' @param ps            A phyloseq object.
#' @param fixed_var     Fixed-effect metadata variable (the predictor).
#' @param random_var    Random-effect variable (e.g. "subject_id").
#' @param rank          Taxonomic rank. Default = "Genus".
#' @param covariates    Additional fixed-effect covariates.
#' @param transform     Abundance transformation. Default = "clr".
#' @param alpha         Significance threshold. Default = 0.05.
#' @return A tidy results data frame.

test_mixed_effects <- function(ps,
                                fixed_var   = "timepoint",
                                random_var  = "subject_id",
                                rank        = "Genus",
                                covariates  = NULL,
                                transform   = "clr",
                                alpha       = 0.05) {

  cat("=== Mixed-effects association model ===\n")

  if (!pkg_available("lme4") || !pkg_available("lmerTest")) {
    cat("  lme4 and lmerTest required.\n")
    cat("  Install: install.packages(c('lme4', 'lmerTest'))\n\n")
    return(NULL)
  }

  library(lme4)
  library(lmerTest)

  ps_agg   <- tax_glom(ps, taxrank = rank, NArm = FALSE)
  new_names <- as.character(tax_table(ps_agg)[, rank])
  new_names[is.na(new_names)] <- paste0("Unknown_", seq_len(sum(is.na(new_names))))
  taxa_names(ps_agg) <- new_names

  otu_mat  <- as.matrix(otu_table(ps_agg))
  if (!taxa_are_rows(ps_agg)) otu_mat <- t(otu_mat)

  # Prevalence filter
  prev     <- rowSums(otu_mat > 0) / ncol(otu_mat)
  otu_mat  <- otu_mat[prev >= 0.20, ]

  if (transform == "clr") {
    feat_mat <- t(apply(otu_mat + 0.5, 2, function(x) log(x) - mean(log(x))))
  } else {
    feat_mat <- t(apply(otu_mat, 2, function(x) x / sum(x)))
  }

  meta_df  <- data.frame(sample_data(ps_agg)) %>%
    rownames_to_column("sample")

  shared   <- intersect(rownames(feat_mat), meta_df$sample)
  feat_sub <- feat_mat[shared, ]
  meta_sub <- meta_df[match(shared, meta_df$sample), ]

  # Build formula
  cov_str  <- if (!is.null(covariates)) paste("+", paste(covariates, collapse = " + ")) else ""
  formula_str <- paste0("abundance ~ ", fixed_var, cov_str,
                         " + (1|", random_var, ")")

  cat("  Formula:", formula_str, "\n")
  cat("  Testing", ncol(feat_sub), "taxa...\n")

  results <- lapply(colnames(feat_sub), function(tx) {
    model_df <- data.frame(
      abundance = feat_sub[, tx],
      meta_sub,
      stringsAsFactors = FALSE
    )

    tryCatch({
      fit  <- lmer(as.formula(formula_str), data = model_df, REML = FALSE)
      coef_tbl <- as.data.frame(summary(fit)$coefficients)

      # Extract row for fixed_var
      row_idx  <- grep(fixed_var, rownames(coef_tbl), value = TRUE)
      if (length(row_idx) == 0) return(NULL)

      lapply(row_idx, function(rn) {
        data.frame(
          taxon        = tx,
          predictor    = rn,
          estimate     = round(coef_tbl[rn, "Estimate"], 4),
          se           = round(coef_tbl[rn, "Std. Error"], 4),
          t_value      = round(coef_tbl[rn, "t value"], 4),
          p_value      = coef_tbl[rn, "Pr(>|t|)"],
          stringsAsFactors = FALSE
        )
      }) %>% bind_rows()
    }, error = function(e) NULL)
  }) %>% bind_rows()

  results <- results %>%
    mutate(
      q_value     = p.adjust(p_value, method = "BH"),
      significant = q_value < alpha,
      sig_label   = case_when(
        q_value < 0.001 ~ "***",
        q_value < 0.01  ~ "**",
        q_value < 0.05  ~ "*",
        TRUE             ~ ""
      )
    ) %>%
    arrange(q_value)

  n_sig <- sum(results$significant)
  cat("  Significant taxa:", n_sig, "\n\n")

  return(results)
}


# =============================================================================
# SECTION 9 — COMPLETE CORRELATION WORKFLOW WRAPPER
# =============================================================================

#' Run the complete correlation and association analysis pipeline.
#'
#' @param ps          A filtered phyloseq object.
#' @param group_var   Grouping variable for stratified scatter plots.
#' @param meta_vars   Numeric metadata variables to correlate against taxa.
#' @param rank        Taxonomic rank. Default = "Genus".
#' @param alpha       Significance threshold. Default = 0.05.
#' @param output_dir  Directory for outputs.
#' @return A named list of all results.

run_correlation_analysis <- function(ps,
                                      group_var  = NULL,
                                      meta_vars  = NULL,
                                      rank       = "Genus",
                                      alpha      = 0.05,
                                      output_dir = "correlation_output") {

  cat("\n", strrep("=", 60), "\n")
  cat("  CORRELATION AND ASSOCIATION PIPELINE\n")
  cat(strrep("=", 60), "\n\n")

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  results <- list()

  # --- Taxa-metadata correlations -------------------------------------------
  if (!is.null(meta_vars) || any(sapply(data.frame(sample_data(ps)), is.numeric))) {
    cat("--- Analysis 1: Taxa-metadata correlations ---\n")
    results$cor_taxa_meta <- correlate_taxa_metadata(
      ps, rank = rank, meta_vars = meta_vars,
      transform = "clr", alpha = alpha
    )

    write.csv(results$cor_taxa_meta,
              file.path(output_dir, "taxa_metadata_correlations.csv"),
              row.names = FALSE)

    cat("--- Plot 1: Correlation heatmap ---\n")
    results$p_heatmap <- plot_taxa_metadata_heatmap(
      results$cor_taxa_meta, top_n = 30, alpha = alpha
    )
    ggsave(file.path(output_dir, "01_taxa_metadata_heatmap.pdf"),
           results$p_heatmap, width = 12, height = 14)

    cat("--- Plot 2: Bubble plot ---\n")
    results$p_bubbles <- plot_association_bubbles(
      results$cor_taxa_meta, top_n = 60, alpha = alpha
    )
    ggsave(file.path(output_dir, "02_association_bubbles.pdf"),
           results$p_bubbles, width = 12, height = 14)

    cat("--- Plot 3: Top scatter plots ---\n")
    results$p_scatter <- plot_top_associations(
      ps, results$cor_taxa_meta,
      rank = rank, top_n = 12, group_var = group_var
    )
    if (!is.null(results$p_scatter)) {
      ggsave(file.path(output_dir, "03_top_association_scatters.pdf"),
             results$p_scatter, width = 16, height = 12)
    }
  }

  # --- Taxa-taxa correlations -----------------------------------------------
  cat("--- Plot 4: Taxa-taxa correlation matrix ---\n")
  tt_res <- compute_taxa_taxa_correlation(
    ps, rank = rank, top_n = 40,
    transform = "clr", alpha = alpha
  )
  results$p_taxa_taxa <- tt_res$plot
  results$taxa_taxa_cor <- tt_res
  ggsave(file.path(output_dir, "04_taxa_taxa_correlation.pdf"),
         tt_res$plot, width = 14, height = 13)

  # --- Metadata-metadata correlations ----------------------------------------
  cat("--- Plot 5: Metadata correlations ---\n")
  mm_res <- plot_metadata_correlations(ps, meta_vars = meta_vars, alpha = alpha)
  if (!is.null(mm_res)) {
    results$p_meta_meta <- mm_res$plot
    ggsave(file.path(output_dir, "05_metadata_correlations.pdf"),
           mm_res$plot, width = 10, height = 9)
  }

  cat("\n", strrep("=", 60), "\n")
  cat("  CORRELATION PIPELINE COMPLETE\n")
  cat(strrep("=", 60), "\n")
  cat("  Output saved to:", output_dir, "\n")
  cat("  Plots:  up to 5 PDF files\n")
  cat("  Tables: taxa-metadata correlation results\n\n")

  return(invisible(results))
}


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# ps <- readRDS("qc_output/phyloseq_qc_filtered.rds")

# --- Option A: Full pipeline -------------------------------------------------
# results <- run_correlation_analysis(
#   ps         = ps,
#   group_var  = "disease_status",
#   meta_vars  = c("age", "bmi", "calprotectin", "crp", "shannon"),
#   rank       = "Genus",
#   alpha      = 0.05,
#   output_dir = "results/correlation"
# )

# --- Option B: Individual steps ----------------------------------------------
# cor_results <- correlate_taxa_metadata(ps,
#                 meta_vars = c("age", "bmi", "calprotectin"),
#                 rank = "Genus", transform = "clr")
#
# cor_results %>% filter(significant) %>%
#   select(taxon, metadata, rho, q_value) %>% head(20)
#
# plot_taxa_metadata_heatmap(cor_results, top_n = 30)
# plot_association_bubbles(cor_results, top_n = 50)
# plot_top_associations(ps, cor_results, top_n = 12, group_var = "disease_status")

# --- Option C: Taxa-taxa only ------------------------------------------------
# tt_res <- compute_taxa_taxa_correlation(ps, rank = "Genus", top_n = 50)
# tt_res$plot
# sum(tt_res$p_matrix[upper.tri(tt_res$p_matrix)] < 0.05)

# --- Option D: Cross-domain (microbiome x metabolomics) ----------------------
# metabolite_mat <- read.csv("metabolomics.csv", row.names = 1)
# cross_res <- correlate_cross_domain(ps, external_mat = metabolite_mat,
#                                     rank = "Genus", alpha = 0.05)
# cross_res$plot
# cross_res$results %>% filter(significant) %>% head(20)

# --- Option E: Mixed-effects (longitudinal data) ----------------------------
# lme_res <- test_mixed_effects(ps,
#              fixed_var  = "timepoint",
#              random_var = "subject_id",
#              covariates = c("age", "sex"))
# lme_res %>% filter(significant) %>% head(10)
