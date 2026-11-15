# =============================================================================
# MICROBIOME TAXONOMIC COMPOSITION ANALYSIS
# =============================================================================
# Description : Comprehensive taxonomic profiling and visualisation pipeline
# Input       : Filtered phyloseq object (output from microbiome_qc.R)
# Output      : Composition plots, abundance tables, diversity summaries
# Author      : Patricia
# Dependencies: phyloseq, ggplot2, dplyr, tidyr, scales, patchwork,
#               RColorBrewer, microbiome, vegan, stringr, forcats
# =============================================================================

# --- 1. LOAD LIBRARIES -------------------------------------------------------

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
  library(RColorBrewer)
  library(microbiome)
  library(vegan)
  library(stringr)
  library(forcats)
  library(tibble)
})


# Shared helpers (fast_tax_glom, group_palette) live in R/utils.R.
# Load them when this file is sourced on its own, e.g. by the test suite.
if (!exists("fast_tax_glom") && file.exists("R/utils.R")) source("R/utils.R")

# Global theme (consistent with QC module)
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

# Colour palette for taxa — extended to 20 colours
TAXA_PALETTE <- c(
  "#3498db", "#e74c3c", "#2ecc71", "#f39c12", "#9b59b6",
  "#1abc9c", "#e67e22", "#34495e", "#e91e63", "#00bcd4",
  "#8bc34a", "#ff5722", "#607d8b", "#795548", "#ffc107",
  "#673ab7", "#009688", "#ff9800", "#4caf50", "#f44336",
  "#b0bec5"   # "Other" always last
)


# =============================================================================
# SECTION 1 — AGGLOMERATE AND TRANSFORM
# =============================================================================

#' Agglomerate taxa to a given taxonomic rank and apply abundance transformation.
#'
#' @param ps        A phyloseq object.
#' @param rank      Taxonomic rank. One of "Phylum","Class","Order","Family","Genus","Species".
#' @param transform Abundance transformation: "relative", "clr", "log10", or "none".
#' @param top_n     Keep only the top N taxa by mean abundance. All others merged as "Other".
#' @return A phyloseq object agglomerated and transformed at the chosen rank.

agglomerate_taxa <- function(ps,
                              rank      = "Genus",
                              transform = "relative",
                              top_n     = 20) {

  cat("=== Agglomerating to", rank, "level ===\n")

  # Validate rank
  available_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  if (!rank %in% available_ranks) {
    stop("rank must be one of: ", paste(available_ranks, collapse = ", "))
  }

  # Agglomerate
  ps_agg <- fast_tax_glom(ps, taxrank = rank, NArm = FALSE)

  # Rename taxa using the chosen rank label
  new_names <- as.character(tax_table(ps_agg)[, rank])
  new_names[is.na(new_names) | new_names == "NA"] <- paste0("Unknown_", rank, "_", seq_len(sum(is.na(new_names) | new_names == "NA")))
  new_names <- make.unique(new_names, sep = "_dup")
  taxa_names(ps_agg) <- new_names

  # Handle NA taxa names




  cat("  Taxa after agglomeration:", ntaxa(ps_agg), "\n")

  # --- Top N filtering -------------------------------------------------------
  if (!is.null(top_n) && ntaxa(ps_agg) > top_n) {
    otu_mat   <- as.matrix(otu_table(ps_agg))
    if (!taxa_are_rows(ps_agg)) otu_mat <- t(otu_mat)

    mean_abund <- rowMeans(otu_mat)
    top_taxa   <- names(sort(mean_abund, decreasing = TRUE))[seq_len(top_n)]
    other_taxa <- setdiff(taxa_names(ps_agg), top_taxa)

    # Sum "Other" taxa into a single row
    if (length(other_taxa) > 0) {
      other_counts <- colSums(otu_mat[other_taxa, , drop = FALSE])
      top_otu      <- otu_mat[top_taxa, , drop = FALSE]
      other_row    <- matrix(other_counts, nrow = 1,
                             dimnames = list("Other", colnames(top_otu)))
      combined_otu <- rbind(top_otu, other_row)

      # Rebuild taxonomy for "Other"
      tax_top         <- as.data.frame(tax_table(ps_agg))[top_taxa, , drop = FALSE]
      other_tax_row   <- as.data.frame(matrix("Other", nrow = 1,
                                              ncol = ncol(tax_top),
                                              dimnames = list("Other", colnames(tax_top))))
      combined_tax    <- rbind(tax_top, other_tax_row)

      ps_agg <- phyloseq(
        otu_table(combined_otu, taxa_are_rows = TRUE),
        tax_table(as.matrix(combined_tax)),
        sample_data(ps_agg)
      )
    } else {
      ps_agg <- prune_taxa(top_taxa, ps_agg)
    }

    cat("  Kept top", top_n, "taxa + Other\n")
  }

  # --- Transformation -------------------------------------------------------
  if (transform == "relative") {
    ps_agg <- transform_sample_counts(ps_agg, function(x) x / sum(x))
    cat("  Applied relative abundance transformation\n")

  } else if (transform == "clr") {
    # Centre log-ratio: requires pseudocount for zeros
    otu_clr <- microbiome::transform(ps_agg, transform = "clr")
    ps_agg  <- otu_clr
    cat("  Applied CLR transformation\n")

  } else if (transform == "log10") {
    ps_agg <- transform_sample_counts(ps_agg, function(x) log10(x + 1))
    cat("  Applied log10(x+1) transformation\n")
  }

  cat("\n")
  return(ps_agg)
}


# =============================================================================
# SECTION 2 — STACKED BAR PLOTS
# =============================================================================

#' Plot stacked bar charts of taxonomic composition per sample or group.
#'
#' @param ps          A phyloseq object (ideally relative-abundance transformed).
#' @param rank        Taxonomic rank to display.
#' @param group_var   Optional metadata variable for grouping samples.
#' @param facet_var   Optional metadata variable for faceting.
#' @param sort_by     How to order samples: "group", "dominant_taxon", or "none".
#' @param show_legend Whether to show the legend. Default = TRUE.
#' @return A ggplot object.

plot_composition_bars <- function(ps,
                                   rank       = "Phylum",
                                   group_var  = NULL,
                                   facet_var  = NULL,
                                   sort_by    = "group",
                                   show_legend = TRUE) {

  cat("=== Stacked bar composition plot ===\n")

  # Melt phyloseq to long format
  ps_melt <- psmelt(ps) %>%
    mutate(taxon = as.character(.data[[rank]])) %>%
    mutate(taxon = as.character(taxon))

  # Handle NA
  ps_melt$taxon[is.na(ps_melt$taxon)] <- "Unknown"

  # Assign colours — "Other" always gets grey (#b0bec5), regardless of top N
  all_taxa     <- unique(ps_melt$taxon)
  non_other    <- all_taxa[all_taxa != "Other"]
  has_other    <- "Other" %in% all_taxa
  taxa_ordered <- c(non_other, if (has_other) "Other")

  n_non_other       <- length(non_other)
  non_other_palette <- TAXA_PALETTE[TAXA_PALETTE != "#b0bec5"]
  # Recycle palette if more taxa than colours
  recycled_palette  <- rep(non_other_palette,
                           length.out = n_non_other)
  col_vals <- setNames(recycled_palette, non_other)
  if (has_other) col_vals["Other"] <- "#b0bec5"

  # --- Sample ordering -------------------------------------------------------
  if (sort_by == "group" && !is.null(group_var)) {
    if (!is.null(facet_var) && facet_var %in% colnames(ps_melt)) {
      # Sort within each facet panel independently
      sample_order <- ps_melt %>%
        select(Sample, all_of(c(group_var, facet_var))) %>%
        distinct() %>%
        arrange(.data[[facet_var]], .data[[group_var]], Sample) %>%
        pull(Sample)
    } else {
      sample_order <- ps_melt %>%
        arrange(.data[[group_var]], Sample) %>%
        pull(Sample) %>% unique()
    }
    ps_melt <- ps_melt %>%
      mutate(Sample = factor(Sample, levels = sample_order))
  } else if (sort_by == "dominant_taxon") {
    dominant <- ps_melt %>%
      group_by(Sample, taxon) %>%
      summarise(abund = sum(Abundance), .groups = "drop") %>%
      group_by(Sample) %>%
      slice_max(abund, n = 1) %>%
      arrange(taxon, desc(abund))
    ps_melt <- ps_melt %>%
      mutate(Sample = factor(Sample, levels = dominant$Sample))
  }

  # Ensure taxa factor order (Other always last)
  ps_melt$taxon <- factor(ps_melt$taxon, levels = rev(taxa_ordered))

  # --- Base plot -------------------------------------------------------------
  p <- ggplot(ps_melt, aes(x = Sample, y = Abundance, fill = taxon)) +
    geom_col(position = "stack", width = 0.9, colour = NA) +
    scale_fill_manual(values = col_vals, name = rank) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
    labs(
      title    = paste0("Taxonomic composition — ", rank, " level"),
      subtitle = paste0(nsamples(ps), " samples"),
      x        = NULL,
      y        = "Relative abundance"
    ) +
    theme_microbiome() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
      panel.border = element_blank(),
      axis.line.y  = element_line(colour = "grey80")
    )

  # Facet: one panel per level, each as wide as its number of samples.
  # The same facet is used for the bars and for the group bar above them,
  # so the two always line up.
  has_facet   <- !is.null(facet_var) && facet_var %in% colnames(ps_melt)
  facet_layer <- if (has_facet) {
    facet_grid(cols     = vars(.data[[facet_var]]),
               scales   = "free_x", space = "free_x",
               labeller = label_wrap_gen(width = 12))
  } else {
    NULL
  }
  p <- p + facet_layer

  # Add group annotation bar if group_var provided
  if (!is.null(group_var) && group_var %in% colnames(ps_melt)) {
    group_df <- ps_melt %>%
      select(Sample, all_of(unique(c(group_var, if (has_facet) facet_var)))) %>%
      distinct()

    p_group <- ggplot(group_df, aes(x = Sample, y = 1,
                                    fill = .data[[group_var]])) +
      geom_col(width = 0.9) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_fill_group(name = group_var) +
      facet_layer +
      theme_void() +
      theme(
        legend.position  = "right",
        axis.text.x      = element_blank(),
        strip.background = element_rect(fill = "#2c3e50", colour = NA),
        strip.text       = element_text(colour = "white", face = "bold",
                                        size = 9,
                                        margin = ggplot2::margin(3, 2, 3, 2))
      )

    # Title and facet labels go at the very top, above the group bar
    p <- p +
      labs(title = NULL, subtitle = NULL) +
      theme(strip.background = element_blank(), strip.text = element_blank())

    p <- p_group / p +
      plot_layout(heights = c(if (has_facet) 0.14 else 0.08, 1), guides = "collect") +
      plot_annotation(
        title    = paste0("Taxonomic composition — ", rank, " level"),
        subtitle = paste0(nsamples(ps), " samples"),
        theme    = theme(
          plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(colour = "grey40", size = 10)
        )
      )
  }

  if (!show_legend) p <- p + theme(legend.position = "none")

  return(p)
}


#' Plot mean composition per group as grouped bars with error bars.
#'
#' @param ps        A phyloseq object (relative abundance).
#' @param rank      Taxonomic rank.
#' @param group_var Metadata variable defining groups.
#' @param top_n     Number of top taxa to display ("Other" is shown in
#'   addition and does not count towards N).
#' @param error_bar What the error bars show: "se" (standard error of the
#'   mean, default), "sd" (standard deviation), "ci95" (95% confidence
#'   interval of the mean) or "none".
#' @return A ggplot object.

plot_mean_composition <- function(ps,
                                   rank      = "Phylum",
                                   group_var = "group",
                                   top_n     = 10,
                                   error_bar = c("se", "sd", "ci95", "none")) {

  cat("=== Mean composition per group ===\n")

  error_bar <- match.arg(error_bar)

  ps_melt <- psmelt(ps) %>%
    mutate(taxon = as.character(.data[[rank]])) %>%
    mutate(taxon = as.character(taxon))

  ps_melt$taxon[is.na(ps_melt$taxon)] <- "Unknown"

  # Select top N taxa overall. "Other" is not a taxon: it is kept, shown
  # last and in grey, and does not count towards the N.
  top_taxa <- ps_melt %>%
    filter(taxon != "Other") %>%
    group_by(taxon) %>%
    summarise(mean_abund = mean(Abundance)) %>%
    arrange(desc(mean_abund)) %>%
    slice_head(n = top_n) %>%
    pull(taxon)

  has_other  <- "Other" %in% ps_melt$taxon
  taxa_shown <- c(top_taxa, if (has_other) "Other")

  ps_melt <- ps_melt %>%
    filter(taxon %in% taxa_shown) %>%
    mutate(taxon = factor(taxon, levels = taxa_shown))

  # Summary stats per group
  summary_df <- ps_melt %>%
    group_by(.data[[group_var]], taxon) %>%
    summarise(
      n_samp     = n(),
      mean_abund = mean(Abundance),
      sd_abund   = sd(Abundance),
      .groups    = "drop"
    ) %>%
    mutate(
      se_abund = sd_abund / sqrt(n_samp),
      err = switch(error_bar,
        se   = se_abund,
        sd   = sd_abund,
        ci95 = suppressWarnings(qt(0.975, df = pmax(n_samp - 1, 1))) * se_abund,
        none = NA_real_
      ),
      err = ifelse(n_samp < 2, NA_real_, err),
      # a proportion cannot be below 0
      ymin = pmax(mean_abund - err, 0),
      ymax = mean_abund + err
    )

  err_label <- switch(error_bar,
    se   = "Error bars: standard error of the mean (SE)",
    sd   = "Error bars: standard deviation (SD)",
    ci95 = "Error bars: 95% confidence interval of the mean",
    none = "No error bars"
  )

  # Group labels with the number of samples
  n_per_group <- ps_melt %>%
    dplyr::distinct(Sample, .data[[group_var]]) %>%
    dplyr::count(.data[[group_var]], name = "n_samp")
  x_labels <- setNames(
    paste0(n_per_group[[group_var]], "\n(n = ", n_per_group$n_samp, ")"),
    n_per_group[[group_var]]
  )

  col_vals <- TAXA_PALETTE[seq_along(top_taxa)]
  names(col_vals) <- top_taxa
  if (has_other) col_vals <- c(col_vals, Other = "#b0bec5")

  p <- ggplot(summary_df,
              aes(x = .data[[group_var]], y = mean_abund, fill = taxon)) +
    geom_col(position = position_dodge(0.85), width = 0.8)

  if (error_bar != "none") {
    p <- p +
      geom_errorbar(
        aes(ymin = ymin, ymax = ymax),
        position = position_dodge(0.85), width = 0.25, linewidth = 0.5,
        na.rm = TRUE
      )
  }

  p <- p +
    scale_fill_manual(values = col_vals, name = rank) +
    scale_x_discrete(labels = x_labels) +
    scale_y_continuous(labels = percent_format(accuracy = 0.1),
                       expand  = expansion(mult = c(0, 0.05))) +
    labs(
      title    = paste0("Mean ", rank, " abundance by group"),
      subtitle = paste0("Top ", length(top_taxa), " taxa | ", err_label),
      x        = group_var,
      y        = "Mean relative abundance"
    ) +
    theme_microbiome()

  return(p)
}


# =============================================================================
# SECTION 3 — HEATMAP
# =============================================================================

#' Hierarchically clustered heatmap of taxon abundance across samples.
#'
#' @param ps          A phyloseq object.
#' @param rank        Taxonomic rank.
#' @param top_n       Number of top taxa to include.
#' @param group_var   Optional metadata variable for the group colour bar
#'   and, with block_by_group, for the blocks.
#' @param transform   Transformation to apply: "clr", "log10", or "relative".
#' @param cluster_samples Whether to cluster samples. Default = TRUE.
#'   Without blocks all samples are clustered together and a column
#'   dendrogram is drawn; with blocks, samples are clustered inside each block.
#' @param cluster_rows Whether to cluster taxa. Default = TRUE.
#' @param block_by_group Keep the samples of each group together, one panel
#'   per group. Default = TRUE when group_var is given.
#' @param scale_rows  Scale each taxon to a z-score across samples.
#'   Default = FALSE.
#' @return A patchwork object (heatmap with optional dendrograms and
#'   group colour bar). Clustering uses Euclidean distance and complete linkage.

plot_abundance_heatmap <- function(ps,
                                    rank            = "Genus",
                                    top_n           = 30,
                                    group_var       = NULL,
                                    transform       = "clr",
                                    cluster_samples = TRUE,
                                    cluster_rows    = TRUE,
                                    block_by_group  = !is.null(group_var),
                                    scale_rows      = FALSE) {

  cat("=== Abundance heatmap ===\n")

  # Agglomerate
  ps_agg  <- agglomerate_taxa(ps, rank = rank, transform = transform,
                               top_n = top_n)

  otu_mat <- as(otu_table(ps_agg), "matrix")
  if (!taxa_are_rows(ps_agg)) otu_mat <- t(otu_mat)

  # Remove "Other" from heatmap
  otu_mat <- otu_mat[rownames(otu_mat) != "Other", , drop = FALSE]

  # Optional: scale each taxon (row z-score) so that colours show how a
  # sample compares with the other samples for that taxon
  if (scale_rows) {
    otu_mat <- t(scale(t(otu_mat)))
    otu_mat[is.na(otu_mat)] <- 0   # taxa with no variation
  }

  # Group of each sample (for the colour bar and the blocks)
  meta_df   <- data.frame(sample_data(ps_agg), check.names = FALSE)
  has_group <- !is.null(group_var) && group_var %in% colnames(meta_df)
  block     <- isTRUE(block_by_group) && has_group

  groups     <- NULL
  grp_levels <- NULL
  if (has_group) {
    grp_col    <- meta_df[colnames(otu_mat), group_var]
    grp_levels <- if (is.factor(grp_col)) levels(droplevels(grp_col))
                  else sort(unique(as.character(grp_col)))
    groups     <- as.character(grp_col)
    if (anyNA(groups)) {
      groups[is.na(groups)] <- "NA"
      grp_levels <- c(grp_levels, "NA")
    }
  }

  # Row (taxa) clustering: Euclidean distance, complete linkage
  row_hclust <- NULL
  if (cluster_rows && nrow(otu_mat) > 2) {
    row_hclust <- hclust(dist(otu_mat))
    otu_mat    <- otu_mat[row_hclust$order, , drop = FALSE]
  }

  # Column (sample) ordering
  #   block           : samples stay inside their group; clustered within it
  #   cluster_samples : all samples clustered together (column dendrogram)
  #   neither         : order of the input data
  col_hclust <- NULL
  if (block) {
    col_order <- unlist(lapply(grp_levels, function(g) {
      samps <- colnames(otu_mat)[groups == g]
      if (!cluster_samples || length(samps) <= 2) return(samps)
      samps[hclust(dist(t(otu_mat[, samps, drop = FALSE])))$order]
    }))
    groups  <- groups[match(col_order, colnames(otu_mat))]
    otu_mat <- otu_mat[, col_order, drop = FALSE]
  } else if (cluster_samples && ncol(otu_mat) > 2) {
    col_hclust <- hclust(dist(t(otu_mat)))
    if (has_group) groups <- groups[col_hclust$order]
    otu_mat    <- otu_mat[, col_hclust$order, drop = FALSE]
  }

  sample_levels <- colnames(otu_mat)
  n_taxa        <- nrow(otu_mat)
  n_samp        <- ncol(otu_mat)

  # Melt to long format
  heat_df <- as.data.frame(otu_mat) %>%
    rownames_to_column("taxon") %>%
    pivot_longer(-taxon, names_to = "sample", values_to = "value") %>%
    mutate(
      taxon  = factor(taxon,  levels = rownames(otu_mat)),
      sample = factor(sample, levels = sample_levels)
    )

  if (has_group) {
    grp_of <- setNames(groups, sample_levels)
    heat_df$grp_fill <- factor(grp_of[as.character(heat_df$sample)],
                               levels = grp_levels)
  }

  # One panel per group, each as wide as its number of samples
  facet_layer <- if (block) {
    facet_grid(cols     = vars(grp_fill),
               scales   = "free_x", space = "free_x",
               labeller = label_wrap_gen(width = 12))
  } else {
    NULL
  }

  value_name <- if (scale_rows) "Z-score"
                else ifelse(transform == "clr", "CLR",
                            ifelse(transform == "log10", "Log10", "Rel. abund."))
  midpoint   <- if (scale_rows) 0 else median(heat_df$value)

  p <- ggplot(heat_df, aes(x = sample, y = taxon, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.15) +
    scale_fill_gradient2(
      low      = "#3498db",
      mid      = "white",
      high     = "#e74c3c",
      midpoint = midpoint,
      name     = value_name
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    facet_layer +
    labs(x = "Samples", y = rank) +
    theme_microbiome() +
    theme(
      axis.text.x      = if (n_samp <= 60) element_text(angle = 45, hjust = 1, size = 7)
                         else element_blank(),
      axis.ticks.x     = if (n_samp <= 60) element_line() else element_blank(),
      axis.text.y      = element_text(size = 8, face = "italic"),
      panel.border     = element_rect(colour = "grey80", fill = NA),
      panel.spacing.x  = unit(4, "pt"),
      strip.background = element_blank(),
      strip.text       = element_blank()
    )

  has_dendro_pkg <- requireNamespace("ggdendro", quietly = TRUE)

  # Group colour bar (carries the block labels when blocking)
  p_annot <- NULL
  if (has_group) {
    annot_df <- data.frame(
      sample   = factor(sample_levels, levels = sample_levels),
      grp_fill = factor(groups, levels = grp_levels)
    )

    p_annot <- ggplot(annot_df, aes(x = sample, y = 1, fill = grp_fill)) +
      geom_tile() +
      scale_fill_group(name = group_var) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      facet_layer +
      theme_void() +
      theme(
        legend.position  = "right",
        plot.margin      = ggplot2::margin(0, 0, 0, 0),
        panel.spacing.x  = unit(4, "pt"),
        strip.background = element_rect(fill = "#2c3e50", colour = NA),
        strip.text       = element_text(colour = "white", face = "bold", size = 8,
                                        margin = ggplot2::margin(3, 2, 3, 2))
      )
  }

  # Row dendrogram (left). Leaf i of the dendrogram is row i of the heatmap.
  p_dendro <- NULL
  if (!is.null(row_hclust) && has_dendro_pkg) {
    dendro_data <- ggdendro::dendro_data(as.dendrogram(row_hclust),
                                         type = "rectangle")

    p_dendro <- ggplot(ggdendro::segment(dendro_data)) +
      geom_segment(aes(x = y, y = x, xend = yend, yend = xend),
                   linewidth = 0.4, colour = "#2c3e50") +
      scale_x_reverse(expand = c(0.02, 0)) +
      scale_y_continuous(limits = c(0.5, n_taxa + 0.5), expand = c(0, 0)) +
      theme_void() +
      theme(plot.margin = ggplot2::margin(0, 0, 0, 0))

    # Move taxa labels to right side when dendrogram is present
    p <- p +
      scale_y_discrete(position = "right", expand = c(0, 0)) +
      theme(
        axis.text.y  = element_text(size = 8, face = "italic",
                                     hjust = 0, margin = ggplot2::margin(l = 4)),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()
      )
  }

  # Column dendrogram (top), only when all samples are clustered together
  p_col_dendro <- NULL
  if (!is.null(col_hclust) && has_dendro_pkg) {
    col_dd       <- ggdendro::dendro_data(as.dendrogram(col_hclust), type = "rectangle")
    p_col_dendro <- ggplot(ggdendro::segment(col_dd)) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend),
                   linewidth = 0.4, colour = "#2c3e50") +
      scale_x_continuous(limits = c(0.5, n_samp + 0.5), expand = c(0, 0)) +
      scale_y_continuous(expand = c(0.02, 0)) +
      theme_void() +
      theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
  }

  # ── Assemble final layout ─────────────────────────────────────────────────
  #   [spacer] column dendrogram
  #   [spacer] group colour bar
  #   row dendrogram | heatmap
  has_left <- !is.null(p_dendro)
  tops     <- Filter(Negate(is.null), list(p_col_dendro, p_annot))
  plots    <- list()
  for (tp in tops) {
    if (has_left) plots <- c(plots, list(patchwork::plot_spacer()))
    plots <- c(plots, list(tp))
  }
  if (has_left) plots <- c(plots, list(p_dendro))
  plots <- c(plots, list(p))

  heights <- c(if (!is.null(p_col_dendro)) 0.15,
               if (!is.null(p_annot)) 0.05,
               1)

  col_label <- if (block) {
    paste0("blocked by ", group_var,
           if (cluster_samples) ", clustered within each block" else "")
  } else if (!is.null(col_hclust)) {
    "all clustered together"
  } else {
    "in the order of the data"
  }

  out <- patchwork::wrap_plots(plots,
                               ncol    = if (has_left) 2 else 1,
                               widths  = if (has_left) c(0.15, 1) else 1,
                               heights = heights,
                               guides  = "collect") +
    patchwork::plot_annotation(
      title    = paste0(rank, "-level abundance heatmap"),
      subtitle = paste0("Top ", n_taxa, " taxa | Values: ", value_name,
                        if (scale_rows) paste0(" of ", transform) else "",
                        " | Samples: ", col_label),
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey40", size = 10)
      )
    )

  return(out)
}


# =============================================================================
# SECTION 4 — TAXONOMIC SUMMARY TABLES
# =============================================================================

#' Generate a tidy abundance summary table at a chosen rank.
#'
#' @param ps          A phyloseq object.
#' @param rank        Taxonomic rank.
#' @param group_var   Optional metadata variable for group-wise summaries.
#' @return A data frame with mean, SD, median, and prevalence per taxon.

make_abundance_table <- function(ps,
                                  rank      = "Genus",
                                  group_var = NULL) {

  cat("=== Generating abundance summary table ===\n")

  ps_rel  <- transform_sample_counts(ps, function(x) x / sum(x))
  ps_agg  <- fast_tax_glom(ps_rel, taxrank = rank, NArm = FALSE)

  ps_melt <- psmelt(ps_agg) %>%
    mutate(taxon = as.character(.data[[rank]])) %>%
    mutate(taxon = ifelse(is.na(taxon), "Unknown", as.character(taxon)))

  # Add taxonomy string
  tax_string_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  available_cols  <- intersect(tax_string_cols, colnames(ps_melt))
  rank_idx        <- which(available_cols == rank)
  lineage_cols    <- available_cols[seq_len(rank_idx)]

  ps_melt <- ps_melt %>%
    mutate(lineage = apply(select(., all_of(lineage_cols)), 1,
                           function(x) paste(na.omit(x), collapse = "; ")))

  # Overall summary
  overall_summary <- ps_melt %>%
    group_by(taxon, lineage) %>%
    summarise(
      mean_relabund  = mean(Abundance),
      sd_relabund    = sd(Abundance),
      median_relabund = median(Abundance),
      prevalence     = mean(Abundance > 0),
      n_samples      = sum(Abundance > 0),
      .groups        = "drop"
    ) %>%
    arrange(desc(mean_relabund)) %>%
    mutate(
      mean_pct   = round(mean_relabund * 100, 3),
      sd_pct     = round(sd_relabund * 100, 3),
      median_pct = round(median_relabund * 100, 3),
      prev_pct   = round(prevalence * 100, 1)
    ) %>%
    select(taxon, lineage, mean_pct, sd_pct, median_pct, prev_pct, n_samples)

  cat("  Unique taxa:", nrow(overall_summary), "\n")

  # Group-wise summary if requested
  if (!is.null(group_var) && group_var %in% colnames(ps_melt)) {
    group_summary <- ps_melt %>%
      group_by(.data[[group_var]], taxon) %>%
      summarise(
        mean_pct = round(mean(Abundance) * 100, 3),
        sd_pct   = round(sd(Abundance) * 100, 3),
        prev_pct = round(mean(Abundance > 0) * 100, 1),
        .groups  = "drop"
      )

    cat("  Group-wise summary produced for:", group_var, "\n\n")
    return(list(overall = overall_summary, by_group = group_summary))
  }

  cat("\n")
  return(list(overall = overall_summary))
}


# =============================================================================
# SECTION 5 — PHYLUM-TO-GENUS DRILL-DOWN
# =============================================================================

#' For a selected phylum, show the genus-level breakdown.
#'
#' @param ps          A phyloseq object.
#' @param phylum      Name of phylum to drill into.
#' @param group_var   Optional metadata variable for grouping.
#' @param top_n       Number of top genera to show. Default = 15.
#' @return A ggplot object.

plot_phylum_drilldown <- function(ps,
                                   phylum    = "Firmicutes",
                                   group_var = NULL,
                                   top_n     = 15) {

  cat("=== Phylum drill-down:", phylum, "===\n")

  # Filter to selected phylum
  ps_rel  <- transform_sample_counts(ps, function(x) x / sum(x))
  ps_phy  <- subset_taxa(ps_rel, Phylum == phylum)

  if (ntaxa(ps_phy) == 0) {
    stop("No taxa found for phylum: ", phylum,
         ". Check spelling and available phyla.")
  }

  cat("  Taxa in", phylum, ":", ntaxa(ps_phy), "\n")

  # Agglomerate to genus
  ps_gen <- fast_tax_glom(ps_phy, taxrank = "Genus", NArm = FALSE)
  genera <- as.character(tax_table(ps_gen)[, "Genus"])
  genera[is.na(genera)] <- "Unknown genus"
  taxa_names(ps_gen) <- make.unique(genera, sep = "_")

  # Keep top N genera
  mean_abund <- rowMeans(as.matrix(otu_table(ps_gen)))
  top_genera <- names(sort(mean_abund, decreasing = TRUE))[seq_len(min(top_n, length(mean_abund)))]
  ps_gen     <- prune_taxa(top_genera, ps_gen)

  ps_melt <- psmelt(ps_gen) %>%
    mutate(Genus = factor(OTU, levels = rev(top_genera)))

  col_vals <- TAXA_PALETTE[seq_len(length(top_genera))]
  names(col_vals) <- top_genera

  if (!is.null(group_var) && group_var %in% colnames(ps_melt)) {
    # Boxplot per genus per group
    p <- ggplot(ps_melt,
                aes(x = .data[[group_var]], y = Abundance,
                    fill = .data[[group_var]])) +
      geom_boxplot(outlier.size = 0.8, alpha = 0.75) +
      geom_jitter(width = 0.15, alpha = 0.4, size = 0.8) +
      facet_wrap(~ Genus, scales = "free_y", ncol = 5) +
      scale_fill_group() +
      scale_y_continuous(labels = percent_format(accuracy = 0.01)) +
      labs(
        title    = paste0("Genus-level breakdown: ", phylum),
        subtitle = paste0("Top ", length(top_genera), " genera"),
        x        = group_var,
        y        = "Relative abundance",
        fill     = group_var
      ) +
      theme_microbiome() +
      theme(
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none"
      )
  } else {
    # Bar plot ordered by abundance
    genus_summary <- ps_melt %>%
      group_by(Genus) %>%
      summarise(
        mean_abund = mean(Abundance),
        se_abund   = sd(Abundance) / sqrt(n()),
        .groups    = "drop"
      ) %>%
      mutate(Genus = fct_reorder(Genus, mean_abund))

    p <- ggplot(genus_summary, aes(x = mean_abund, y = Genus, fill = Genus)) +
      geom_col(width = 0.75, show.legend = FALSE) +
      geom_errorbarh(
        aes(xmin = mean_abund - se_abund, xmax = mean_abund + se_abund),
        height = 0.3, linewidth = 0.5
      ) +
      scale_fill_manual(values = col_vals) +
      scale_x_continuous(labels = percent_format(accuracy = 0.01),
                         expand  = expansion(mult = c(0, 0.1))) +
      labs(
        title    = paste0("Genus-level breakdown: ", phylum),
        subtitle = paste0("Top ", nrow(genus_summary), " genera | Error bars = ±1 SE"),
        x        = "Mean relative abundance",
        y        = NULL
      ) +
      theme_microbiome() +
      theme(axis.text.y = element_text(face = "italic"))
  }

  return(p)
}


# =============================================================================
# SECTION 6 — FIRMICUTES:BACTEROIDOTA RATIO
# =============================================================================

#' Calculate and visualise the Firmicutes:Bacteroidota ratio per sample.
#'
#' This ratio is widely reported in gut microbiome studies as a marker of
#' dysbiosis, though its clinical significance remains debated. It was
#' designed for gut data and is not meaningful for other environments.
#'
#' Samples in which either phylum is absent have no defined ratio and are
#' excluded from the plot and the tests.
#'
#' Group comparisons are non-parametric and run on the log2 ratio:
#'   unpaired: Kruskal-Wallis overall, Dunn's test for the pairs
#'             (Wilcoxon rank-sum when only two groups can be tested)
#'   paired:   Wilcoxon signed-rank for the pairs, on subjects sampled in both
#'             groups; Friedman overall only for a repeated-measures design
#'             (more than two groups, subjects sampled in every group)
#'
#' @param ps         A phyloseq object.
#' @param group_var  Metadata variable for group comparisons.
#' @param design     "unpaired" (default) or "paired".
#' @param pair_var   Metadata variable identifying the subject (paired only).
#' @param show_pairs Comparisons to draw as brackets, as "A vs B" keys taken
#'                   from the pairwise table. Default = NULL (none).
#' @param p_adjust   P-value adjustment method. Default = "BH".
#' @return A list: plot, ratio data frame, overall test, pairwise table,
#'         test label and number of samples without a ratio.

calculate_fb_ratio <- function(ps, group_var = NULL,
                               design     = c("unpaired", "paired"),
                               pair_var   = NULL,
                               show_pairs = NULL,
                               p_adjust   = "BH") {

  design <- match.arg(design)

  cat("=== Firmicutes:Bacteroidota ratio ===\n")

  tt <- tax_table(ps, errorIfNULL = FALSE)
  if (is.null(tt) || !"Phylum" %in% colnames(tt)) {
    stop("The F:B ratio needs a 'Phylum' column in the taxonomy table.")
  }

  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
  ps_phy <- fast_tax_glom(ps_rel, taxrank = "Phylum", NArm = FALSE)
  ps_melt <- psmelt(ps_phy) %>%
    filter(Phylum %in% c("Firmicutes", "Bacillota",    # both names used
                          "Bacteroidota", "Bacteroidetes")) %>%
    mutate(phylum_group = case_when(
      Phylum %in% c("Firmicutes", "Bacillota")          ~ "Firmicutes",
      Phylum %in% c("Bacteroidota", "Bacteroidetes")    ~ "Bacteroidota"
    ))

  if (nrow(ps_melt) == 0) {
    stop("Neither Firmicutes (Bacillota) nor Bacteroidota (Bacteroidetes) was found in this dataset.")
  }

  ratio_df <- ps_melt %>%
    group_by(Sample, phylum_group) %>%
    summarise(abund = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = phylum_group, values_from = abund,
                values_fill = 0)
  for (ph in c("Firmicutes", "Bacteroidota")) {
    if (!ph %in% colnames(ratio_df)) ratio_df[[ph]] <- 0
  }
  ratio_df <- ratio_df %>%
    mutate(
      FB_ratio = ifelse(Firmicutes > 0 & Bacteroidota > 0,
                        Firmicutes / Bacteroidota, NA),
      log_FB   = log2(FB_ratio)
    )

  # Add metadata
  meta_df <- data.frame(sample_data(ps), check.names = FALSE) %>%
    rownames_to_column("Sample")
  meta_df <- meta_df[, !duplicated(colnames(meta_df)), drop = FALSE]
  ratio_df <- left_join(ratio_df, meta_df, by = "Sample")

  n_excluded <- sum(is.na(ratio_df$FB_ratio))
  cat("  Median F:B ratio:", round(median(ratio_df$FB_ratio, na.rm = TRUE), 2), "\n")
  cat("  Samples without a ratio (one phylum absent):", n_excluded, "\n\n")

  fill_var <- if (!is.null(group_var) && group_var %in% colnames(ratio_df)) group_var else NULL
  if (!is.null(pair_var) && (!nzchar(pair_var) || !pair_var %in% colnames(ratio_df))) {
    pair_var <- NULL
  }

  stats <- fb_ratio_tests(ratio_df, fill_var, design = design,
                          pair_var = pair_var, p_adjust = p_adjust)

  p <- plot_fb_ratio(ratio_df, fill_var, pairwise = stats$pairwise,
                     show_pairs = show_pairs, test_label = stats$test_label)

  return(list(plot = p, data = ratio_df, group_var = fill_var,
              overall = stats$overall, pairwise = stats$pairwise,
              test_label = stats$test_label, design = design,
              n_excluded = n_excluded))
}


#' Format a p-value with three decimals for plots and tables.
#'
#' @param p        Numeric vector of p-values.
#' @param prefix   If TRUE returns "p = 0.012", otherwise "0.012".
#' @return A character vector.

fb_format_p <- function(p, prefix = FALSE) {
  out <- ifelse(is.na(p), "",
         ifelse(p < 0.001,
                if (prefix) "p < 0.001" else "< 0.001",
                paste0(if (prefix) "p = " else "", formatC(p, format = "f", digits = 3))))
  out
}


#' Non-parametric group comparisons of the log2 F:B ratio.
#'
#' @param ratio_df  Data frame returned by calculate_fb_ratio()$data.
#' @param group_var Grouping variable (column of ratio_df), or NULL.
#' @param design    "unpaired" or "paired".
#' @param pair_var  Subject variable (paired only).
#' @param p_adjust  P-value adjustment method.
#' @param min_n     Minimum samples per group for an unpaired test. Default = 3.
#' @param min_pairs Minimum shared subjects for a paired test. Default = 6.
#' @return A list: overall (test, statistic, p, n, note), pairwise table and
#'         test_label.

fb_ratio_tests <- function(ratio_df, group_var = NULL,
                           design    = c("unpaired", "paired"),
                           pair_var  = NULL,
                           p_adjust  = "BH",
                           min_n     = 3,
                           min_pairs = 6) {

  design <- match.arg(design)
  empty  <- list(overall = NULL, pairwise = NULL, test_label = NULL)
  if (is.null(group_var)) return(empty)

  df <- ratio_df[!is.na(ratio_df$log_FB) & !is.na(ratio_df[[group_var]]), , drop = FALSE]
  if (nrow(df) == 0) return(empty)
  lev <- if (is.factor(df[[group_var]])) levels(droplevels(df[[group_var]])) else
    sort(unique(as.character(df[[group_var]])))
  df$.grp <- factor(as.character(df[[group_var]]), levels = lev)
  if (length(lev) < 2) return(empty)

  pairs <- t(utils::combn(lev, 2))
  pw <- data.frame(group1 = pairs[, 1], group2 = pairs[, 2],
                   stringsAsFactors = FALSE)
  pw$key <- paste(pw$group1, "vs", pw$group2)
  pw$n1 <- NA_integer_; pw$n2 <- NA_integer_; pw$n_pairs <- NA_integer_
  pw$p  <- NA_real_;    pw$p_adj <- NA_real_
  pw$testable <- FALSE; pw$note <- ""
  overall <- list(test = NA_character_, statistic = NA_real_, p = NA_real_,
                  n = NA_integer_, note = "")

  if (design == "unpaired") {

    n_grp <- table(df$.grp)
    pw$n1 <- as.integer(n_grp[pw$group1])
    pw$n2 <- as.integer(n_grp[pw$group2])
    pw$testable <- pw$n1 >= min_n & pw$n2 >= min_n
    pw$note[!pw$testable] <- "too few samples to test"
    ok_lev <- names(n_grp)[n_grp >= min_n]

    if (length(ok_lev) == 2) {
      d2  <- droplevels(df[df$.grp %in% ok_lev, , drop = FALSE])
      wt  <- suppressWarnings(stats::wilcox.test(log_FB ~ .grp, data = d2))
      idx <- which(pw$testable)
      pw$p[idx] <- wt$p.value
      overall <- list(test = "Wilcoxon rank-sum test",
                      statistic = unname(wt$statistic), p = wt$p.value,
                      n = nrow(d2), note = "")
      pair_test <- "Wilcoxon rank-sum test"
    } else if (length(ok_lev) > 2) {
      d2 <- droplevels(df[df$.grp %in% ok_lev, , drop = FALSE])
      kt <- stats::kruskal.test(log_FB ~ .grp, data = d2)
      overall <- list(test = "Kruskal-Wallis test",
                      statistic = unname(kt$statistic), p = kt$p.value,
                      n = nrow(d2), note = "")
      dt <- rstatix::dunn_test(d2, log_FB ~ .grp, p.adjust.method = "none")
      k1 <- paste(dt$group1, "vs", dt$group2)
      k2 <- paste(dt$group2, "vs", dt$group1)
      m  <- match(pw$key, k1); m[is.na(m)] <- match(pw$key, k2)[is.na(m)]
      pw$p[pw$testable] <- dt$p[m[pw$testable]]
      pair_test <- "Dunn's test"
    } else {
      overall$note <- paste0("No test: fewer than two groups have at least ",
                             min_n, " samples.")
      pair_test <- "Dunn's test"
    }
    if (length(ok_lev) < length(lev) && length(ok_lev) >= 2) {
      overall$note <- paste0("Groups with fewer than ", min_n,
                             " samples are left out of the tests.")
    }
    test_label <- if (length(ok_lev) == 2) {
      "Unpaired samples: Wilcoxon rank-sum test"
    } else {
      paste0("Unpaired samples: Kruskal-Wallis test, then ", pair_test,
             " for the pairs")
    }

  } else {

    if (is.null(pair_var)) {
      overall$note <- "Choose the variable that identifies the subject."
      pw$note <- "no subject variable chosen"
      return(list(overall = overall, pairwise = pw,
                  test_label = "Paired samples: Wilcoxon signed-rank test"))
    }
    df$.subj <- as.character(df[[pair_var]])
    df <- df[!is.na(df$.subj), , drop = FALSE]
    # one value per subject and group (mean when a subject has several samples)
    wide <- tapply(df$log_FB, list(df$.subj, df$.grp), mean)

    for (i in seq_len(nrow(pw))) {
      x <- wide[, pw$group1[i]]; y <- wide[, pw$group2[i]]
      both <- !is.na(x) & !is.na(y)
      pw$n_pairs[i] <- sum(both)
      pw$n1[i] <- sum(both); pw$n2[i] <- sum(both)
      if (sum(both) >= min_pairs) {
        wt <- suppressWarnings(stats::wilcox.test(x[both], y[both], paired = TRUE))
        pw$p[i] <- wt$p.value
        pw$testable[i] <- !is.na(wt$p.value)
        if (is.na(wt$p.value)) pw$note[i] <- "no difference between the paired values"
      } else {
        pw$note[i] <- "too few samples to test"
      }
    }

    complete <- wide[stats::complete.cases(wide), , drop = FALSE]
    if (length(lev) == 2) {
      overall <- list(test = "Wilcoxon signed-rank test",
                      statistic = NA_real_, p = pw$p[1], n = pw$n_pairs[1],
                      note = if (!pw$testable[1]) pw$note[1] else "")
    } else if (nrow(complete) >= min_pairs) {
      ft <- stats::friedman.test(complete)
      overall <- list(test = "Friedman test (repeated measures)",
                      statistic = unname(ft$statistic), p = ft$p.value,
                      n = nrow(complete),
                      note = paste0("Uses only the ", nrow(complete),
                                    " subjects sampled in every group."))
    } else {
      overall$note <- paste0(
        "No overall test: the Friedman test needs a repeated-measures design ",
        "(the same subjects sampled in every group); only ", nrow(complete),
        " subject(s) meet this. Read the pairwise tests instead.")
    }
    test_label <- "Paired samples: Wilcoxon signed-rank test for the pairs"
  }

  n_tested <- sum(pw$testable)
  if (n_tested > 0) {
    pw$p_adj[pw$testable] <- stats::p.adjust(pw$p[pw$testable], method = p_adjust)
  }
  test_label <- paste0(test_label, " | ", p_adjust, "-adjusted p across ",
                       n_tested, " comparison", if (n_tested != 1) "s")

  list(overall = overall, pairwise = pw, test_label = test_label)
}


#' Boxplot of the log2 F:B ratio, with optional comparison brackets.
#'
#' @param ratio_df   Data frame returned by calculate_fb_ratio()$data.
#' @param group_var  Grouping variable, or NULL.
#' @param pairwise   Pairwise table from fb_ratio_tests(), or NULL.
#' @param show_pairs Keys ("A vs B") of the comparisons to draw.
#' @param test_label Text describing the test, added to the caption.
#' @return A ggplot object.

plot_fb_ratio <- function(ratio_df, group_var = NULL, pairwise = NULL,
                          show_pairs = NULL, test_label = NULL) {

  df <- ratio_df[!is.na(ratio_df$log_FB), , drop = FALSE]
  if (!is.null(group_var)) {
    df  <- df[!is.na(df[[group_var]]), , drop = FALSE]
    lev <- if (is.factor(df[[group_var]])) levels(droplevels(df[[group_var]])) else
      sort(unique(as.character(df[[group_var]])))
    df$.grp <- factor(as.character(df[[group_var]]), levels = lev)
  } else {
    lev <- "All"
    df$.grp <- factor("All")
  }
  if (nrow(df) == 0) stop("No sample has both Firmicutes and Bacteroidota, so no ratio can be shown.")

  n_lab <- table(df$.grp)
  x_lab <- stats::setNames(paste0(names(n_lab), "\n(n=", as.integer(n_lab), ")"),
                           names(n_lab))

  caption <- "Note: clinical significance of F:B ratio remains debated"
  if (!is.null(test_label)) {
    caption <- paste0(gsub(" | ", "\n", test_label, fixed = TRUE), "\n", caption)
  }

  p <- ggplot(df, aes(x = .grp, y = log_FB, fill = .grp)) +
    geom_boxplot(alpha = 0.75, outlier.shape = NA) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.7) +
    annotate("text", x = Inf, y = 0, vjust = -0.5,
             label = "F:B = 1", hjust = 1.1, size = 3, colour = "grey50") +
    scale_fill_group() +
    scale_x_discrete(labels = x_lab) +
    labs(
      title    = "Firmicutes : Bacteroidota ratio",
      subtitle = "log2 ratio | 0 = equal abundance | >0 = Firmicutes dominant",
      x        = if (!is.null(group_var)) group_var else "",
      y        = "log2(F:B ratio)",
      fill     = group_var,
      caption  = caption
    ) +
    theme_microbiome() +
    theme(legend.position = "none")

  if (length(lev) > 4) {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  # Brackets for the chosen comparisons (adjusted p-value)
  if (!is.null(pairwise) && length(show_pairs) > 0) {
    br <- pairwise[pairwise$key %in% show_pairs & pairwise$testable, , drop = FALSE]
    if (nrow(br) > 0) {
      br$x1 <- match(br$group1, lev); br$x2 <- match(br$group2, lev)
      br <- br[order(abs(br$x2 - br$x1), pmin(br$x1, br$x2)), , drop = FALSE]
      rng  <- range(df$log_FB)
      span <- if (diff(rng) > 0) diff(rng) else 1
      step <- span * 0.09
      tick <- span * 0.015
      br$y <- max(rng) + step * seq_len(nrow(br))
      br$label <- fb_format_p(br$p_adj, prefix = TRUE)
      for (i in seq_len(nrow(br))) {
        p <- p +
          annotate("segment", x = br$x1[i], xend = br$x2[i],
                   y = br$y[i], yend = br$y[i], linewidth = 0.4) +
          annotate("segment", x = br$x1[i], xend = br$x1[i],
                   y = br$y[i], yend = br$y[i] - tick, linewidth = 0.4) +
          annotate("segment", x = br$x2[i], xend = br$x2[i],
                   y = br$y[i], yend = br$y[i] - tick, linewidth = 0.4) +
          annotate("text", x = (br$x1[i] + br$x2[i]) / 2, y = br$y[i],
                   label = br$label[i], vjust = -0.4, size = 3.3)
      }
      p <- p + expand_limits(y = max(br$y) + step * 0.7)
    }
  }

  p
}


# =============================================================================
# SECTION 7 — CORE MICROBIOME
# =============================================================================

#' Identify the core microbiome — taxa present in a high fraction of samples.
#'
#' @param ps              A phyloseq object.
#' @param prevalence_cuts Vector of prevalence thresholds to evaluate. Default = c(0.5, 0.75, 0.9).
#' @param min_abundance   Minimum relative abundance to count as "present". Default = 0.001.
#' @param group_var       Optional metadata variable for group-specific core.
#' @return A list: plot, core taxa list, and summary table.

identify_core_microbiome <- function(ps,
                                      prevalence_cuts = c(0.5, 0.75, 0.9),
                                      min_abundance   = 0.001,
                                      group_var       = NULL) {

  cat("=== Core microbiome identification ===\n")

  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
  ps_gen <- fast_tax_glom(ps_rel, taxrank = "Genus", NArm = FALSE)
  gen_names <- as.character(tax_table(ps_gen)[, "Genus"])
  gen_names[is.na(gen_names) | gen_names == "NA"] <- "Unknown_Genus"
  taxa_names(ps_gen) <- make.unique(gen_names, sep = "_")

  compute_core <- function(ps_sub, label = "All") {
    otu_mat <- as.matrix(otu_table(ps_sub))
    if (!taxa_are_rows(ps_sub)) otu_mat <- t(otu_mat)

    prev_df <- data.frame(
      taxon      = rownames(otu_mat),
      prevalence = rowSums(otu_mat >= min_abundance) / ncol(otu_mat),
      mean_abund = rowMeans(otu_mat),
      group      = label,
      stringsAsFactors = FALSE
    )
    return(prev_df)
  }

  if (!is.null(group_var) && group_var %in% sample_variables(ps)) {
    groups    <- unique(as.character(sample_data(ps_gen)[[group_var]]))
    prev_list <- lapply(groups, function(g) {
      keep   <- as.character(sample_data(ps_gen)[[group_var]]) == g
      ps_sub <- prune_samples(keep, ps_gen)
      compute_core(ps_sub, label = g)
    })
    prev_df <- bind_rows(prev_list)
  } else {
    prev_df <- compute_core(ps_gen, label = "All")
  }

  # Core at each threshold
  core_summary <- lapply(prevalence_cuts, function(cut) {
    core <- prev_df %>%
      filter(prevalence >= cut) %>%
      group_by(group) %>%
      summarise(n_core_taxa = n(), .groups = "drop") %>%
      mutate(prevalence_threshold = cut)
    core
  }) %>% bind_rows()

  cat("  Core taxa at different thresholds:\n")
  print(core_summary)
  cat("\n")

  # Most stringent core
  strictest_core <- prev_df %>%
    filter(prevalence >= max(prevalence_cuts)) %>%
    arrange(desc(mean_abund)) %>%
    pull(taxon) %>% unique()

  cat("  Core taxa (>= ", max(prevalence_cuts) * 100, "% prevalence): ",
      length(strictest_core), "\n\n", sep = "")

  # --- Plot: prevalence vs abundance scatter --------------------------------
  prev_plot_df <- prev_df %>%
    mutate(
      is_core = prevalence >= 0.75,
      log_abund = log10(mean_abund + 1e-5)
    )

  p <- ggplot(prev_plot_df, aes(x = prevalence, y = log_abund)) +
    {
      if (!is.null(group_var) && group_var %in% colnames(prev_plot_df)) {
        list(
          geom_point(aes(colour = group, shape = is_core), alpha = 0.6, size = 1.8),
          facet_wrap(~ group)
        )
      } else {
        list(
          geom_point(aes(colour = is_core), alpha = 0.6, size = 1.8),
          scale_colour_manual(values = c(`FALSE` = "#95a5a6", `TRUE` = "#e74c3c"),
                              labels = c("Peripheral", "Core (≥75%)"),
                              name   = "Status")
        )
      }
    } +
    # Threshold lines
    geom_vline(xintercept = prevalence_cuts,
               linetype = "dashed", colour = "#e74c3c", alpha = 0.6, linewidth = 0.6) +
    geom_text(data = data.frame(x = prevalence_cuts + 0.01,
                                y = max(prev_plot_df$log_abund) * 0.95,
                                label = paste0(prevalence_cuts * 100, "%")),
              aes(x = x, y = y, label = label),
              hjust = 0, size = 2.5, colour = "#e74c3c", inherit.aes = FALSE) +
    geom_text(
      data = prev_plot_df %>% filter(is_core) %>%
        arrange(desc(mean_abund)) %>% slice_head(n = 10),
      aes(label = taxon), size = 2.5, hjust = -0.1, fontface = "italic"
    ) +
    scale_x_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, 1.05)) +
    labs(
      title    = "Core microbiome — prevalence vs abundance",
      subtitle = paste0("Min abundance threshold: ", min_abundance * 100, "%"),
      x        = "Prevalence (fraction of samples)",
      y        = "Log10 mean relative abundance"
    ) +
    theme_microbiome()

  return(list(
    plot         = p,
    core_taxa    = strictest_core,
    prevalence   = prev_df,
    core_summary = core_summary
  ))
}


# =============================================================================
# SECTION 8 — COMPLETE COMPOSITION WORKFLOW WRAPPER
# =============================================================================

#' Run the complete taxonomic composition pipeline.
#'
#' @param ps          A filtered phyloseq object (from microbiome_qc.R).
#' @param group_var   Metadata variable for group comparisons.
#' @param output_dir  Directory to save all outputs.
#' @return A named list of all plots and tables.

run_composition_analysis <- function(ps,
                                      group_var  = NULL,
                                      output_dir = "composition_output") {

  cat("\n", strrep("=", 60), "\n")
  cat("  TAXONOMIC COMPOSITION PIPELINE\n")
  cat(strrep("=", 60), "\n\n")

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  results <- list()

  # --- 1. Phylum-level stacked bars ----------------------------------------
  cat("--- Plot 1: Phylum composition bars ---\n")
  ps_phy   <- agglomerate_taxa(ps, rank = "Phylum", transform = "relative",
                                top_n = 12)
  results$p_phylum_bars <- plot_composition_bars(ps_phy, rank = "Phylum",
                                                  group_var = group_var,
                                                  sort_by   = "group")
  ggsave(file.path(output_dir, "01_phylum_composition_bars.pdf"),
         results$p_phylum_bars, width = 14, height = 7)

  # --- 2. Genus-level stacked bars -----------------------------------------
  cat("--- Plot 2: Genus composition bars ---\n")
  ps_gen   <- agglomerate_taxa(ps, rank = "Genus", transform = "relative",
                                top_n = 20)
  results$p_genus_bars <- plot_composition_bars(ps_gen, rank = "Genus",
                                                 group_var = group_var,
                                                 sort_by   = "group")
  ggsave(file.path(output_dir, "02_genus_composition_bars.pdf"),
         results$p_genus_bars, width = 16, height = 8)

  # --- 3. Mean composition by group ----------------------------------------
  if (!is.null(group_var)) {
    cat("--- Plot 3: Mean phylum composition by group ---\n")
    results$p_mean_comp <- plot_mean_composition(ps_phy, rank = "Phylum",
                                                  group_var = group_var,
                                                  top_n = 8)
    ggsave(file.path(output_dir, "03_mean_phylum_by_group.pdf"),
           results$p_mean_comp, width = 10, height = 7)
  }

  # --- 4. Genus heatmap -----------------------------------------------------
  cat("--- Plot 4: Genus-level heatmap ---\n")
  results$p_heatmap <- plot_abundance_heatmap(ps, rank = "Genus",
                                               top_n     = 30,
                                               group_var = group_var,
                                               transform = "clr")
  ggsave(file.path(output_dir, "04_genus_heatmap.pdf"),
         results$p_heatmap, width = 14, height = 10)

  # --- 5. Firmicutes drill-down --------------------------------------------
  cat("--- Plot 5: Firmicutes drill-down ---\n")
  tryCatch({
    results$p_firmicutes <- plot_phylum_drilldown(ps, phylum = "Firmicutes",
                                                   group_var = group_var,
                                                   top_n = 15)
    ggsave(file.path(output_dir, "05_firmicutes_drilldown.pdf"),
           results$p_firmicutes, width = 14, height = 8)
  }, error = function(e) cat("  Skipping Firmicutes drill-down:", e$message, "\n"))

  # --- 6. F:B ratio ---------------------------------------------------------
  cat("--- Plot 6: F:B ratio ---\n")
  fb_res <- calculate_fb_ratio(ps, group_var = group_var)
  results$p_fb_ratio <- fb_res$plot
  results$fb_data    <- fb_res$data
  ggsave(file.path(output_dir, "06_fb_ratio.pdf"),
         results$p_fb_ratio, width = 8, height = 6)

  # --- 7. Core microbiome ---------------------------------------------------
  cat("--- Plot 7: Core microbiome ---\n")
  core_res <- identify_core_microbiome(ps, group_var = group_var)
  results$p_core    <- core_res$plot
  results$core_taxa <- core_res$core_taxa
  ggsave(file.path(output_dir, "07_core_microbiome.pdf"),
         results$p_core, width = 10, height = 7)

  # --- 8. Abundance tables --------------------------------------------------
  cat("--- Table 1: Abundance summary tables ---\n")
  abundance_tables <- make_abundance_table(ps, rank = "Genus",
                                            group_var = group_var)
  write.csv(abundance_tables$overall,
            file.path(output_dir, "abundance_table_overall.csv"),
            row.names = FALSE)
  if (!is.null(abundance_tables$by_group)) {
    write.csv(abundance_tables$by_group,
              file.path(output_dir, "abundance_table_by_group.csv"),
              row.names = FALSE)
  }
  results$abundance_tables <- abundance_tables

  cat("\n", strrep("=", 60), "\n")
  cat("  COMPOSITION PIPELINE COMPLETE\n")
  cat(strrep("=", 60), "\n")
  cat("  Output saved to:", output_dir, "\n")
  cat("  Plots:  7 PDF files\n")
  cat("  Tables: abundance_table_overall.csv\n\n")

  return(invisible(results))
}


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# Load filtered phyloseq from QC module
# ps <- readRDS("qc_output/phyloseq_qc_filtered.rds")

# --- Option A: Complete pipeline --------------------------------------------
# results <- run_composition_analysis(
#   ps         = ps,
#   group_var  = "disease_status",
#   output_dir = "results/composition"
# )

# --- Option B: Individual plots ---------------------------------------------
# ps_phy     <- agglomerate_taxa(ps, rank = "Phylum", transform = "relative")
# p_bars     <- plot_composition_bars(ps_phy, rank = "Phylum",
#                                     group_var = "disease_status")
# p_heatmap  <- plot_abundance_heatmap(ps, rank = "Genus", group_var = "disease_status")
# core_res   <- identify_core_microbiome(ps, group_var = "disease_status")
# core_res$core_taxa
# core_res$plot
