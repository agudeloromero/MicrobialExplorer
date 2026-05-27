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
  ps_agg <- tax_glom(ps, taxrank = rank, NArm = FALSE)

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

  # Add group annotation bar if group_var provided
  if (!is.null(group_var) && group_var %in% colnames(ps_melt)) {
    group_df <- ps_melt %>%
      select(Sample, all_of(group_var)) %>%
      distinct()

    p_group <- ggplot(group_df, aes(x = Sample, y = 1,
                                    fill = .data[[group_var]])) +
      geom_col(width = 0.9) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_fill_brewer(palette = "Set2", name = group_var) +
      theme_void() +
      theme(
        legend.position = "right",
        axis.text.x     = element_blank()
      )

    p <- p_group / p +
      plot_layout(heights = c(0.08, 1), guides = "collect")
  }

  # Facet
  if (!is.null(facet_var) && facet_var %in% colnames(ps_melt)) {
    p <- p + facet_wrap(~ .data[[facet_var]], scales = "free_x", nrow = 1)
  }

  if (!show_legend) p <- p + theme(legend.position = "none")

  return(p)
}


#' Plot mean composition per group as grouped bars with error bars.
#'
#' @param ps        A phyloseq object (relative abundance).
#' @param rank      Taxonomic rank.
#' @param group_var Metadata variable defining groups.
#' @param top_n     Number of top taxa to display.
#' @return A ggplot object.

plot_mean_composition <- function(ps,
                                   rank      = "Phylum",
                                   group_var = "group",
                                   top_n     = 10) {

  cat("=== Mean composition per group ===\n")

  ps_melt <- psmelt(ps) %>%
    mutate(taxon = as.character(.data[[rank]])) %>%
    mutate(taxon = as.character(taxon))

  ps_melt$taxon[is.na(ps_melt$taxon)] <- "Unknown"

  # Select top N taxa overall
  top_taxa <- ps_melt %>%
    group_by(taxon) %>%
    summarise(mean_abund = mean(Abundance)) %>%
    arrange(desc(mean_abund)) %>%
    slice_head(n = top_n) %>%
    pull(taxon)

  ps_melt <- ps_melt %>%
    filter(taxon %in% top_taxa) %>%
    mutate(taxon = factor(taxon, levels = rev(top_taxa)))

  # Summary stats per group
  summary_df <- ps_melt %>%
    group_by(.data[[group_var]], taxon) %>%
    summarise(
      mean_abund = mean(Abundance),
      se_abund   = sd(Abundance) / sqrt(n()),
      .groups    = "drop"
    )

  n_taxa   <- length(top_taxa)
  col_vals <- TAXA_PALETTE[seq_len(n_taxa)]
  names(col_vals) <- rev(top_taxa)

  p <- ggplot(summary_df,
              aes(x = .data[[group_var]], y = mean_abund, fill = taxon)) +
    geom_col(position = position_dodge(0.85), width = 0.8) +
    geom_errorbar(
      aes(ymin = mean_abund - se_abund, ymax = mean_abund + se_abund),
      position = position_dodge(0.85), width = 0.25, linewidth = 0.5
    ) +
    scale_fill_manual(values = col_vals, name = rank) +
    scale_y_continuous(labels = percent_format(accuracy = 0.1),
                       expand  = expansion(mult = c(0, 0.05))) +
    labs(
      title    = paste0("Mean ", rank, " abundance by group"),
      subtitle = paste0("Top ", top_n, " taxa | Error bars = ±1 SE"),
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
#' @param group_var   Optional metadata variable for column annotation.
#' @param transform   Transformation to apply: "clr", "log10", or "relative".
#' @param cluster_samples Whether to cluster samples. Default = TRUE.
#' @return A ggplot heatmap.

plot_abundance_heatmap <- function(ps,
                                    rank            = "Genus",
                                    top_n           = 30,
                                    group_var       = NULL,
                                    transform       = "clr",
                                    cluster_samples = TRUE,
                                    cluster_rows    = TRUE) {

  cat("=== Abundance heatmap ===\n")

  # Agglomerate
  ps_agg  <- agglomerate_taxa(ps, rank = rank, transform = transform,
                               top_n = top_n)

  otu_mat <- as.matrix(otu_table(ps_agg))
  if (!taxa_are_rows(ps_agg)) otu_mat <- t(otu_mat)

  # Remove "Other" from heatmap
  otu_mat <- otu_mat[rownames(otu_mat) != "Other", ]

  # Row (taxa) clustering
  # When a group variable is selected, cluster on group-mean profiles so
  # the dendrogram reflects which taxa differ most between groups
  row_hclust <- NULL
  if (cluster_rows && nrow(otu_mat) > 2) {
    if (!is.null(group_var)) {
      meta_df <- data.frame(sample_data(ps_agg))
      if (group_var %in% colnames(meta_df)) {
        groups     <- as.character(meta_df[[group_var]])
        # Build group-mean matrix (taxa x groups)
        grp_levels <- unique(groups)
        mean_mat   <- sapply(grp_levels, function(g) {
          idx <- which(groups == g)
          if (length(idx) == 1) otu_mat[, idx]
          else rowMeans(otu_mat[, idx, drop = FALSE])
        })
        row_hclust <- hclust(dist(mean_mat))
      } else {
        row_hclust <- hclust(dist(otu_mat))
      }
    } else {
      row_hclust <- hclust(dist(otu_mat))
    }
    otu_mat <- otu_mat[row_hclust$order, ]
  }

  # Column (sample) ordering
  if (cluster_samples && !is.null(group_var) && ncol(otu_mat) > 2) {
    # Cluster within each group separately, then concatenate groups
    meta_df <- data.frame(sample_data(ps_agg))
    if (group_var %in% colnames(meta_df)) {
      groups      <- unique(as.character(meta_df[[group_var]]))
      col_order   <- unlist(lapply(groups, function(g) {
        samps <- rownames(meta_df)[meta_df[[group_var]] == g]
        samps <- intersect(samps, colnames(otu_mat))
        if (length(samps) <= 2) return(samps)
        sub_mat   <- otu_mat[, samps, drop = FALSE]
        samps[hclust(dist(t(sub_mat)))$order]
      }))
      otu_mat <- otu_mat[, col_order, drop = FALSE]
    } else {
      # Fall back to global clustering if group_var not in metadata
      col_order <- hclust(dist(t(otu_mat)))$order
      otu_mat   <- otu_mat[, col_order]
    }
  } else if (cluster_samples && is.null(group_var) && ncol(otu_mat) > 2) {
    # Global hierarchical clustering — no group sorting, shows column dendrogram
    col_order <- hclust(dist(t(otu_mat)))$order
    otu_mat   <- otu_mat[, col_order, drop = FALSE]
  } else if (!cluster_samples && ncol(otu_mat) > 2) {
    # No clustering: sort samples by group variable if available
    if (!is.null(group_var)) {
      meta_df <- data.frame(sample_data(ps_agg))
      if (group_var %in% colnames(meta_df)) {
        col_order <- rownames(meta_df)[order(as.character(meta_df[[group_var]]))]
        col_order <- intersect(col_order, colnames(otu_mat))
        otu_mat   <- otu_mat[, col_order, drop = FALSE]
      }
    }
  }

  # Store final column order BEFORE melting — used for both plot and annotation
  sample_levels <- colnames(otu_mat)

  # Melt to long format
  heat_df <- as.data.frame(otu_mat) %>%
    rownames_to_column("taxon") %>%
    pivot_longer(-taxon, names_to = "sample", values_to = "value") %>%
    mutate(
      taxon  = factor(taxon,  levels = rownames(otu_mat)),
      sample = factor(sample, levels = sample_levels)
    )

  # Add group annotation
  if (!is.null(group_var) && group_var %in% sample_variables(ps)) {
    group_df <- data.frame(sample_data(ps)) %>%
      select(all_of(group_var)) %>%
      rownames_to_column("sample")
    heat_df <- left_join(heat_df, group_df, by = "sample")
  }

  # Colour scale midpoint
  midpoint <- median(heat_df$value)

  p <- ggplot(heat_df, aes(x = sample, y = taxon, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.15) +
    scale_fill_gradient2(
      low      = "#3498db",
      mid      = "white",
      high     = "#e74c3c",
      midpoint = midpoint,
      name     = ifelse(transform == "clr", "CLR",
                        ifelse(transform == "log10", "Log10", "Rel. abund."))
    ) +
    scale_x_discrete(limits = sample_levels, expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      title    = paste0(rank, "-level abundance heatmap"),
      subtitle = paste0("Top ", nrow(otu_mat), " taxa | Transform: ", transform),
      x        = "Samples",
      y        = rank
    ) +
    theme_microbiome() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y  = element_text(size = 8, face = "italic"),
      panel.border = element_rect(colour = "grey80", fill = NA)
    )

  # ── Collect optional components before assembly ──────────────────────────
  has_annot  <- !is.null(group_var) && group_var %in% colnames(heat_df)
  has_dendro <- !is.null(row_hclust) && requireNamespace("ggdendro", quietly = TRUE)

  # Build annotation strip
  p_annot <- NULL
  if (has_annot) {
    annot_df <- data.frame(sample_data(ps_agg)) %>%
      rownames_to_column("sample") %>%
      select(sample, all_of(group_var)) %>%
      filter(sample %in% sample_levels) %>%
      mutate(
        sample   = factor(sample, levels = sample_levels),
        annot_y  = "group",
        grp_fill = factor(as.character(.data[[group_var]]))
      )

    n_groups <- nlevels(annot_df$grp_fill)
    pal_20   <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
                  "#a65628","#f781bf","#999999","#66c2a5","#fc8d62",
                  "#8da0cb","#e78ac3","#a6d854","#ffd92f","#e5c494",
                  "#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e")

    grp_cols <- if (n_groups <= 8) {
      scale_fill_brewer(palette = "Set2", name = group_var)
    } else if (n_groups <= 20) {
      scale_fill_manual(
        values = setNames(pal_20[seq_len(n_groups)], levels(annot_df$grp_fill)),
        name   = group_var)
    } else {
      scale_fill_viridis_d(name = group_var)
    }

    p_annot <- ggplot(annot_df,
                      aes(x = sample, y = annot_y, fill = grp_fill)) +
      geom_tile() +
      grp_cols +
      scale_x_discrete(limits = sample_levels, expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      labs(fill = group_var) +
      theme_void() +
      theme(legend.position  = "right",
            plot.margin      = ggplot2::margin(0, 0, 0, 0))
  }

  # Build row dendrogram
  p_dendro <- NULL
  if (has_dendro) {
    dendro_data <- ggdendro::dendro_data(as.dendrogram(row_hclust),
                                         type = "rectangle")
    n_taxa <- nrow(otu_mat)

    p_dendro <- ggplot(ggdendro::segment(dendro_data)) +
      geom_segment(aes(x = y, y = n_taxa + 1 - x,
                       xend = yend, yend = n_taxa + 1 - xend),
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
  # Build column dendrogram (when cluster_samples=TRUE, group_var=NULL)
  p_col_dendro <- NULL
  if (cluster_samples && is.null(group_var) &&
      requireNamespace("ggdendro", quietly = TRUE) && ncol(otu_mat) > 2) {
    col_hclust   <- hclust(dist(t(otu_mat)))
    col_dd       <- ggdendro::dendro_data(as.dendrogram(col_hclust), type = "rectangle")
    n_samp       <- ncol(otu_mat)
    p_col_dendro <- ggplot(ggdendro::segment(col_dd)) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend),
                   linewidth = 0.4, colour = "#2c3e50") +
      scale_x_continuous(limits = c(0.5, n_samp + 0.5), expand = c(0, 0)) +
      scale_y_continuous(expand = c(0.02, 0)) +
      theme_void() +
      theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
  }

  # ── Assemble final layout ─────────────────────────────────────────────────
  has_col_dendro <- !is.null(p_col_dendro)
  if (has_dendro && has_annot) {
    p <- (patchwork::plot_spacer() + p_annot + p_dendro + p) +
      patchwork::plot_layout(
        design  = "AB
CD",
        widths  = c(0.15, 1),
        heights = c(0.05, 1),
        guides  = "collect"
      )
  } else if (has_col_dendro && has_dendro) {
    p <- (patchwork::plot_spacer() + p_col_dendro + p_dendro + p) +
      patchwork::plot_layout(
        design  = "AB
CD",
        widths  = c(0.15, 1),
        heights = c(0.15, 1),
        guides  = "collect"
      )
  } else if (has_col_dendro) {
    p <- (p_col_dendro / p) +
      patchwork::plot_layout(heights = c(0.15, 1), guides = "collect")
  } else if (has_dendro) {
    p <- (p_dendro | p) +
      patchwork::plot_layout(widths = c(0.15, 1), guides = "collect")
  } else if (has_annot) {
    p <- (p_annot / p) +
      patchwork::plot_layout(heights = c(0.05, 1), guides = "collect")
  }
  return(p)
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
  ps_agg  <- tax_glom(ps_rel, taxrank = rank, NArm = FALSE)

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
  ps_gen <- tax_glom(ps_phy, taxrank = "Genus", NArm = FALSE)
  genera <- as.character(tax_table(ps_gen)[, "Genus"])
  genera[is.na(genera)] <- "Unknown genus"
  taxa_names(ps_gen) <- genera

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
      scale_fill_brewer(palette = "Set2") +
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
#' dysbiosis, though its clinical significance remains debated.
#'
#' @param ps        A phyloseq object.
#' @param group_var Metadata variable for group comparisons.
#' @return A list: plot and ratio data frame.

calculate_fb_ratio <- function(ps, group_var = NULL) {

  cat("=== Firmicutes:Bacteroidota ratio ===\n")

  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
  ps_phy <- tax_glom(ps_rel, taxrank = "Phylum", NArm = FALSE)
  ps_melt <- psmelt(ps_phy) %>%
    filter(Phylum %in% c("Firmicutes", "Bacillota",    # both names used
                          "Bacteroidota", "Bacteroidetes")) %>%
    mutate(phylum_group = case_when(
      Phylum %in% c("Firmicutes", "Bacillota")          ~ "Firmicutes",
      Phylum %in% c("Bacteroidota", "Bacteroidetes")    ~ "Bacteroidota"
    ))

  ratio_df <- ps_melt %>%
    group_by(Sample, phylum_group) %>%
    summarise(abund = sum(Abundance), .groups = "drop") %>%
    pivot_wider(names_from = phylum_group, values_from = abund,
                values_fill = 0) %>%
    mutate(
      FB_ratio = ifelse(Bacteroidota > 0,
                        Firmicutes / Bacteroidota, NA),
      log_FB   = log2(FB_ratio + 0.001)
    )

  # Add metadata
  meta_df <- data.frame(sample_data(ps)) %>% rownames_to_column("Sample")
  ratio_df <- left_join(ratio_df, meta_df, by = "Sample")

  cat("  Median F:B ratio:", round(median(ratio_df$FB_ratio, na.rm = TRUE), 2), "\n\n")

  fill_var <- if (!is.null(group_var) && group_var %in% colnames(ratio_df)) group_var else NULL

  p <- ggplot(ratio_df,
              aes(x = if (!is.null(fill_var)) .data[[fill_var]] else "All",
                  y = log_FB,
                  fill = if (!is.null(fill_var)) .data[[fill_var]] else "All")) +
    geom_boxplot(alpha = 0.75, outlier.size = 1) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.7) +
    annotate("text", x = 0.55, y = 0.1,
             label = "F:B = 1", hjust = 0, size = 3, colour = "grey50") +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title    = "Firmicutes : Bacteroidota ratio",
      subtitle = "log2 ratio | 0 = equal abundance | >0 = Firmicutes dominant",
      x        = if (!is.null(fill_var)) fill_var else "",
      y        = "log2(F:B ratio)",
      fill     = fill_var,
      caption  = "Note: clinical significance of F:B ratio remains debated"
    ) +
    theme_microbiome()

  return(list(plot = p, data = ratio_df))
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
  ps_gen <- tax_glom(ps_rel, taxrank = "Genus", NArm = FALSE)
  taxa_names(ps_gen) <- as.character(tax_table(ps_gen)[, "Genus"])

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
    annotate("text", x = prevalence_cuts + 0.01,
             y = max(prev_plot_df$log_abund) * 0.95,
             label = paste0(prevalence_cuts * 100, "%"),
             hjust = 0, size = 2.5, colour = "#e74c3c") +
    # Label top core taxa
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
