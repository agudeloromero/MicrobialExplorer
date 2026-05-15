# =============================================================================
# MICROBIOME BETA DIVERSITY ANALYSIS
# =============================================================================
# Description : Comprehensive beta diversity analysis pipeline
# Input       : Filtered phyloseq object (output from microbiome_qc.R)
# Output      : Ordination plots, distance matrices, statistical tests,
#               dispersion analysis, and differential abundance results
# Author      : Patricia
# Dependencies: phyloseq, vegan, ggplot2, dplyr, tidyr, patchwork,
#               ape, rstatix, ggpubr, scales, tibble, RColorBrewer
# =============================================================================

# --- 1. LOAD LIBRARIES -------------------------------------------------------

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ape)
  library(rstatix)
  library(scales)
  library(tibble)
  library(RColorBrewer)
  library(stringr)
})

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
# SECTION 1 — COMPUTE DISTANCE MATRICES
# =============================================================================

#' Compute one or more beta diversity distance matrices from a phyloseq object.
#'
#' Supported distances:
#'   - bray          : Bray-Curtis dissimilarity (abundance-based)
#'   - jaccard       : Jaccard distance (presence/absence)
#'   - unifrac       : Unweighted UniFrac (requires tree)
#'   - wunifrac      : Weighted UniFrac   (requires tree)
#'   - aitchison     : Aitchison distance (CLR-transformed Euclidean)
#'   - euclidean     : Euclidean distance on CLR-transformed data
#'
#' @param ps          A phyloseq object (raw counts recommended).
#' @param methods     Character vector of distance methods. Default = c("bray","jaccard").
#' @param rarefaction Whether to rarefy before computing distances. Default = TRUE.
#' @param rare_depth  Rarefaction depth. Default = minimum sample depth.
#' @param seed        Random seed. Default = 42.
#' @return A named list of dist objects.

compute_distances <- function(ps,
                               methods     = c("bray", "jaccard"),
                               rarefaction = TRUE,
                               rare_depth  = NULL,
                               seed        = 42) {

  cat("=== Computing beta diversity distances ===\n")

  # --- Rarefaction ----------------------------------------------------------
  if (rarefaction) {
    if (is.null(rare_depth)) rare_depth <- min(sample_sums(ps))
    cat("  Rarefying to depth:", rare_depth, "\n")
    set.seed(seed)
    ps <- rarefy_even_depth(ps, sample.size = rare_depth,
                             rngseed = seed, replace = FALSE,
                             verbose = FALSE)
  }

  otu_mat <- as.matrix(otu_table(ps))
  if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  otu_t   <- matrix(t(otu_mat), nrow=ncol(otu_mat), ncol=nrow(otu_mat),
                   dimnames=list(colnames(otu_mat), rownames(otu_mat)))   # samples as rows for vegan

  dist_list <- list()

  for (method in methods) {

    cat("  Computing:", method, "... ")

    d <- tryCatch({

      if (method == "bray") {
        vegdist(otu_t, method = "bray")

      } else if (method == "jaccard") {
        vegdist(otu_t, method = "jaccard", binary = TRUE)

      } else if (method == "unifrac") {
        if (is.null(phy_tree(ps, errorIfNULL = FALSE))) {
          cat("SKIPPED (no tree)\n"); next
        }
        UniFrac(ps, weighted = FALSE, normalized = TRUE)

      } else if (method == "wunifrac") {
        if (is.null(phy_tree(ps, errorIfNULL = FALSE))) {
          cat("SKIPPED (no tree)\n"); next
        }
        UniFrac(ps, weighted = TRUE, normalized = TRUE)

      } else if (method %in% c("aitchison", "euclidean")) {
        # CLR transformation with pseudocount
        clr_mat <- log(otu_t + 0.5) -
          rowMeans(log(otu_t + 0.5))
        dist(clr_mat, method = "euclidean")

      } else {
        vegdist(otu_t, method = method)
      }

    }, error = function(e) {
      cat("ERROR:", e$message, "\n")
      return(NULL)
    })

    if (!is.null(d)) {
      dist_list[[method]] <- d
      cat("done\n")
    }
  }

  cat("  Distance matrices computed:", paste(names(dist_list), collapse = ", "), "\n\n")
  return(dist_list)
}


# =============================================================================
# SECTION 2 — ORDINATION
# =============================================================================

#' Run ordination on a distance matrix.
#'
#' @param dist_obj    A dist object.
#' @param method      Ordination method: "PCoA", "NMDS", or "PCA".
#' @param k           Number of dimensions to compute. Default = 3.
#' @return A list: ordination result and variance explained (if applicable).

run_ordination <- function(dist_obj,
                            method = "PCoA",
                            k      = 3) {

  cat("  Running", method, "ordination...\n")

  if (method == "PCoA") {
    ord     <- cmdscale(dist_obj, k = k, eig = TRUE)
    eigenvalues <- ord$eig
    # Keep only positive eigenvalues for variance calculation
    pos_eig <- eigenvalues[eigenvalues > 0]
    var_exp <- round(pos_eig / sum(pos_eig) * 100, 1)
    coords  <- as.data.frame(ord$points)
    colnames(coords) <- paste0("PC", seq_len(ncol(coords)))

    return(list(
      coords     = coords,
      eigenvalues = eigenvalues,
      var_exp    = var_exp,
      method     = "PCoA"
    ))

  } else if (method == "NMDS") {
    set.seed(42)
    ord <- suppressMessages(
      metaMDS(dist_obj, k = 2, trymax = 100, trace = FALSE)
    )
    coords <- as.data.frame(ord$points)
    colnames(coords) <- c("NMDS1", "NMDS2")

    cat("    NMDS stress:", round(ord$stress, 4))
    if (ord$stress > 0.2) cat(" WARNING: High stress - interpret with caution")
    cat("\n")

    return(list(
      coords  = coords,
      stress  = ord$stress,
      method  = "NMDS",
      var_exp = NULL
    ))

  } else if (method == "PCA") {
    pca     <- prcomp(dist_obj, scale. = TRUE)
    var_exp <- round(summary(pca)$importance[2, ] * 100, 1)
    coords  <- as.data.frame(pca$x[, seq_len(k)])

    return(list(
      coords  = coords,
      var_exp = var_exp,
      method  = "PCA"
    ))
  }
}


# =============================================================================
# SECTION 3 — ORDINATION PLOTS
# =============================================================================

#' Plot ordination results with metadata overlays and ellipses.
#'
#' @param ord_result  Output from run_ordination().
#' @param ps          A phyloseq object (for metadata access).
#' @param group_var   Primary grouping variable for colours.
#' @param shape_var   Optional second variable for point shapes.
#' @param label_var   Optional variable for labelling points.
#' @param ellipse     Whether to draw 95% confidence ellipses. Default = TRUE.
#' @param axes        Which axes to plot. Default = c(1, 2).
#' @param point_size  Size of points. Default = 3.
#' @param spider      Whether to draw spider lines to group centroids. Default = FALSE.
#' @return A ggplot object.

plot_ordination <- function(ord_result,
                             ps,
                             group_var  = "group",
                             shape_var  = NULL,
                             label_var  = NULL,
                             ellipse    = TRUE,
                             axes       = c(1, 2),
                             point_size = 3,
                             spider     = FALSE) {

  cat("=== Ordination plot ===\n")

  # Combine coordinates with metadata
  coords  <- ord_result$coords
  meta_df <- data.frame(sample_data(ps)) %>% rownames_to_column("sample")

  plot_df <- coords %>%
    rownames_to_column("sample") %>%
    left_join(meta_df, by = "sample")

  # Select axes
  axis_cols <- colnames(coords)[axes]
  x_col     <- axis_cols[1]
  y_col     <- axis_cols[2]

  # Axis labels
  if (!is.null(ord_result$var_exp)) {
    x_lab <- paste0(x_col, " (", ord_result$var_exp[axes[1]], "%)")
    y_lab <- paste0(y_col, " (", ord_result$var_exp[axes[2]], "%)")
  } else if (ord_result$method == "NMDS") {
    x_lab <- paste0("NMDS1 | Stress = ", round(ord_result$stress, 3))
    y_lab <- "NMDS2"
  } else {
    x_lab <- x_col
    y_lab <- y_col
  }

  groups   <- unique(plot_df[[group_var]])
  n_groups <- length(groups)
  colours  <- brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)]
  names(colours) <- groups

  # Spider lines (centroids)
  if (spider) {
    centroids <- plot_df %>%
      group_by(.data[[group_var]]) %>%
      summarise(
        cx = mean(.data[[x_col]]),
        cy = mean(.data[[y_col]]),
        .groups = "drop"
      )
    plot_df <- left_join(plot_df, centroids, by = group_var)
  }

  p <- ggplot(plot_df, aes(x = .data[[x_col]], y = .data[[y_col]]))

  # Spider lines
  if (spider) {
    p <- p + geom_segment(
      aes(xend = cx, yend = cy, colour = .data[[group_var]]),
      alpha = 0.3, linewidth = 0.4
    )
  }

  # Ellipses
  if (ellipse && n_groups > 1) {
    p <- p + stat_ellipse(
      aes(colour = .data[[group_var]], fill = .data[[group_var]]),
      geom  = "polygon", alpha = 0.08, level = 0.95,
      linewidth = 0.6
    )
  }

  # Points
  p <- p + {
    if (!is.null(shape_var) && shape_var %in% colnames(plot_df)) {
      geom_point(aes(colour = .data[[group_var]],
                     shape  = .data[[shape_var]]),
                 size = point_size, alpha = 0.85)
    } else {
      geom_point(aes(colour = .data[[group_var]]),
                 size = point_size, alpha = 0.85)
    }
  }

  # Labels
  if (!is.null(label_var) && label_var %in% colnames(plot_df)) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = .data[[label_var]]),
      size = 2.5, max.overlaps = 20
    )
  }

  p <- p +
    scale_colour_manual(values = colours, name = group_var) +
    scale_fill_manual(values = colours, name = group_var, guide = "none") +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey80", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey80", linewidth = 0.4) +
    labs(
      title    = paste0(ord_result$method, " - ", group_var),
      subtitle = paste0(nrow(plot_df), " samples | ",
                        n_groups, " groups"),
      x        = x_lab,
      y        = y_lab
    ) +
    theme_microbiome() +
    coord_fixed()

  return(p)
}


#' Plot all three axis pairs from a PCoA (PC1-2, PC1-3, PC2-3) in a panel.
#'
#' @param ord_result  Output from run_ordination() with method = "PCoA".
#' @param ps          A phyloseq object.
#' @param group_var   Grouping variable.
#' @return A patchwork of three ordination plots.

plot_pcoa_panel <- function(ord_result, ps, group_var = "group") {

  p12 <- plot_ordination(ord_result, ps, group_var, axes = c(1, 2), ellipse = TRUE)
  p13 <- plot_ordination(ord_result, ps, group_var, axes = c(1, 3), ellipse = TRUE)
  p23 <- plot_ordination(ord_result, ps, group_var, axes = c(2, 3), ellipse = TRUE)

  combined <- (p12 | p13 | p23) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = paste0("PCoA panel — ", group_var),
      subtitle = paste0("Variance explained: PC1 = ", ord_result$var_exp[1],
                        "%, PC2 = ", ord_result$var_exp[2],
                        "%, PC3 = ", ord_result$var_exp[3], "%"),
      theme    = theme(plot.title    = element_text(size = 14, face = "bold"),
                       plot.subtitle = element_text(size = 10, colour = "grey40"))
    )

  return(combined)
}


# =============================================================================
# SECTION 4 — PERMANOVA
# =============================================================================

#' Run PERMANOVA (adonis2) to test group differences in community composition.
#'
#' @param dist_obj    A dist object.
#' @param ps          A phyloseq object (for metadata access).
#' @param formula_rhs Right-hand side of the formula as a string.
#'                    E.g. "group", "group + age + sex", "group * timepoint".
#' @param permutations Number of permutations. Default = 999.
#' @param strata      Optional stratification variable (for paired designs).
#' @return A list: PERMANOVA result, pairwise results, and betadisper test.

run_permanova <- function(dist_obj,
                           ps,
                           formula_rhs  = "group",
                           permutations = 999,
                           strata       = NULL) {

  cat("=== PERMANOVA ===\n")

  meta_df   <- data.frame(sample_data(ps))

  # Ensure sample order matches distance matrix
  dist_samples <- attr(dist_obj, "Labels")
  if (!is.null(dist_samples)) {
    meta_df <- meta_df[dist_samples, , drop = FALSE]
  }

  formula_full <- as.formula(paste("dist_obj ~", formula_rhs))

  # --- Main PERMANOVA -------------------------------------------------------
  set.seed(42)
  perm_result <- adonis2(
    formula      = formula_full,
    data         = meta_df,
    permutations = permutations,
    by           = "margin",
    strata       = if (!is.null(strata)) meta_df[[strata]] else NULL
  )

  cat("  PERMANOVA results:\n")
  print(perm_result)
  cat("\n")

  # --- Pairwise PERMANOVA ---------------------------------------------------
  # Extract first variable from formula for pairwise tests
  first_var    <- trimws(strsplit(formula_rhs, "\\+|\\*")[[1]][1])
  group_levels <- unique(meta_df[[first_var]])
  n_groups     <- length(group_levels)

  pairwise_results <- NULL

  if (n_groups > 2) {
    cat("  Running pairwise PERMANOVA (", first_var, ")...\n")

    pairs     <- combn(group_levels, 2, simplify = FALSE)
    pair_list <- lapply(pairs, function(pair) {
      idx      <- meta_df[[first_var]] %in% pair
      sub_dist <- as.dist(as.matrix(dist_obj)[idx, idx])
      sub_meta <- meta_df[idx, , drop = FALSE]

      set.seed(42)
      res <- adonis2(
        formula      = as.formula(paste("sub_dist ~", first_var)),
        data         = sub_meta,
        permutations = permutations,
        by           = "margin"
      )

      data.frame(
        group1    = pair[1],
        group2    = pair[2],
        R2        = round(res$R2[1], 4),
        F_stat    = round(res$F[1], 3),
        p_value   = res$`Pr(>F)`[1],
        stringsAsFactors = FALSE
      )
    })

    pairwise_results <- bind_rows(pair_list) %>%
      mutate(p_adjusted = p.adjust(p_value, method = "BH"),
             significance = case_when(
               p_adjusted < 0.001 ~ "***",
               p_adjusted < 0.01  ~ "**",
               p_adjusted < 0.05  ~ "*",
               p_adjusted < 0.1   ~ ".",
               TRUE               ~ "ns"
             ))

    cat("  Pairwise results:\n")
    print(pairwise_results)
    cat("\n")
  }

  # --- Betadisper (homogeneity of dispersion) --------------------------------
  cat("  Testing homogeneity of dispersion (betadisper)...\n")

  betadisp_result <- tryCatch({
    bd  <- betadisper(dist_obj, meta_df[[first_var]])
    bdt <- permutest(bd, permutations = permutations)
    list(betadisper = bd, test = bdt)
  }, error = function(e) {
    cat("  betadisper error:", e$message, "\n")
    NULL
  })

  if (!is.null(betadisp_result)) {
    bd_p <- betadisp_result$test$tab$`Pr(>F)`[1]
    cat("  Betadisper p-value:", round(bd_p, 4))
    if (bd_p < 0.05) {
      cat(" WARNING: Significant dispersion differences — PERMANOVA results may",
          "reflect dispersion rather than location\n")
    } else {
      cat(" OK: No significant dispersion differences\n")
    }
    cat("\n")
  }

  return(list(
    permanova  = perm_result,
    pairwise   = pairwise_results,
    betadisper = betadisp_result
  ))
}


# =============================================================================
# SECTION 5 — BETADISPER VISUALISATION
# =============================================================================

#' Plot dispersion (distance to centroid) per group from betadisper.
#'
#' @param dist_obj  A dist object.
#' @param ps        A phyloseq object.
#' @param group_var Grouping variable.
#' @return A patchwork of betadisper PCoA and boxplot.

plot_betadisper <- function(dist_obj, ps, group_var = "group") {
  
  cat("=== Betadisper visualisation ===\n")
  
  meta_df <- data.frame(sample_data(ps))
  
  # Align metadata to distance matrix samples
  dist_mat <- as.matrix(dist_obj)
  dist_samples <- rownames(dist_mat)
  meta_df <- meta_df[dist_samples, , drop = FALSE]
  
  if (!group_var %in% colnames(meta_df)) {
    stop("group_var '", group_var, "' not found in sample metadata.")
  }
  
  groups <- factor(meta_df[[group_var]])
  names(groups) <- rownames(meta_df)
  
  if (length(groups) != attr(dist_obj, "Size")) {
    stop(
      "Distance matrix has ", attr(dist_obj, "Size"),
      " samples but group vector has ", length(groups), "."
    )
  }
  
  bd <- vegan::betadisper(dist_obj, groups)
  
  if (ncol(bd$vectors) < 2) {
    stop("Betadisper ordination has fewer than 2 axes; cannot plot PC1/PC2.")
  }
  
  bd_pcoa <- data.frame(
    sample = rownames(bd$vectors),
    PC1 = as.numeric(bd$vectors[, 1]),
    PC2 = as.numeric(bd$vectors[, 2]),
    group = groups,
    stringsAsFactors = FALSE
  )
  
  bd_centroids <- data.frame(
    group = rownames(bd$centroids),
    cx = as.numeric(bd$centroids[, 1]),
    cy = as.numeric(bd$centroids[, 2]),
    stringsAsFactors = FALSE
  )
  
  bd_df <- data.frame(
    sample = names(bd$distances),
    distance = as.numeric(bd$distances),
    group = groups,
    stringsAsFactors = FALSE
  )
  
  segment_df <- bd_pcoa %>%
    left_join(bd_centroids, by = "group")
  
  n_groups <- dplyr::n_distinct(bd_pcoa$group)
  colours <- RColorBrewer::brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)]
  names(colours) <- levels(groups)
  
  p1 <- ggplot(bd_pcoa, aes(x = PC1, y = PC2, colour = group)) +
    geom_point(alpha = 0.75, size = 2.5) +
    geom_point(
      data = bd_centroids,
      aes(x = cx, y = cy, colour = group),
      shape = 8,
      size = 5,
      stroke = 1.5,
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = segment_df,
      aes(x = PC1, y = PC2, xend = cx, yend = cy, colour = group),
      alpha = 0.25,
      linewidth = 0.4,
      inherit.aes = FALSE
    ) +
    scale_colour_manual(values = colours, name = group_var) +
    labs(
      title = "Betadisper - distance to centroid",
      subtitle = "Stars = group centroids",
      x = "PC1",
      y = "PC2"
    ) +
    theme_microbiome() +
    coord_fixed()
  
  p2 <- ggplot(bd_df, aes(x = group, y = distance, fill = group)) +
    geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.55) +
    geom_jitter(aes(colour = group), width = 0.15, alpha = 0.5, size = 1.8) +
    scale_fill_manual(values = colours, guide = "none") +
    scale_colour_manual(values = colours, guide = "none") +
    labs(
      title = "Within-group dispersion",
      subtitle = "Distance of each sample to its group centroid",
      x = group_var,
      y = "Distance to centroid"
    ) +
    theme_microbiome() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  combined <- p1 | p2
  return(combined)
}

# =============================================================================
# SECTION 6 — DISTANCE MATRIX HEATMAP
# =============================================================================

#' Visualise the pairwise distance matrix as a hierarchically clustered heatmap.
#'
#' @param dist_obj  A dist object.
#' @param ps        A phyloseq object.
#' @param group_var Grouping variable for annotation strip.
#' @return A ggplot heatmap.

plot_distance_heatmap <- function(dist_obj, ps, group_var = "group") {
  
  cat("=== Distance matrix heatmap ===\n")
  
  dist_mat <- as.matrix(dist_obj)
  meta_df  <- data.frame(sample_data(ps))
  
  dist_samples <- rownames(dist_mat)
  meta_df <- meta_df[dist_samples, , drop = FALSE]
  
  if (!group_var %in% colnames(meta_df)) {
    stop("group_var '", group_var, "' not found in sample metadata.")
  }
  
  # Hierarchical clustering for ordering
  hc <- hclust(dist_obj, method = "ward.D2")
  sample_order <- hc$labels[hc$order]
  
  dist_mat <- dist_mat[sample_order, sample_order]
  
  heat_df <- as.data.frame(dist_mat) %>%
    rownames_to_column("sample1") %>%
    pivot_longer(-sample1, names_to = "sample2", values_to = "distance") %>%
    mutate(
      sample1 = factor(sample1, levels = sample_order),
      sample2 = factor(sample2, levels = rev(sample_order))
    )
  
  group_df <- meta_df %>%
    rownames_to_column("sample") %>%
    transmute(
      sample = sample,
      group1 = .data[[group_var]]
    )
  
  group_df_for_join <- group_df
  names(group_df_for_join)[names(group_df_for_join) == "sample"] <- "sample1"
  
  heat_df <- heat_df %>%
    left_join(group_df_for_join, by = "sample1")
  
  p <- ggplot(heat_df, aes(x = sample1, y = sample2, fill = distance)) +
    geom_tile() +
    scale_fill_distiller(
      palette = "RdYlBu",
      direction = -1,
      name = "Distance"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      title = "Pairwise distance matrix",
      subtitle = paste0("Hierarchical clustering - Ward.D2 | ",
                        nrow(dist_mat), " samples"),
      x = NULL,
      y = NULL
    ) +
    theme_microbiome() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      panel.border = element_rect(colour = "grey80", fill = NA)
    )
  
  return(p)
}


# =============================================================================
# SECTION 7 — ENVFIT: METADATA VECTORS ON ORDINATION
# =============================================================================

#' Fit metadata variables onto an ordination and visualise as vectors/centroids.
#'
#' @param ord_result  Output from run_ordination().
#' @param ps          A phyloseq object.
#' @param group_var   Grouping variable for point colouring.
#' @param env_vars    Character vector of metadata variables to fit.
#' @param p_threshold Significance threshold for displaying vectors. Default = 0.05.
#' @param permutations Number of permutations. Default = 999.
#' @return A list: plot and envfit results.

plot_envfit <- function(ord_result,
                         ps,
                         group_var     = "group",
                         env_vars      = NULL,
                         p_threshold   = 0.05,
                         permutations  = 999) {

  cat("=== envfit — metadata overlay on ordination ===\n")

  meta_df   <- data.frame(sample_data(ps))

  if (is.null(env_vars)) {
    # Auto-select numeric metadata variables
    env_vars <- names(meta_df)[sapply(meta_df, is.numeric)]
    cat("  Auto-detected numeric variables:", paste(env_vars, collapse = ", "), "\n")
  }

  env_vars  <- intersect(env_vars, colnames(meta_df))
  if (length(env_vars) == 0) {
    cat("  No valid numeric variables found for envfit. Skipping.\n")
    return(NULL)
  }

  env_df    <- meta_df[, env_vars, drop = FALSE]
  coords    <- ord_result$coords[, 1:2]

  # Align sample order
  shared    <- intersect(rownames(coords), rownames(env_df))
  coords    <- coords[shared, ]
  env_df    <- env_df[shared, , drop = FALSE]

  set.seed(42)
  ef        <- envfit(coords, env_df, permutations = permutations, na.rm = TRUE)

  # Extract significant vectors
  ef_df     <- as.data.frame(scores(ef, "vectors"))
  ef_df$variable <- rownames(ef_df)
  ef_df$r2       <- ef$vectors$r
  ef_df$p        <- ef$vectors$pvals

  ef_sig    <- ef_df %>% filter(p <= p_threshold)
  cat("  Significant variables (p <", p_threshold, "):",
      nrow(ef_sig), "of", nrow(ef_df), "\n\n")

  # Build ordination base plot
  axis_cols <- colnames(coords)
  x_col     <- axis_cols[1]
  y_col     <- axis_cols[2]

  plot_df   <- coords %>%
    rownames_to_column("sample") %>%
    left_join(meta_df %>% rownames_to_column("sample"), by = "sample")

  groups    <- unique(plot_df[[group_var]])
  n_groups  <- length(groups)
  colours   <- brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)]
  names(colours) <- groups

  # Scale vectors for plotting
  scale_fac <- max(abs(coords)) * 0.8
  ef_sig    <- ef_sig %>%
    mutate(
      x_end = .data[[x_col]] * scale_fac / max(abs(.data[[x_col]])),
      y_end = .data[[y_col]] * scale_fac / max(abs(.data[[y_col]]))
    )

  p <- ggplot(plot_df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    stat_ellipse(aes(colour = .data[[group_var]], fill = .data[[group_var]]),
                 geom = "polygon", alpha = 0.07, level = 0.95) +
    geom_point(aes(colour = .data[[group_var]]), size = 2.5, alpha = 0.8) +
    # Significant metadata vectors
    geom_segment(data = ef_sig,
                 aes(x = 0, y = 0, xend = x_end, yend = y_end),
                 arrow     = arrow(length = unit(0.2, "cm")),
                 colour    = "#2c3e50", linewidth = 0.8,
                 inherit.aes = FALSE) +
    geom_text(data = ef_sig,
              aes(x = x_end * 1.1, y = y_end * 1.1, label = variable),
              size = 3.5, fontface = "bold", colour = "#2c3e50",
              inherit.aes = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey80") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey80") +
    scale_colour_manual(values = colours, name = group_var) +
    scale_fill_manual(values = colours, guide = "none") +
    labs(
      title    = "Ordination with metadata vectors (envfit)",
      subtitle = paste0(nrow(ef_sig), " significant variables (p < ",
                        p_threshold, ")"),
      x        = x_col, y = y_col
    ) +
    theme_microbiome() + coord_fixed()

  return(list(plot = p, envfit = ef, vectors = ef_df))
}


# =============================================================================
# SECTION 8 — MANTEL TEST
# =============================================================================

#' Run Mantel test to compare two distance matrices.
#'
#' Useful for testing whether microbiome similarity correlates with
#' geographic distance, host genetic distance, or another biological metric.
#'
#' @param dist1       Primary dist object (e.g. microbiome).
#' @param dist2       Secondary dist object (e.g. geography, genetics).
#' @param method      Correlation method: "pearson" or "spearman". Default = "pearson".
#' @param permutations Number of permutations. Default = 999.
#' @return A list: Mantel result and scatter plot.

run_mantel_test <- function(dist1,
                              dist2,
                              method        = "pearson",
                              permutations  = 999) {

  cat("=== Mantel test ===\n")

  set.seed(42)
  mantel_res <- mantel(dist1, dist2,
                       method       = method,
                       permutations = permutations)

  cat("  Mantel r =", round(mantel_res$statistic, 4),
      "| p =", mantel_res$signif, "\n\n")

  # Scatter plot of pairwise distances
  d1_vec <- as.vector(as.matrix(dist1))
  d2_vec <- as.vector(as.matrix(dist2))

  # Remove diagonal (self-comparisons = 0)
  n      <- attr(dist1, "Size")
  off_diag <- !as.logical(diag(n))
  d1_vec <- d1_vec[off_diag]
  d2_vec <- d2_vec[off_diag]

  scatter_df <- data.frame(dist1 = d1_vec, dist2 = d2_vec)

  p <- ggplot(scatter_df, aes(x = dist2, y = dist1)) +
    geom_point(alpha = 0.15, size = 0.8, colour = "#3498db") +
    geom_smooth(method = "lm", se = TRUE,
                colour = "#e74c3c", fill = "#e74c3c", alpha = 0.15) +
    labs(
      title    = "Mantel test — pairwise distance correlation",
      subtitle = paste0("r = ", round(mantel_res$statistic, 4),
                        " | p = ", mantel_res$signif,
                        " | Method: ", method,
                        " | Permutations: ", permutations),
      x        = "Distance matrix 2",
      y        = "Microbiome distance"
    ) +
    theme_microbiome()

  return(list(result = mantel_res, plot = p))
}


# =============================================================================
# SECTION 9 — COMPLETE BETA DIVERSITY WORKFLOW WRAPPER
# =============================================================================

#' Run the complete beta diversity pipeline.
#'
#' @param ps          A filtered phyloseq object (raw counts).
#' @param group_var   Primary metadata variable for group comparisons.
#' @param formula_rhs PERMANOVA formula right-hand side. Default = group_var.
#' @param methods     Distance methods to compute. Default = c("bray","jaccard").
#' @param rare_depth  Rarefaction depth. Default = minimum sample depth.
#' @param env_vars    Numeric metadata variables for envfit. Default = auto-detect.
#' @param output_dir  Directory to save all outputs.
#' @return A named list of all results.

run_beta_diversity <- function(ps,
                                group_var   = NULL,
                                formula_rhs = NULL,
                                methods     = c("bray", "jaccard"),
                                rare_depth  = NULL,
                                env_vars    = NULL,
                                output_dir  = "beta_output") {

  cat("\n", strrep("=", 60), "\n")
  cat("  BETA DIVERSITY PIPELINE\n")
  cat(strrep("=", 60), "\n\n")

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (is.null(formula_rhs)) formula_rhs <- group_var

  results <- list()

  # --- Step 1: Compute distances --------------------------------------------
  results$distances <- compute_distances(
    ps          = ps,
    methods     = methods,
    rarefaction = TRUE,
    rare_depth  = rare_depth
  )

  # Use Bray-Curtis as primary distance (most common in literature)
  primary_dist  <- results$distances[[methods[1]]]
  primary_label <- methods[1]

  # --- Step 2: PCoA ordination ----------------------------------------------
  cat("--- Ordination (PCoA) ---\n")
  results$pcoa <- run_ordination(primary_dist, method = "PCoA", k = 3)

  # PCoA panel (3 axis pairs)
  results$p_pcoa_panel <- plot_pcoa_panel(results$pcoa, ps, group_var)
  ggsave(file.path(output_dir, "01_pcoa_panel.pdf"),
         results$p_pcoa_panel, width = 18, height = 7)

  # Single PCoA (PC1-2) for main figure
  results$p_pcoa_main <- plot_ordination(results$pcoa, ps, group_var,
                                          ellipse = TRUE, spider = FALSE)
  ggsave(file.path(output_dir, "02_pcoa_main.pdf"),
         results$p_pcoa_main, width = 8, height = 7)

  # --- Step 3: NMDS ordination ----------------------------------------------
  cat("--- Ordination (NMDS) ---\n")
  results$nmds  <- run_ordination(primary_dist, method = "NMDS")
  results$p_nmds <- plot_ordination(results$nmds, ps, group_var,
                                     ellipse = TRUE, spider = TRUE)
  ggsave(file.path(output_dir, "03_nmds.pdf"),
         results$p_nmds, width = 8, height = 7)

  # --- Step 4: PERMANOVA ----------------------------------------------------
  cat("--- PERMANOVA ---\n")
  results$permanova <- run_permanova(
    dist_obj     = primary_dist,
    ps           = ps,
    formula_rhs  = formula_rhs,
    permutations = 999
  )

  # Save PERMANOVA result
  capture.output(
    print(results$permanova$permanova),
    file = file.path(output_dir, "permanova_result.txt")
  )
  if (!is.null(results$permanova$pairwise)) {
    write.csv(results$permanova$pairwise,
              file.path(output_dir, "permanova_pairwise.csv"),
              row.names = FALSE)
  }

  # --- Step 5: Betadisper ---------------------------------------------------
  cat("--- Betadisper ---\n")
  results$p_betadisper <- plot_betadisper(primary_dist, ps, group_var)
  ggsave(file.path(output_dir, "04_betadisper.pdf"),
         results$p_betadisper, width = 14, height = 6)

  # --- Step 6: Distance heatmap ---------------------------------------------
  cat("--- Distance heatmap ---\n")
  results$p_dist_heatmap <- plot_distance_heatmap(primary_dist, ps, group_var)
  ggsave(file.path(output_dir, "05_distance_heatmap.pdf"),
         results$p_dist_heatmap, width = 12, height = 10)

  # --- Step 7: Envfit -------------------------------------------------------
  cat("--- Envfit ---\n")
  envfit_res <- tryCatch({
    plot_envfit(results$pcoa, ps, group_var,
                env_vars = env_vars, p_threshold = 0.05)
  }, error = function(e) {
    cat("  Envfit skipped:", e$message, "\n")
    NULL
  })

  if (!is.null(envfit_res)) {
    results$p_envfit   <- envfit_res$plot
    results$envfit_res <- envfit_res$envfit
    ggsave(file.path(output_dir, "06_envfit.pdf"),
           results$p_envfit, width = 9, height = 8)
    write.csv(envfit_res$vectors,
              file.path(output_dir, "envfit_vectors.csv"),
              row.names = FALSE)
  }

  # --- Step 8: Multi-distance comparison plot ------------------------------
  if (length(results$distances) > 1) {
    cat("--- Multi-distance PCoA comparison ---\n")
    p_multi_list <- lapply(names(results$distances), function(m) {
      ord_tmp <- run_ordination(results$distances[[m]], method = "PCoA", k = 2)
      p_tmp   <- plot_ordination(ord_tmp, ps, group_var, ellipse = TRUE)
      p_tmp + labs(title = paste0("PCoA — ", m))
    })

    p_multi <- wrap_plots(p_multi_list, ncol = 2) +
      plot_layout(guides = "collect") +
      plot_annotation(
        title = "Beta diversity — distance method comparison",
        theme = theme(plot.title = element_text(size = 14, face = "bold"))
      )

    results$p_multi <- p_multi
    ggsave(file.path(output_dir, "07_multi_distance_comparison.pdf"),
           results$p_multi, width = 16, height = 7)
  }

  cat("\n", strrep("=", 60), "\n")
  cat("  BETA DIVERSITY PIPELINE COMPLETE\n")
  cat(strrep("=", 60), "\n")
  cat("  Output saved to:", output_dir, "\n")
  cat("  Plots:  up to 7 PDF files\n")
  cat("  Tables: PERMANOVA results, envfit vectors\n\n")

  return(invisible(results))
}


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# Load filtered phyloseq from QC module
# ps <- readRDS("qc_output/phyloseq_qc_filtered.rds")

# --- Option A: Full pipeline -------------------------------------------------
# results <- run_beta_diversity(
#   ps          = ps,
#   group_var   = "disease_status",
#   formula_rhs = "disease_status + age + sex",
#   methods     = c("bray", "jaccard", "wunifrac"),
#   rare_depth  = 10000,
#   env_vars    = c("age", "bmi", "calprotectin"),
#   output_dir  = "results/beta"
# )

# --- Option B: Individual steps ----------------------------------------------
# dists     <- compute_distances(ps, methods = c("bray", "jaccard"))
# pcoa_res  <- run_ordination(dists$bray, method = "PCoA")
# p_pcoa    <- plot_ordination(pcoa_res, ps, group_var = "disease_status")
# perm_res  <- run_permanova(dists$bray, ps, formula_rhs = "disease_status + age")
# perm_res$permanova
# perm_res$pairwise

# --- Option C: Mantel test ---------------------------------------------------
# geo_dist  <- as.dist(geosphere::distm(coords_matrix))
# mantel_res <- run_mantel_test(dists$bray, geo_dist)
# mantel_res$result
# mantel_res$plot
