# =============================================================================
# MICROBIOME LONGITUDINAL ANALYSIS
# =============================================================================
# Description : Comprehensive longitudinal microbiome analysis pipeline
#               covering trajectory modelling, time-series diversity,
#               intervention response, stability analysis, enterotype
#               dynamics, mixed-effects models, and changepoint detection
# Input       : Filtered phyloseq object with time and subject metadata
# Output      : Trajectory plots, stability indices, response curves,
#               mixed-effects model results, changepoint annotations
# Author      : Patricia
# Dependencies: phyloseq, ggplot2, dplyr, tidyr, patchwork, vegan,
#               lme4, lmerTest, scales, tibble, RColorBrewer, stringr,
#               forcats, broom.mixed, ggrepel, zoo
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
  library(ggrepel)
  library(zoo)           # Rolling statistics
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
# SECTION 1 — VALIDATE AND PREPARE LONGITUDINAL DATA
# =============================================================================

#' Validate and prepare a phyloseq object for longitudinal analysis.
#'
#' Checks that required time and subject variables are present,
#' computes per-sample metadata summaries, and returns a structured
#' data object for downstream functions.
#'
#' @param ps          A phyloseq object (raw counts).
#' @param time_var    Metadata column for time points. Must be numeric or ordered factor.
#' @param subject_var Metadata column for subject/individual identifiers.
#' @param group_var   Optional metadata column for treatment/disease group.
#' @param rank        Taxonomic rank to agglomerate to. Default = "Genus".
#' @param transform   Abundance transformation: "clr", "relative". Default = "clr".
#' @param min_timepoints Minimum time points per subject. Default = 2.
#' @return A named list with phyloseq, metadata, diversity, and feature matrix.

prepare_longitudinal_data <- function(ps,
                                       time_var       = "timepoint",
                                       subject_var    = "subject_id",
                                       group_var      = NULL,
                                       rank           = "Genus",
                                       transform      = "clr",
                                       min_timepoints = 2) {

  cat("=== Preparing longitudinal data ===\n")

  # --- Validate metadata ---------------------------------------------------
  meta_df <- data.frame(sample_data(ps))

  required_vars <- c(time_var, subject_var)
  missing_vars  <- setdiff(required_vars, colnames(meta_df))
  if (length(missing_vars) > 0) {
    stop("Missing required metadata columns: ", paste(missing_vars, collapse = ", "))
  }

  # Ensure time is numeric
  meta_df[[time_var]] <- as.numeric(as.character(meta_df[[time_var]]))

  # Filter subjects with insufficient time points
  subject_counts <- table(meta_df[[subject_var]])
  valid_subjects <- names(subject_counts[subject_counts >= min_timepoints])
  keep_samples   <- meta_df[[subject_var]] %in% valid_subjects

  if (sum(keep_samples) < nsamples(ps)) {
    n_removed <- sum(!keep_samples)
    cat("  Removed", n_removed, "samples from subjects with <",
        min_timepoints, "time points\n")
    ps <- prune_samples(keep_samples, ps)
    meta_df <- meta_df[keep_samples, , drop = FALSE]
  }

  # Summarise study design
  n_subjects  <- n_distinct(meta_df[[subject_var]])
  time_points <- sort(unique(meta_df[[time_var]]))
  n_times     <- length(time_points)

  cat("  Subjects:", n_subjects, "\n")
  cat("  Time points:", paste(time_points, collapse = ", "), "\n")
  cat("  Total samples:", nsamples(ps), "\n")

  if (!is.null(group_var) && group_var %in% colnames(meta_df)) {
    cat("  Groups:", paste(unique(meta_df[[group_var]]), collapse = ", "), "\n")
  }

  # --- Agglomerate and transform -------------------------------------------
  ps_agg    <- tax_glom(ps, taxrank = rank, NArm = FALSE)
  new_names <- as.character(tax_table(ps_agg)[, rank])
  new_names[is.na(new_names)] <- paste0("Unknown_", seq_len(sum(is.na(new_names))))
  dup_idx   <- duplicated(new_names)
  new_names <- make.unique(new_names, sep = "_dup")
  taxa_names(ps_agg) <- new_names

  otu_mat <- as.matrix(otu_table(ps_agg))
  if (!taxa_are_rows(ps_agg)) otu_mat <- t(otu_mat)

  if (transform == "clr") {
    feat_mat <- t(apply(otu_mat + 0.5, 2, function(x) log(x) - mean(log(x))))
  } else {
    feat_mat <- t(apply(otu_mat, 2, function(x) x / sum(x)))
  }

  # --- Alpha diversity at each time point -----------------------------------
  otu_counts <- as.matrix(otu_table(ps))
  if (!taxa_are_rows(ps)) otu_counts <- t(otu_counts)
  otu_t      <- t(otu_counts)

  div_df <- data.frame(
    sample      = rownames(otu_t),
    observed    = rowSums(otu_t > 0),
    shannon     = round(vegan::diversity(otu_t, index = "shannon"), 4),
    simpson     = round(vegan::diversity(otu_t, index = "simpson"), 4),
    pielou      = NA,
    stringsAsFactors = FALSE
  )
  div_df$pielou <- ifelse(div_df$observed > 1,
                           round(div_df$shannon / log(div_df$observed), 4), 0)

  # Add metadata
  div_full <- div_df %>%
    left_join(meta_df %>% rownames_to_column("sample"), by = "sample")

  cat("\n")

  return(list(
    ps           = ps_agg,
    ps_raw       = ps,
    metadata     = meta_df,
    diversity    = div_full,
    features     = feat_mat,
    time_var     = time_var,
    subject_var  = subject_var,
    group_var    = group_var,
    time_points  = time_points,
    n_subjects   = n_subjects,
    rank         = rank
  ))
}


# =============================================================================
# SECTION 2 — LONGITUDINAL ALPHA DIVERSITY TRAJECTORIES
# =============================================================================

#' Plot alpha diversity trajectories over time per subject and group.
#'
#' @param long_data   Output from prepare_longitudinal_data().
#' @param metrics     Diversity metrics to plot. Default = c("shannon","observed").
#' @param test        Whether to run Wilcoxon signed-rank test. Default = TRUE.
#' @return A list: patchwork plot and statistical test results.

plot_diversity_trajectories <- function(long_data,
                                         metrics = c("shannon", "observed",
                                                     "pielou"),
                                         test    = TRUE) {

  cat("=== Alpha diversity trajectories ===\n")

  div_df    <- long_data$diversity
  time_var  <- long_data$time_var
  subj_var  <- long_data$subject_var
  group_var <- long_data$group_var

  metrics   <- intersect(metrics, colnames(div_df))
  n_groups  <- if (!is.null(group_var)) n_distinct(div_df[[group_var]]) else 1
  colours   <- brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)]

  metric_labels <- c(
    shannon  = "Shannon entropy (H')",
    observed = "Observed richness",
    pielou   = "Pielou's evenness (J')",
    simpson  = "Simpson diversity (1-D)"
  )

  stat_results <- list()
  plots        <- list()

  for (metric in metrics) {

    # Group summary for ribbon
    if (!is.null(group_var)) {
      grp_summary <- div_df %>%
        group_by(.data[[time_var]], .data[[group_var]]) %>%
        summarise(
          mean_val = mean(.data[[metric]], na.rm = TRUE),
          se_val   = sd(.data[[metric]], na.rm = TRUE) / sqrt(n()),
          .groups  = "drop"
        )
    } else {
      grp_summary <- div_df %>%
        group_by(.data[[time_var]]) %>%
        summarise(
          mean_val = mean(.data[[metric]], na.rm = TRUE),
          se_val   = sd(.data[[metric]], na.rm = TRUE) / sqrt(n()),
          .groups  = "drop"
        ) %>%
        mutate(group_dummy = "All")
    }

    p <- ggplot(div_df, aes(x = .data[[time_var]], y = .data[[metric]])) +
      # Individual trajectories
      geom_line(
        aes(group = .data[[subj_var]],
            colour = if (!is.null(group_var)) .data[[group_var]] else "All"),
        alpha = 0.2, linewidth = 0.4
      ) +
      geom_point(
        aes(colour = if (!is.null(group_var)) .data[[group_var]] else "All"),
        alpha = 0.3, size = 1.2
      )

    # Group mean + SE ribbon
    if (!is.null(group_var)) {
      p <- p +
        geom_ribbon(
          data = grp_summary,
          aes(x    = .data[[time_var]],
              ymin = mean_val - se_val,
              ymax = mean_val + se_val,
              fill = .data[[group_var]]),
          alpha = 0.2, inherit.aes = FALSE
        ) +
        geom_line(
          data      = grp_summary,
          aes(x = .data[[time_var]], y = mean_val,
              colour = .data[[group_var]], group = .data[[group_var]]),
          linewidth = 1.3, inherit.aes = FALSE
        ) +
        geom_point(
          data = grp_summary,
          aes(x = .data[[time_var]], y = mean_val,
              colour = .data[[group_var]]),
          size = 3, inherit.aes = FALSE
        ) +
        scale_colour_manual(values = colours, name = group_var) +
        scale_fill_manual(values = colours, guide = "none")
    } else {
      p <- p +
        geom_ribbon(
          data = grp_summary,
          aes(x = .data[[time_var]],
              ymin = mean_val - se_val,
              ymax = mean_val + se_val),
          fill = "#3498db", alpha = 0.2, inherit.aes = FALSE
        ) +
        geom_line(
          data = grp_summary,
          aes(x = .data[[time_var]], y = mean_val),
          colour = "#3498db", linewidth = 1.3, inherit.aes = FALSE
        ) +
        scale_colour_manual(values = "#3498db", guide = "none")
    }

    # Statistical annotation
    if (test && !is.null(group_var) && n_groups == 2) {
      group_levels <- unique(div_df[[group_var]])
      time_pts     <- sort(unique(div_df[[time_var]]))

      test_res_list <- lapply(time_pts, function(tp) {
        g1 <- div_df[[metric]][div_df[[time_var]] == tp &
                                 div_df[[group_var]] == group_levels[1]]
        g2 <- div_df[[metric]][div_df[[time_var]] == tp &
                                 div_df[[group_var]] == group_levels[2]]
        if (length(g1) < 3 || length(g2) < 3) return(NULL)
        wt <- tryCatch(wilcox.test(g1, g2, exact = FALSE), error = function(e) NULL)
        if (is.null(wt)) return(NULL)
        data.frame(timepoint = tp, p_value = wt$p.value, stringsAsFactors = FALSE)
      })

      test_res <- bind_rows(test_res_list) %>%
        mutate(
          p_adj = p.adjust(p_value, method = "BH"),
          sig   = case_when(
            p_adj < 0.001 ~ "***",
            p_adj < 0.01  ~ "**",
            p_adj < 0.05  ~ "*",
            TRUE           ~ ""
          )
        )
      stat_results[[metric]] <- test_res

      sig_df <- test_res %>% filter(sig != "")
      if (nrow(sig_df) > 0) {
        y_max <- max(div_df[[metric]], na.rm = TRUE)
        p <- p +
          geom_text(
            data = sig_df,
            aes(x = timepoint, y = y_max * 1.05, label = sig),
            inherit.aes = FALSE, size = 4, colour = "black", fontface = "bold"
          )
      }
    }

    p <- p +
      scale_x_continuous(breaks = long_data$time_points) +
      labs(
        title = metric_labels[[metric]] %||% metric,
        x     = time_var,
        y     = metric_labels[[metric]] %||% metric
      ) +
      theme_microbiome()

    plots[[metric]] <- p
  }

  n_cols   <- min(length(plots), 3)
  combined <- wrap_plots(plots, ncol = n_cols) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = "Alpha diversity over time",
      subtitle = paste0(long_data$n_subjects, " subjects | ",
                        length(long_data$time_points), " time points"),
      theme    = theme(plot.title    = element_text(size = 14, face = "bold"),
                       plot.subtitle = element_text(size = 10, colour = "grey40"))
    )

  cat("  Metrics plotted:", paste(metrics, collapse = ", "), "\n\n")

  return(list(plot = combined, stats = stat_results))
}


# =============================================================================
# SECTION 3 — BETA DIVERSITY TRAJECTORIES
# =============================================================================

#' Track beta diversity trajectories over time using PCoA ordination.
#'
#' Connects sequential samples from the same subject with arrows
#' to show directional movement through community space.
#'
#' @param long_data     Output from prepare_longitudinal_data().
#' @param distance      Distance metric: "bray", "jaccard", "aitchison".
#' @param rarefaction   Whether to rarefy. Default = TRUE.
#' @param rare_depth    Rarefaction depth. Default = minimum sample depth.
#' @param show_centroid Whether to show group centroids. Default = TRUE.
#' @return A list: ordination plot with trajectories and distance-to-baseline.

plot_beta_trajectories <- function(long_data,
                                    distance    = "bray",
                                    rarefaction = TRUE,
                                    rare_depth  = NULL,
                                    show_centroid = TRUE) {

  cat("=== Beta diversity trajectories ===\n")

  ps        <- long_data$ps_raw
  time_var  <- long_data$time_var
  subj_var  <- long_data$subject_var
  group_var <- long_data$group_var
  meta_df   <- long_data$metadata

  # Rarefy
  if (rarefaction) {
    if (is.null(rare_depth)) rare_depth <- min(sample_sums(ps))
    set.seed(42)
    ps <- rarefy_even_depth(ps, sample.size = rare_depth,
                             rngseed = 42, replace = FALSE, verbose = FALSE)
  }

  # Distance matrix
  otu_mat <- as.matrix(otu_table(ps))
  if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  otu_t   <- matrix(t(otu_mat), nrow=ncol(otu_mat), ncol=nrow(otu_mat),
                   dimnames=list(colnames(otu_mat), rownames(otu_mat)))

  if (distance == "aitchison") {
    clr_mat  <- t(apply(otu_mat + 0.5, 2, function(x) log(x) - mean(log(x))))
    dist_obj <- dist(clr_mat)
  } else {
    dist_obj <- vegdist(otu_t, method = distance)
  }

  # PCoA
  pcoa     <- cmdscale(dist_obj, k = 2, eig = TRUE)
  eig_pos  <- pcoa$eig[pcoa$eig > 0]
  var_exp  <- round(eig_pos / sum(eig_pos) * 100, 1)

  pcoa_df  <- as.data.frame(pcoa$points) %>%
    setNames(c("PC1", "PC2")) %>%
    rownames_to_column("sample") %>%
    left_join(meta_df %>% rownames_to_column("sample"), by = "sample") %>%
    arrange(.data[[subj_var]], .data[[time_var]])

  n_groups <- if (!is.null(group_var)) n_distinct(pcoa_df[[group_var]]) else 1
  colours  <- brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)]

  # Build trajectory segments (arrows between consecutive time points)
  seg_df <- pcoa_df %>%
    group_by(.data[[subj_var]]) %>%
    arrange(.data[[time_var]]) %>%
    mutate(
      PC1_end = lead(PC1),
      PC2_end = lead(PC2),
      time_end = lead(.data[[time_var]])
    ) %>%
    filter(!is.na(PC1_end)) %>%
    ungroup()

  # Time colour scale
  time_vals <- sort(unique(pcoa_df[[time_var]]))
  n_times   <- length(time_vals)
  time_pal  <- colorRampPalette(c("#3498db", "#e74c3c"))(n_times)
  names(time_pal) <- as.character(time_vals)

  p <- ggplot(pcoa_df, aes(x = PC1, y = PC2))

  # Trajectory arrows
  if (!is.null(group_var)) {
    p <- p +
      geom_segment(
        data = seg_df,
        aes(x = PC1, y = PC2, xend = PC1_end, yend = PC2_end,
            colour = .data[[group_var]]),
        arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
        alpha = 0.35, linewidth = 0.5
      )
  } else {
    p <- p +
      geom_segment(
        data = seg_df,
        aes(x = PC1, y = PC2, xend = PC1_end, yend = PC2_end),
        colour = "grey60",
        arrow  = arrow(length = unit(0.12, "cm"), type = "closed"),
        alpha  = 0.3, linewidth = 0.5
      )
  }

  # Points coloured by time
  p <- p +
    geom_point(
      aes(fill  = as.character(.data[[time_var]]),
          shape = if (!is.null(group_var)) .data[[group_var]] else NULL),
      size = 3, alpha = 0.85, colour = "white", stroke = 0.4
    ) +
    scale_fill_manual(values = time_pal, name = paste0(time_var, "\n(colour)")) +
    scale_shape_manual(values = c(21, 22, 23, 24, 25)[seq_len(n_groups)],
                       name = group_var)

  # Group centroids per time point
  if (show_centroid && !is.null(group_var)) {
    centroids <- pcoa_df %>%
      group_by(.data[[group_var]], .data[[time_var]]) %>%
      summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")

    centroid_segs <- centroids %>%
      group_by(.data[[group_var]]) %>%
      arrange(.data[[time_var]]) %>%
      mutate(PC1_end = lead(PC1), PC2_end = lead(PC2)) %>%
      filter(!is.na(PC1_end))

    p <- p +
      geom_path(data = centroids,
                aes(group = .data[[group_var]], colour = .data[[group_var]]),
                linewidth = 1.5, alpha = 0.8) +
      geom_point(data = centroids,
                 aes(colour = .data[[group_var]]),
                 size = 5, shape = 8, stroke = 1.5) +
      scale_colour_manual(values = colours, name = group_var)
  }

  p <- p +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey80") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey80") +
    labs(
      title    = paste0("Community trajectory - ", distance, " PCoA"),
      subtitle = paste0("Arrows connect consecutive time points per subject"),
      x        = paste0("PC1 (", var_exp[1], "%)"),
      y        = paste0("PC2 (", var_exp[2], "%)")
    ) +
    theme_microbiome() +
    coord_fixed()

  # --- Distance to baseline ------------------------------------------------
  # Calculate Bray-Curtis distance of each sample from that subject's baseline
  dist_mat  <- as.matrix(dist_obj)
  meta_ord  <- meta_df[rownames(dist_mat), , drop = FALSE]

  baseline_dist <- lapply(unique(meta_ord[[subj_var]]), function(subj) {
    subj_samples <- rownames(meta_ord)[meta_ord[[subj_var]] == subj]
    time_vals_s  <- meta_ord[[time_var]][meta_ord[[subj_var]] == subj]
    baseline_s   <- subj_samples[which.min(time_vals_s)]

    lapply(subj_samples, function(s) {
      data.frame(
        sample   = s,
        subject  = subj,
        time     = meta_ord[[time_var]][rownames(meta_ord) == s],
        dist_to_baseline = dist_mat[s, baseline_s],
        stringsAsFactors = FALSE
      )
    }) %>% bind_rows()
  }) %>% bind_rows()

  baseline_dist <- left_join(baseline_dist,
                              meta_df %>% rownames_to_column("sample") %>%
                                select(sample, any_of(group_var)),
                              by = "sample")

  grp_summary_bd <- if (!is.null(group_var)) {
    baseline_dist %>%
      group_by(time, .data[[group_var]]) %>%
      summarise(
        mean_dist = mean(dist_to_baseline),
        se_dist   = sd(dist_to_baseline) / sqrt(n()),
        .groups   = "drop"
      )
  } else {
    baseline_dist %>%
      group_by(time) %>%
      summarise(
        mean_dist = mean(dist_to_baseline),
        se_dist   = sd(dist_to_baseline) / sqrt(n()),
        .groups   = "drop"
      ) %>%
      mutate(group_dummy = "All")
  }

  p_dist <- ggplot(baseline_dist, aes(x = time, y = dist_to_baseline)) +
    geom_line(
      aes(group  = subject,
          colour = if (!is.null(group_var)) .data[[group_var]] else "All"),
      alpha = 0.2, linewidth = 0.4
    ) +
    {
      if (!is.null(group_var)) {
        list(
          geom_ribbon(data = grp_summary_bd,
                      aes(x = time,
                          ymin  = mean_dist - se_dist,
                          ymax  = mean_dist + se_dist,
                          fill  = .data[[group_var]]),
                      alpha = 0.2, inherit.aes = FALSE),
          geom_line(data = grp_summary_bd,
                    aes(x = time, y = mean_dist,
                        colour = .data[[group_var]],
                        group  = .data[[group_var]]),
                    linewidth = 1.3, inherit.aes = FALSE),
          scale_colour_manual(values = colours, name = group_var),
          scale_fill_manual(values = colours, guide = "none")
        )
      } else {
        list(
          geom_ribbon(data = grp_summary_bd,
                      aes(x = time,
                          ymin  = mean_dist - se_dist,
                          ymax  = mean_dist + se_dist),
                      fill = "#3498db", alpha = 0.2, inherit.aes = FALSE),
          geom_line(data = grp_summary_bd,
                    aes(x = time, y = mean_dist),
                    colour = "#3498db", linewidth = 1.3, inherit.aes = FALSE)
        )
      }
    } +
    scale_x_continuous(breaks = long_data$time_points) +
    labs(
      title    = "Distance to baseline over time",
      subtitle = paste0(distance, " dissimilarity from T0 per subject"),
      x        = time_var,
      y        = paste0(str_to_title(distance), " distance from baseline")
    ) +
    theme_microbiome()

  combined <- p / p_dist +
    plot_layout(heights = c(2, 1), guides = "collect") +
    plot_annotation(title = "Beta Diversity Trajectories")

  cat("  PCoA variance explained: PC1 =", var_exp[1],
      "%, PC2 =", var_exp[2], "%\n\n")

  return(list(
    plot           = combined,
    p_pcoa         = p,
    p_dist_baseline = p_dist,
    pcoa_coords    = pcoa_df,
    dist_baseline  = baseline_dist
  ))
}


# =============================================================================
# SECTION 4 — INTRA-INDIVIDUAL STABILITY ANALYSIS
# =============================================================================

#' Measure within-subject microbiome stability over time.
#'
#' Computes consecutive-timepoint Bray-Curtis dissimilarity per subject.
#' Lower values indicate a more stable microbiome. Compares stability
#' across groups and identifies unstable subjects.
#'
#' @param long_data   Output from prepare_longitudinal_data().
#' @param distance    Distance metric. Default = "bray".
#' @return A list: stability table, violin plot, and trajectory plot.

analyse_stability <- function(long_data, distance = "bray") {

  cat("=== Intra-individual stability analysis ===\n")

  ps        <- long_data$ps_raw
  time_var  <- long_data$time_var
  subj_var  <- long_data$subject_var
  group_var <- long_data$group_var
  meta_df   <- long_data$metadata

  otu_mat   <- as.matrix(otu_table(ps))
  if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  otu_t     <- t(otu_mat)

  dist_mat  <- as.matrix(vegdist(otu_t, method = distance))

  # Consecutive-timepoint dissimilarity per subject
  subjects  <- unique(meta_df[[subj_var]])

  stability_list <- lapply(subjects, function(subj) {
    subj_rows <- which(meta_df[[subj_var]] == subj)
    subj_meta <- meta_df[subj_rows, , drop = FALSE] %>%
      rownames_to_column("sample") %>%
      arrange(.data[[time_var]])

    if (nrow(subj_meta) < 2) return(NULL)

    pairs <- lapply(seq_len(nrow(subj_meta) - 1), function(i) {
      s1     <- subj_meta$sample[i]
      s2     <- subj_meta$sample[i + 1]
      t1     <- subj_meta[[time_var]][i]
      t2     <- subj_meta[[time_var]][i + 1]
      d      <- dist_mat[s1, s2]
      data.frame(
        subject   = subj,
        from_time = t1,
        to_time   = t2,
        interval  = t2 - t1,
        distance  = round(d, 4),
        stringsAsFactors = FALSE
      )
    })
    bind_rows(pairs)
  })

  stability_df <- bind_rows(stability_list)

  # Add group variable
  subj_group <- meta_df %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::group_by(.data[[subj_var]]) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      subject = dplyr::all_of(subj_var),
      dplyr::any_of(group_var)
    )

  if (!is.null(group_var)) {
    stability_df <- dplyr::left_join(
      stability_df,
      subj_group,
      by = "subject"
    )
  }

  # Per-subject mean stability
  subj_stability <- stability_df %>%
    group_by(subject) %>%
    summarise(
      mean_distance = round(mean(distance), 4),
      sd_distance   = round(sd(distance), 4),
      n_intervals   = n(),
      .groups = "drop"
    ) %>%
    arrange(mean_distance)

  # Add group
  if (!is.null(group_var)) {
    subj_stability <- dplyr::left_join(
      subj_stability,
      subj_group,
      by = "subject"
    )
  }

  cat("  Overall median consecutive distance:", round(median(stability_df$distance), 4), "\n")

  if (!is.null(group_var)) {
    group_stab <- stability_df %>%
      group_by(.data[[group_var]]) %>%
      summarise(
        median_dist = round(median(distance), 4),
        mean_dist   = round(mean(distance), 4),
        .groups = "drop"
      )
    cat("  Group stability:\n")
    print(group_stab)
  }
  cat("\n")

  n_groups <- if (!is.null(group_var)) n_distinct(stability_df[[group_var]]) else 1
  colours  <- brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)]

  # --- Plot 1: Violin plot of stability per group --------------------------
  p_violin <- ggplot(
    stability_df,
    aes(
      x    = if (!is.null(group_var)) .data[[group_var]] else "All",
      y    = distance,
      fill = if (!is.null(group_var)) .data[[group_var]] else "All"
    )
  ) +
    geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1.2) +
    scale_fill_manual(values = colours, guide = "none") +
    labs(
      title    = "Microbiome stability - consecutive distance",
      subtitle = paste0(distance, " dissimilarity between adjacent time points"),
      x        = if (!is.null(group_var)) group_var else "",
      y        = paste0("Consecutive ", distance, " distance")
    ) +
    theme_microbiome()

  # --- Plot 2: Stability per subject as ranked bar -------------------------
  p_rank <- ggplot(
    subj_stability,
    aes(
      x    = reorder(subject, mean_distance),
      y    = mean_distance,
      fill = if (!is.null(group_var)) .data[[group_var]] else "All"
    )
  ) +
    geom_col(width = 0.8, alpha = 0.85) +
    geom_errorbar(
      aes(ymin = mean_distance - sd_distance,
          ymax = mean_distance + sd_distance),
      width = 0.3, linewidth = 0.5
    ) +
    scale_fill_manual(values = colours,
                      name   = if (!is.null(group_var)) group_var else NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
      title    = "Per-subject microbiome stability",
      subtitle = "Mean +/- SD of consecutive-timepoint distances",
      x        = "Subject (ranked by stability)",
      y        = "Mean distance"
    ) +
    theme_microbiome() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  combined <- p_violin | p_rank

  return(list(
    plot             = combined,
    p_violin         = p_violin,
    p_rank           = p_rank,
    stability_df     = stability_df,
    subject_stability = subj_stability
  ))
}


# =============================================================================
# SECTION 5 — TAXONOMIC COMPOSITION OVER TIME
# =============================================================================

#' Plot how taxonomic composition changes over time.
#'
#' @param long_data   Output from prepare_longitudinal_data().
#' @param rank        Taxonomic rank for composition. Default = "Phylum".
#' @param top_n       Top taxa to display. Default = 10.
#' @param facet_by    Whether to facet by group or subject. Default = "group".
#' @return A list: stacked bar plot and relative change plot.

plot_composition_over_time <- function(long_data,
                                        rank     = "Phylum",
                                        top_n    = 10,
                                        facet_by = "group") {

  cat("=== Composition over time ===\n")

  ps        <- long_data$ps_raw
  time_var  <- long_data$time_var
  group_var <- long_data$group_var

  # Agglomerate to chosen rank
  ps_phy    <- tax_glom(ps, taxrank = rank, NArm = FALSE)
  ps_rel    <- transform_sample_counts(ps_phy, function(x) x / sum(x))
  taxa_names(ps_rel) <- as.character(tax_table(ps_rel)[, rank])

  # Top taxa
  mean_abund <- rowMeans(as.matrix(otu_table(ps_rel)))
  top_taxa   <- names(sort(mean_abund, decreasing = TRUE))[seq_len(min(top_n, length(mean_abund)))]

  ps_top    <- prune_taxa(top_taxa, ps_rel)
  ps_melt   <- psmelt(ps_top) %>%
    mutate(
      taxon      = OTU,
      time       = as.numeric(as.character(.data[[time_var]])),
      taxon      = factor(taxon, levels = rev(top_taxa))
    )

  # Group-level mean per time point
  grp_vars <- c("time", "taxon",
                if (!is.null(group_var)) group_var else NULL)

  comp_summary <- ps_melt %>%
    group_by(across(all_of(grp_vars))) %>%
    summarise(
      mean_abund = mean(Abundance),
      se_abund   = sd(Abundance) / sqrt(n()),
      .groups    = "drop"
    )

  # Colour palette
  n_taxa   <- length(top_taxa)
  pal      <- c(brewer.pal(min(8, n_taxa), "Set2"),
                colorRampPalette(brewer.pal(8, "Set2"))(max(0, n_taxa - 8)))
  names(pal) <- rev(top_taxa)

  # --- Plot 1: Stacked area chart ------------------------------------------
  p_area <- ggplot(
    comp_summary,
    aes(x = time, y = mean_abund, fill = taxon)
  ) +
    geom_area(position = "stack", alpha = 0.9, colour = "white",
              linewidth = 0.15) +
    scale_fill_manual(values = pal, name = rank) +
    scale_x_continuous(breaks = long_data$time_points) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = c(0, 0)) +
    labs(
      title    = paste0(rank, " composition over time"),
      subtitle = paste0("Mean relative abundance | top ", top_n, " taxa"),
      x        = time_var,
      y        = "Mean relative abundance"
    ) +
    theme_microbiome()

  if (!is.null(group_var) && facet_by == "group") {
    p_area <- p_area + facet_wrap(as.formula(paste("~", group_var)))
  }

  # --- Plot 2: Relative change from baseline per taxon ---------------------
  if (!is.null(group_var)) {
    baseline_time <- min(long_data$time_points)
    baseline_comp <- comp_summary %>%
      filter(time == baseline_time) %>%
      select(taxon, all_of(group_var), baseline_abund = mean_abund)

    rel_change <- comp_summary %>%
      left_join(baseline_comp, by = c("taxon", group_var)) %>%
      mutate(
        rel_change = (mean_abund - baseline_abund) / (baseline_abund + 0.001)
      )

    p_change <- ggplot(
      rel_change %>% filter(time != baseline_time),
      aes(x = time, y = rel_change,
          colour = taxon, group = taxon)
    ) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_line(linewidth = 0.9, alpha = 0.8) +
      geom_point(size = 2.5, alpha = 0.8) +
      scale_colour_manual(values = pal, name = rank) +
      scale_x_continuous(breaks = long_data$time_points) +
      scale_y_continuous(labels = label_percent(accuracy = 1)) +
      facet_wrap(as.formula(paste("~", group_var))) +
      labs(
        title    = "Relative change in taxon abundance from baseline",
        subtitle = paste0("Baseline = T", baseline_time),
        x        = time_var,
        y        = "Relative change from baseline"
      ) +
      theme_microbiome()

    combined <- p_area / p_change +
      plot_layout(heights = c(1, 1), guides = "collect")
  } else {
    combined <- p_area
  }

  cat("  Taxa plotted:", top_n, "| Time points:", length(long_data$time_points), "\n\n")

  return(list(
    plot     = combined,
    p_area   = p_area,
    data     = comp_summary
  ))
}


# =============================================================================
# SECTION 6 — MIXED-EFFECTS MODELLING
# =============================================================================

#' Run linear mixed-effects models for all taxa over time.
#'
#' Tests whether each taxon changes significantly over time,
#' controlling for covariates and accounting for repeated measures
#' within subjects.
#'
#' @param long_data   Output from prepare_longitudinal_data().
#' @param fixed_vars  Additional fixed-effect covariates. Default = NULL.
#' @param interaction Whether to test time × group interaction. Default = TRUE.
#' @param alpha       Significance threshold. Default = 0.05.
#' @param top_n       Top significant taxa to visualise. Default = 20.
#' @return A list: results table and forest plot.

run_lme_over_time <- function(long_data,
                               fixed_vars  = NULL,
                               interaction = TRUE,
                               alpha       = 0.05,
                               top_n       = 20) {

  cat("=== Mixed-effects model over time ===\n")

  if (!pkg_available("lme4") || !pkg_available("lmerTest")) {
    cat("  Install lme4 and lmerTest for mixed-effects models.\n\n")
    return(NULL)
  }

  library(lme4)
  library(lmerTest)

  feat_mat  <- long_data$features
  meta_df   <- long_data$metadata %>% rownames_to_column("sample")
  time_var  <- long_data$time_var
  subj_var  <- long_data$subject_var
  group_var <- long_data$group_var

  # Align
  shared   <- intersect(rownames(feat_mat), meta_df$sample)
  feat_sub <- feat_mat[shared, , drop = FALSE]
  meta_sub <- meta_df[match(shared, meta_df$sample), ]

  # Build formula
  cov_str  <- if (!is.null(fixed_vars))
    paste("+", paste(fixed_vars, collapse = " + ")) else ""

  if (!is.null(group_var) && interaction) {
    formula_str <- paste0("abundance ~ ", time_var, " * ", group_var,
                           cov_str, " + (1|", subj_var, ")")
  } else if (!is.null(group_var)) {
    formula_str <- paste0("abundance ~ ", time_var, " + ", group_var,
                           cov_str, " + (1|", subj_var, ")")
  } else {
    formula_str <- paste0("abundance ~ ", time_var,
                           cov_str, " + (1|", subj_var, ")")
  }

  cat("  Formula:", formula_str, "\n")
  cat("  Testing", ncol(feat_sub), "taxa...\n")

  # Filter to prevalent taxa
  prev     <- colMeans(feat_sub != 0)
  feat_sub <- feat_sub[, prev >= 0.20, drop = FALSE]
  cat("  After prevalence filter:", ncol(feat_sub), "taxa\n")

  results <- lapply(colnames(feat_sub), function(tx) {
    model_df <- cbind(
      data.frame(abundance = feat_sub[, tx]),
      meta_sub
    )

    tryCatch({
      fit  <- lmer(as.formula(formula_str), data = model_df,
                   REML = FALSE,
                   control = lmerControl(optimizer = "bobyqa"))
      coef_tbl <- as.data.frame(summary(fit)$coefficients)

      # Extract time effect
      time_row <- grep(paste0("^", time_var), rownames(coef_tbl), value = TRUE)

      lapply(time_row, function(rn) {
        data.frame(
          taxon     = tx,
          term      = rn,
          estimate  = round(coef_tbl[rn, "Estimate"], 4),
          se        = round(coef_tbl[rn, "Std. Error"], 4),
          t_value   = round(coef_tbl[rn, "t value"], 4),
          p_value   = coef_tbl[rn, "Pr(>|t|)"],
          stringsAsFactors = FALSE
        )
      }) %>% bind_rows()
    }, error = function(e) NULL)
  }) %>% bind_rows()

  if (is.null(results) || nrow(results) == 0) {
    cat("  No results returned.\n\n")
    return(NULL)
  }

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

  # --- Forest plot ---------------------------------------------------------
  forest_df <- results %>%
    filter(significant) %>%
    arrange(desc(estimate)) %>%
    slice_head(n = top_n) %>%
    mutate(
      taxon     = fct_reorder(taxon, estimate),
      direction = ifelse(estimate > 0, "Increasing", "Decreasing")
    )

  colour_map <- c("Increasing" = "#e74c3c", "Decreasing" = "#3498db")

  p_forest <- if (nrow(forest_df) > 0) {
    ggplot(forest_df,
           aes(x = estimate, y = taxon,
               colour = direction)) +
      geom_vline(xintercept = 0, linetype = "dashed",
                 colour = "grey50", linewidth = 0.8) +
      geom_errorbarh(
        aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se),
        height = 0.3, linewidth = 0.6
      ) +
      geom_point(size = 3) +
      geom_text(aes(label = sig_label), nudge_y = 0.35,
                size = 3, colour = "black", fontface = "bold") +
      scale_colour_manual(values = colour_map, name = "Trend") +
      labs(
        title    = paste0("Taxa changing over time - mixed-effects model"),
        subtitle = paste0(n_sig, " significant taxa (q < ", alpha, ") | ",
                          "Top ", nrow(forest_df), " shown"),
        x        = paste0("Fixed effect estimate (", time_var, ")"),
        y        = NULL
      ) +
      theme_microbiome() +
      theme(axis.text.y = element_text(face = "italic", size = 8))
  } else {
    ggplot() + annotate("text", x = 0.5, y = 0.5,
                        label = "No significant taxa", size = 6) + theme_void()
  }

  return(list(
    results  = results,
    plot     = p_forest,
    n_sig    = n_sig
  ))
}


# =============================================================================
# SECTION 7 — INTERVENTION RESPONSE ANALYSIS
# =============================================================================

#' Analyse microbiome response to an intervention at a defined time point.
#'
#' Compares pre- vs post-intervention diversity and composition,
#' and identifies taxa that respond to the intervention.
#'
#' @param long_data        Output from prepare_longitudinal_data().
#' @param pre_timepoints   Numeric vector of pre-intervention time points.
#' @param post_timepoints  Numeric vector of post-intervention time points.
#' @param intervention_label Name of the intervention (for plot titles).
#' @return A list: response plots and responder/non-responder classification.

analyse_intervention_response <- function(long_data,
                                           pre_timepoints,
                                           post_timepoints,
                                           intervention_label = "Intervention") {

  cat("=== Intervention response analysis ===\n")
  cat("  Pre time points: ", paste(pre_timepoints, collapse = ", "), "\n")
  cat("  Post time points:", paste(post_timepoints, collapse = ", "), "\n")

  div_df    <- long_data$diversity
  time_var  <- long_data$time_var
  subj_var  <- long_data$subject_var
  group_var <- long_data$group_var

  # Label pre/post
  div_df <- div_df %>%
    mutate(period = case_when(
      .data[[time_var]] %in% pre_timepoints  ~ "Pre",
      .data[[time_var]] %in% post_timepoints ~ "Post",
      TRUE                                    ~ "Other"
    )) %>%
    filter(period != "Other") %>%
    mutate(period = factor(period, levels = c("Pre", "Post")))

  n_groups <- if (!is.null(group_var)) n_distinct(div_df[[group_var]]) else 1
  colours  <- brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)]

  # --- Per-subject pre vs post paired plot ---------------------------------
  paired_df <- div_df %>%
    group_by(.data[[subj_var]], period) %>%
    summarise(shannon = mean(shannon), observed = mean(observed), .groups = "drop")

  if (!is.null(group_var)) {
    grp_df <- div_df %>%
      select(.data[[subj_var]], .data[[group_var]]) %>%
      distinct()
    paired_df <- left_join(paired_df, grp_df, by = subj_var)
  }

  build_paired_plot <- function(metric, label) {
    wt_p <- tryCatch({
      pre_vals  <- paired_df[[metric]][paired_df$period == "Pre"]
      post_vals <- paired_df[[metric]][paired_df$period == "Post"]
      subj_paired <- intersect(
        paired_df[[subj_var]][paired_df$period == "Pre"],
        paired_df[[subj_var]][paired_df$period == "Post"]
      )
      pre_ord  <- paired_df[[metric]][paired_df[[subj_var]] %in% subj_paired &
                                        paired_df$period == "Pre"]
      post_ord <- paired_df[[metric]][paired_df[[subj_var]] %in% subj_paired &
                                        paired_df$period == "Post"]
      wt <- wilcox.test(pre_ord, post_ord, paired = TRUE, exact = FALSE)
      paste0("Wilcoxon p = ", signif(wt$p.value, 3))
    }, error = function(e) "")

    p <- ggplot(paired_df,
                aes(x = period, y = .data[[metric]],
                    group = .data[[subj_var]],
                    colour = if (!is.null(group_var)) .data[[group_var]] else "All")) +
      geom_line(alpha = 0.4, linewidth = 0.6) +
      geom_point(size = 2.5, alpha = 0.7) +
      stat_summary(aes(group = 1), fun = mean, geom = "line",
                   colour = "black", linewidth = 1.5, linetype = "dashed") +
      stat_summary(aes(group = 1), fun = mean, geom = "point",
                   colour = "black", size = 4, shape = 18) +
      annotate("text", x = 1.5, y = Inf, vjust = 1.5,
               label = wt_p, size = 3, colour = "grey30") +
      scale_colour_manual(values = colours,
                          name = if (!is.null(group_var)) group_var else NULL) +
      labs(
        title = label,
        x     = NULL, y = label
      ) +
      theme_microbiome()

    p
  }

  p_shannon  <- build_paired_plot("shannon", "Shannon entropy")
  p_observed <- build_paired_plot("observed", "Observed richness")

  combined <- (p_shannon | p_observed) +
    plot_annotation(
      title    = paste0(intervention_label, " - diversity response"),
      subtitle = paste0("Pre: T", paste(pre_timepoints, collapse="/"),
                        " to Post: T", paste(post_timepoints, collapse="/")),
      theme    = theme(plot.title = element_text(size = 13, face = "bold"))
    )

  # --- Responder classification based on Shannon change -------------------
  resp_df <- paired_df %>%
    select(.data[[subj_var]], period, shannon,
           any_of(group_var)) %>%
    pivot_wider(names_from = period, values_from = shannon,
                names_prefix = "shannon_") %>%
    mutate(
      delta_shannon = shannon_Post - shannon_Pre,
      responder     = delta_shannon > 0
    ) %>%
    filter(!is.na(delta_shannon))

  n_resp    <- sum(resp_df$responder)
  n_nonresp <- sum(!resp_df$responder)
  cat("  Responders (Shannon increased):     ", n_resp, "\n")
  cat("  Non-responders (Shannon decreased):", n_nonresp, "\n\n")

  p_resp <- ggplot(resp_df,
                   aes(x = delta_shannon,
                       fill = responder)) +
    geom_histogram(bins = 15, colour = "white", alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey30", linewidth = 0.9) +
    scale_fill_manual(
      values = c(`TRUE` = "#27ae60", `FALSE` = "#e74c3c"),
      labels = c(`TRUE` = "Responder", `FALSE` = "Non-responder"),
      name   = NULL
    ) +
    labs(
      title    = "Response classification",
      subtitle = paste0("Based on change in Shannon entropy | ",
                        n_resp, " responders, ", n_nonresp, " non-responders"),
      x        = "Delta Shannon (Post - Pre)",
      y        = "Number of subjects"
    ) +
    theme_microbiome()

  return(list(
    plot       = combined / p_resp +
                 plot_layout(heights = c(2, 1)),
    p_paired   = combined,
    p_response = p_resp,
    responders = resp_df
  ))
}


# =============================================================================
# SECTION 8 — CHANGEPOINT DETECTION
# =============================================================================

#' Detect changepoints in community diversity trajectories.
#'
#' Uses rolling window statistics to identify time points where
#' the microbiome undergoes abrupt shifts.
#'
#' @param long_data   Output from prepare_longitudinal_data().
#' @param metric      Diversity metric to analyse. Default = "shannon".
#' @param window      Rolling window size. Default = 2.
#' @return A list: changepoint plot and summary table.

detect_changepoints <- function(long_data,
                                 metric = "shannon",
                                 window = 2) {

  cat("=== Changepoint detection ===\n")

  div_df    <- long_data$diversity
  time_var  <- long_data$time_var
  subj_var  <- long_data$subject_var
  group_var <- long_data$group_var

  # Per-subject per-timepoint mean
  time_summary <- div_df %>%
    group_by(.data[[time_var]],
             if (!is.null(group_var)) .data[[group_var]] else NULL) %>%
    summarise(
      mean_val = mean(.data[[metric]], na.rm = TRUE),
      sd_val   = sd(.data[[metric]], na.rm = TRUE),
      n        = n(),
      .groups  = "drop"
    )

  if (!is.null(group_var)) {
    colnames(time_summary)[2] <- group_var
  }

  # Rolling mean and SD to detect abrupt changes
  time_summary <- time_summary %>%
    arrange(if (!is.null(group_var)) .data[[group_var]] else NULL,
            .data[[time_var]]) %>%
    group_by(if (!is.null(group_var)) .data[[group_var]] else NULL) %>%
    mutate(
      rolling_mean = rollmean(mean_val, k = window, fill = NA, align = "right"),
      abs_change   = abs(mean_val - lag(mean_val)),
      z_change     = (abs_change - mean(abs_change, na.rm = TRUE)) /
                     (sd(abs_change, na.rm = TRUE) + 1e-10),
      is_changepoint = z_change > 1.5
    ) %>%
    ungroup()

  n_changepoints <- sum(time_summary$is_changepoint, na.rm = TRUE)
  cat("  Changepoints detected:", n_changepoints, "\n\n")

  n_groups <- if (!is.null(group_var)) n_distinct(time_summary[[group_var]]) else 1
  colours  <- brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)]

  metric_label <- c(
    shannon  = "Shannon entropy (H')",
    observed = "Observed richness",
    pielou   = "Pielou's evenness"
  )[[metric]] %||% metric

  p <- ggplot(time_summary,
              aes(x = .data[[time_var]],
                  y = mean_val,
                  colour = if (!is.null(group_var)) .data[[group_var]] else "All",
                  group  = if (!is.null(group_var)) .data[[group_var]] else "All")) +
    geom_ribbon(
      aes(ymin  = mean_val - sd_val / sqrt(n),
          ymax  = mean_val + sd_val / sqrt(n),
          fill  = if (!is.null(group_var)) .data[[group_var]] else "All"),
      alpha = 0.15, colour = NA
    ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    # Highlight changepoints
    geom_point(
      data = time_summary %>% filter(is_changepoint == TRUE),
      aes(x = .data[[time_var]], y = mean_val),
      shape = 23, size = 5, fill = "#f39c12", colour = "#e67e22",
      stroke = 1.5, inherit.aes = FALSE
    ) +
    geom_text(
      data = time_summary %>% filter(is_changepoint == TRUE),
      aes(x = .data[[time_var]],
          y = mean_val,
          label = "WARNING: Shift"),
      vjust = -1.3, size = 3, colour = "#e67e22", fontface = "bold",
      inherit.aes = FALSE
    ) +
    scale_colour_manual(values = colours,
                        name = if (!is.null(group_var)) group_var else NULL) +
    scale_fill_manual(values = colours, guide = "none") +
    scale_x_continuous(breaks = long_data$time_points) +
    labs(
      title    = paste0("Changepoint detection - ", metric_label),
      subtitle = paste0("Orange diamonds = detected shifts (|z| > 1.5) | ",
                        n_changepoints, " changepoints"),
      x        = time_var,
      y        = metric_label,
      caption  = "Shaded area = +/-1 SE of the mean"
    ) +
    theme_microbiome()

  return(list(
    plot  = p,
    data  = time_summary,
    n_changepoints = n_changepoints
  ))
}


# =============================================================================
# SECTION 9 — COMPLETE LONGITUDINAL WORKFLOW WRAPPER
# =============================================================================

#' Run the complete longitudinal microbiome analysis pipeline.
#'
#' @param ps              A filtered phyloseq object (raw counts).
#' @param time_var        Metadata column for time. Must be numeric.
#' @param subject_var     Metadata column for subject ID.
#' @param group_var       Optional grouping variable.
#' @param rank            Taxonomic rank. Default = "Genus".
#' @param pre_timepoints  Time points before intervention (optional).
#' @param post_timepoints Time points after intervention (optional).
#' @param intervention_label Label for intervention (optional).
#' @param output_dir      Directory for outputs.
#' @return A named list of all results.

run_longitudinal_analysis <- function(ps,
                                       time_var           = "timepoint",
                                       subject_var        = "subject_id",
                                       group_var          = NULL,
                                       rank               = "Genus",
                                       pre_timepoints     = NULL,
                                       post_timepoints    = NULL,
                                       intervention_label = "Intervention",
                                       output_dir         = "longitudinal_output") {

  cat("\n", strrep("=", 60), "\n")
  cat("  LONGITUDINAL ANALYSIS PIPELINE\n")
  cat(strrep("=", 60), "\n\n")

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  results <- list()

  # --- Prepare data ---------------------------------------------------------
  long_data <- prepare_longitudinal_data(
    ps, time_var = time_var, subject_var = subject_var,
    group_var = group_var, rank = rank
  )
  results$long_data <- long_data

  # --- Alpha diversity trajectories ----------------------------------------
  cat("--- Plot 1: Alpha diversity trajectories ---\n")
  alpha_res <- plot_diversity_trajectories(long_data, test = TRUE)
  results$p_alpha <- alpha_res$plot
  ggsave(file.path(output_dir, "01_alpha_diversity_trajectories.pdf"),
         alpha_res$plot, width = 16, height = 6)

  # --- Beta diversity trajectories -----------------------------------------
  cat("--- Plot 2: Beta diversity trajectories ---\n")
  beta_res <- plot_beta_trajectories(long_data, distance = "bray")
  results$p_beta  <- beta_res$plot
  ggsave(file.path(output_dir, "02_beta_diversity_trajectories.pdf"),
         beta_res$plot, width = 12, height = 14)

  # --- Stability -----------------------------------------------------------
  cat("--- Plot 3: Stability analysis ---\n")
  stab_res <- analyse_stability(long_data)
  results$p_stability <- stab_res$plot
  write.csv(stab_res$subject_stability,
            file.path(output_dir, "subject_stability.csv"),
            row.names = FALSE)
  ggsave(file.path(output_dir, "03_stability.pdf"),
         stab_res$plot, width = 14, height = 6)

  # --- Composition over time -----------------------------------------------
  cat("--- Plot 4: Composition over time ---\n")
  comp_res <- plot_composition_over_time(long_data, rank = "Phylum", top_n = 10)
  results$p_composition <- comp_res$plot
  ggsave(file.path(output_dir, "04_composition_over_time.pdf"),
         comp_res$plot, width = 12, height = 12)

  # --- Mixed-effects models ------------------------------------------------
  cat("--- Analysis 1: Mixed-effects models ---\n")
  lme_res <- run_lme_over_time(long_data, interaction = !is.null(group_var))
  if (!is.null(lme_res)) {
    results$p_lme <- lme_res$plot
    write.csv(lme_res$results,
              file.path(output_dir, "lme_results.csv"),
              row.names = FALSE)
    ggsave(file.path(output_dir, "05_lme_forest_plot.pdf"),
           lme_res$plot, width = 10, height = max(6, lme_res$n_sig * 0.4))
  }

  # --- Changepoint detection -----------------------------------------------
  cat("--- Plot 6: Changepoint detection ---\n")
  cp_res <- detect_changepoints(long_data, metric = "shannon")
  results$p_changepoints <- cp_res$plot
  ggsave(file.path(output_dir, "06_changepoints.pdf"),
         cp_res$plot, width = 10, height = 6)

  # --- Intervention response (optional) ------------------------------------
  if (!is.null(pre_timepoints) && !is.null(post_timepoints)) {
    cat("--- Plot 7: Intervention response ---\n")
    int_res <- analyse_intervention_response(
      long_data,
      pre_timepoints     = pre_timepoints,
      post_timepoints    = post_timepoints,
      intervention_label = intervention_label
    )
    results$p_intervention <- int_res$plot
    write.csv(int_res$responders,
              file.path(output_dir, "intervention_responders.csv"),
              row.names = FALSE)
    ggsave(file.path(output_dir, "07_intervention_response.pdf"),
           int_res$plot, width = 12, height = 12)
  }

  cat("\n", strrep("=", 60), "\n")
  cat("  LONGITUDINAL PIPELINE COMPLETE\n")
  cat(strrep("=", 60), "\n")
  cat("  Output saved to:", output_dir, "\n")
  cat("  Plots:  up to 7 PDF files\n")
  cat("  Tables: stability, LME results, responders\n\n")

  return(invisible(results))
}


# Helper operator
`%||%` <- function(a, b) if (!is.null(a)) a else b


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# ps <- readRDS("qc_output/phyloseq_qc_filtered.rds")

# --- Option A: Full pipeline -------------------------------------------------
# results <- run_longitudinal_analysis(
#   ps                 = ps,
#   time_var           = "week",
#   subject_var        = "patient_id",
#   group_var          = "treatment",
#   rank               = "Genus",
#   pre_timepoints     = c(0, 1),
#   post_timepoints    = c(4, 8, 12),
#   intervention_label = "Antibiotic treatment",
#   output_dir         = "results/longitudinal"
# )

# --- Option B: Step by step --------------------------------------------------
# long_data <- prepare_longitudinal_data(ps, time_var = "week",
#               subject_var = "patient_id", group_var = "treatment")
# plot_diversity_trajectories(long_data)$plot
# plot_beta_trajectories(long_data)$plot
# analyse_stability(long_data)$plot
# run_lme_over_time(long_data)$plot

# --- Option C: Changepoint only ----------------------------------------------
# long_data <- prepare_longitudinal_data(ps, time_var = "week",
#               subject_var = "patient_id")
# cp_res    <- detect_changepoints(long_data, metric = "shannon")
# cp_res$plot
# cp_res$data %>% filter(is_changepoint)
