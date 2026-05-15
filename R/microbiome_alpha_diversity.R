# =============================================================================
# MICROBIOME ALPHA DIVERSITY ANALYSIS
# =============================================================================
# Description : Comprehensive alpha diversity analysis pipeline
# Input       : Filtered phyloseq object (output from microbiome_qc.R)
# Output      : Diversity metric plots, statistical tests, rarefaction analysis
# Author      : Patricia
# Dependencies: phyloseq, vegan, ggplot2, dplyr, tidyr, patchwork,
#               rstatix, ggpubr, iNEXT, picante (optional, for Faith's PD)
# =============================================================================

# --- 1. LOAD LIBRARIES -------------------------------------------------------

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(rstatix)      # Tidy statistical tests
  library(ggpubr)       # Significance brackets on plots
  library(scales)
  library(tibble)
  library(forcats)
  library(RColorBrewer)
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
# SECTION 1 — CALCULATE ALPHA DIVERSITY METRICS
# =============================================================================

#' Calculate a comprehensive set of alpha diversity metrics for each sample.
#'
#' Metrics included:
#'   - Observed richness   : number of ASVs/OTUs detected
#'   - Chao1               : estimated true richness including rare taxa
#'   - ACE                 : abundance-based coverage estimator
#'   - Shannon entropy     : richness + evenness (log base e)
#'   - Simpson (1-D)       : dominance-corrected diversity
#'   - InvSimpson          : inverse Simpson index
#'   - Pielou's evenness   : Shannon / log(richness)
#'   - Fisher's alpha      : log-series diversity parameter
#'   - Faith's PD          : phylogenetic diversity (requires tree)
#'
#' @param ps          A phyloseq object. Counts (not relative abundance).
#' @param rarefaction Whether to rarefy before calculating metrics. Default = TRUE.
#' @param rare_depth  Rarefaction depth. Default = minimum sample depth.
#' @param seed        Random seed for rarefaction reproducibility. Default = 42.
#' @param n_rare      Number of rarefaction iterations for averaging. Default = 10.
#' @return A data frame with one row per sample and all diversity metrics.

calculate_alpha_diversity <- function(ps,
                                       rarefaction = TRUE,
                                       rare_depth  = NULL,
                                       seed        = 42,
                                       n_rare      = 10) {

  cat("=== Calculating alpha diversity metrics ===\n")

  # Ensure we are working with counts not relative abundance
  if (max(sample_sums(ps)) <= 1) {
    stop("Alpha diversity requires raw count data, not relative abundance. ",
         "Please provide the unrarefied count phyloseq object.")
  }

  otu_mat <- as.matrix(otu_table(ps))
  if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  otu_t   <- matrix(t(otu_mat), nrow=ncol(otu_mat), ncol=nrow(otu_mat),
                   dimnames=list(colnames(otu_mat), rownames(otu_mat)))   # samples as rows for vegan

  # --- Rarefaction ----------------------------------------------------------
  if (rarefaction) {
    if (is.null(rare_depth)) {
      rare_depth <- min(rowSums(otu_t))
      cat("  Rarefaction depth set to minimum sample size:", rare_depth, "\n")
    }

    # Warn if rarefaction removes samples
    n_drop <- sum(rowSums(otu_t) < rare_depth)
    if (n_drop > 0) {
      warning(n_drop, " sample(s) have fewer reads than the rarefaction depth ",
              "and will be excluded.")
      otu_t <- otu_t[rowSums(otu_t) >= rare_depth, ]
    }

    cat("  Averaging over", n_rare, "rarefaction iterations...\n")

    # Rarefy n_rare times and average
    rare_matrices <- lapply(seq_len(n_rare), function(i) {
      set.seed(seed + i)
      rrarefy(otu_t[rowSums(otu_t) >= rare_depth, , drop=FALSE], sample = rare_depth)
    })

    otu_t_rare <- Reduce("+", rare_matrices) / n_rare
  } else {
    otu_t_rare <- otu_t
    cat("  Skipping rarefaction (rarefaction = FALSE)\n")
  }

  # --- Standard vegan metrics -----------------------------------------------
  observed    <- rowSums(otu_t_rare > 0)
  shannon     <- vegan::diversity(otu_t_rare, index = "shannon")
  simpson     <- vegan::diversity(otu_t_rare, index = "simpson")
  invsimpson  <- vegan::diversity(otu_t_rare, index = "invsimpson")
  fisher_a    <- fisher.alpha(round(otu_t_rare))

  # Chao1 and ACE via estimateR — requires integer counts
  cat("  Computing Chao1 and ACE...\n")
  chao_ace <- tryCatch({
    t(estimateR(round(otu_t)))
  }, error = function(e) {
    cat("  Warning: Chao1/ACE estimation failed:", e$message, "\n")
    matrix(NA, nrow = nrow(otu_t_rare), ncol = 4,
           dimnames = list(rownames(otu_t_rare),
                           c("S.obs", "S.chao1", "S.ACE", "se.ACE")))
  })

  chao1 <- chao_ace[rownames(otu_t_rare), "S.chao1"]
  ace   <- chao_ace[rownames(otu_t_rare), "S.ACE"]

  # Pielou's evenness = Shannon / log(richness)
  pielou <- ifelse(observed > 1, shannon / log(observed), 0)

  # --- Faith's phylogenetic diversity (requires tree) -----------------------
  faith_pd <- rep(NA_real_, nrow(otu_t_rare))
  names(faith_pd) <- rownames(otu_t_rare)

  if (!is.null(phy_tree(ps, errorIfNULL = FALSE))) {
    if (requireNamespace("picante", quietly = TRUE)) {
      cat("  Computing Faith's PD...\n")
      tree     <- phy_tree(ps)
      pd_res   <- picante::pd(round(otu_t_rare), tree, include.root = TRUE)
      faith_pd <- pd_res[rownames(otu_t_rare), "PD"]
    } else {
      cat("  picante not installed — skipping Faith's PD\n")
    }
  } else {
    cat("  No phylogenetic tree found — skipping Faith's PD\n")
    faith_pd <- rep(NA_real_, nrow(otu_t_rare))
  }

  # --- Compile into data frame ---------------------------------------------
  diversity_df <- data.frame(
    sample      = rownames(otu_t_rare),
    observed    = as.integer(observed),
    chao1       = round(chao1, 1),
    ace         = round(ace, 1),
    shannon     = round(shannon, 4),
    simpson     = round(simpson, 4),
    invsimpson  = round(invsimpson, 4),
    pielou      = round(pielou, 4),
    fisher      = round(fisher_a, 4),
    faith_pd    = round(faith_pd, 4),
    stringsAsFactors = FALSE
  )

  # Add sample metadata
  meta_df <- data.frame(sample_data(ps)) %>% rownames_to_column("sample")
  diversity_df <- left_join(diversity_df, meta_df, by = "sample")

  cat("  Metrics calculated for", nrow(diversity_df), "samples\n")
  cat("  Metrics: observed, chao1, ace, shannon, simpson, invsimpson,",
      "pielou, fisher, faith_pd\n\n")

  return(diversity_df)
}


# =============================================================================
# SECTION 2 — ALPHA DIVERSITY PLOTS
# =============================================================================

#' Plot alpha diversity distributions per group with statistical testing.
#'
#' @param diversity_df  Data frame from calculate_alpha_diversity().
#' @param metrics       Character vector of metrics to plot.
#' @param group_var     Metadata variable defining groups.
#' @param test          Statistical test: "wilcox", "kruskal", or "anova".
#' @param p_adjust      P-value adjustment method. Default = "BH" (Benjamini-Hochberg).
#' @param show_points   Whether to overlay individual sample points. Default = TRUE.
#' @param colour_var    Metadata variable for colouring points. Default = group_var.
#' @return A patchwork plot combining all metric panels.

plot_alpha_diversity <- function(diversity_df,
                                  metrics    = c("observed", "shannon",
                                                 "simpson",  "pielou"),
                                  group_var  = "group",
                                  test       = "wilcox",
                                  p_adjust   = "BH",
                                  show_points = TRUE,
                                  colour_var  = NULL) {

  cat("=== Alpha diversity plots ===\n")

  if (is.null(colour_var)) colour_var <- group_var

  # Validate metrics
  available <- c("observed", "chao1", "ace", "shannon",
                 "simpson", "invsimpson", "pielou", "fisher", "faith_pd")
  metrics <- intersect(metrics, colnames(diversity_df))
  if (length(metrics) == 0) stop("No valid metrics found in diversity_df.")

  # Friendly metric labels
  metric_labels <- c(
    observed   = "Observed richness",
    chao1      = "Chao1 richness",
    ace        = "ACE richness",
    shannon    = "Shannon entropy (H')",
    simpson    = "Simpson diversity (1-D)",
    invsimpson = "Inverse Simpson (1/D)",
    pielou     = "Pielou's evenness (J')",
    fisher     = "Fisher's alpha",
    faith_pd   = "Faith's phylogenetic diversity"
  )

  groups   <- unique(diversity_df[[group_var]])
  n_groups <- length(groups)
  colours  <- brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)]
  names(colours) <- groups

  # --- Statistical test per metric -----------------------------------------
  stat_results <- lapply(metrics, function(m) {

    formula_str <- as.formula(paste(m, "~", group_var))

    if (test == "wilcox" && n_groups == 2) {
      stat_test <- diversity_df %>%
        wilcox_test(formula_str) %>%
        adjust_pvalue(method = p_adjust) %>%
        add_significance("p.adj")

    } else if (test == "kruskal" || (test == "wilcox" && n_groups > 2)) {
      # Kruskal overall + pairwise Dunn
      overall <- kruskal_test(diversity_df, formula_str)
      stat_test <- diversity_df %>%
        dunn_test(formula_str, p.adjust.method = p_adjust) %>%
        add_significance("p.adj")
      attr(stat_test, "kruskal_p") <- overall$p

    } else if (test == "anova") {
      stat_test <- diversity_df %>%
        tukey_hsd(formula_str) %>%
        add_significance("p.adj")
    }

    stat_test$.metric <- m
    stat_test
  })

  names(stat_results) <- metrics

  # --- Build individual metric plots ----------------------------------------
  plots <- lapply(metrics, function(m) {

    y_vals  <- diversity_df[[m]]
    y_range <- range(y_vals, na.rm = TRUE)
    y_pad   <- diff(y_range) * 0.15

    p <- ggplot(diversity_df,
                aes(x = .data[[group_var]],
                    y = .data[[m]],
                    fill = .data[[group_var]])) +
      geom_boxplot(
        alpha         = 0.75,
        outlier.shape = NA,
        width         = 0.55,
        linewidth     = 0.5
      ) +
      {
        if (show_points) {
          geom_jitter(
            aes(colour = .data[[colour_var]]),
            width = 0.15, alpha = 0.6, size = 1.8
          )
        }
      } +
      scale_fill_manual(values = colours, guide = "none") +
      scale_colour_manual(values = colours, guide = "none") +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      labs(
        title = metric_labels[[m]],
        x     = NULL,
        y     = metric_labels[[m]]
      ) +
      theme_microbiome() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    # Add significance brackets
    stat_df <- stat_results[[m]]
    if (nrow(stat_df) > 0) {
      sig_df <- stat_df %>% filter(p.adj.signif != "ns")
      if (nrow(sig_df) > 0) {
        p <- p + stat_pvalue_manual(
          sig_df,
          label      = "p.adj.signif",
          y.position = seq(y_range[2] + y_pad,
                           by = y_pad,
                           length.out = nrow(sig_df)),
          tip.length = 0.01,
          size       = 3
        )
      }
    }

    p
  })

  # Combine all panels
  n_cols  <- min(4, length(plots))
  combined <- wrap_plots(plots, ncol = n_cols) +
    plot_annotation(
      title    = "Alpha Diversity",
      subtitle = paste0("n = ", nrow(diversity_df), " samples | ",
                        n_groups, " groups | Test: ", test,
                        " + ", p_adjust, " correction"),
      theme    = theme(
        plot.title    = element_text(size = 15, face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey40")
      )
    )

  cat("  Plots generated for:", paste(metrics, collapse = ", "), "\n\n")
  return(list(plot = combined, stats = stat_results))
}


# =============================================================================
# SECTION 3 — RAREFACTION CURVES PER SAMPLE
# =============================================================================

#' Generate individual rarefaction curves showing how richness accumulates.
#'
#' @param ps          A phyloseq object (raw counts).
#' @param step        Rarefaction step size. Default = 500.
#' @param group_var   Metadata variable for colouring curves.
#' @param show_se     Whether to show standard error ribbon. Default = TRUE.
#' @param n_iter      Number of iterations per depth. Default = 5.
#' @return A list: plot and rarefaction data frame.

plot_rarefaction_curves <- function(ps,
                                     step     = 500,
                                     group_var = NULL,
                                     show_se   = TRUE,
                                     n_iter    = 5) {

  cat("=== Rarefaction curves per sample ===\n")

  otu_mat <- as.matrix(otu_table(ps))
  if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  otu_t   <- matrix(t(otu_mat), nrow=ncol(otu_mat), ncol=nrow(otu_mat),
                   dimnames=list(colnames(otu_mat), rownames(otu_mat)))

  max_depth <- max(rowSums(otu_t))
  depths    <- c(seq(step, max_depth, by = step))
  if (tail(depths, 1) < max_depth) depths <- c(depths, max_depth)

  cat("  Samples:", nrow(otu_t), "| Max depth:", max_depth,
      "| Steps:", length(depths), "\n")

  # Compute rarefaction with multiple iterations for SE
  rare_list <- lapply(depths, function(d) {
    iter_list <- lapply(seq_len(n_iter), function(i) {
      set.seed(i)
      subs <- otu_t[rowSums(otu_t) >= d, , drop = FALSE]
      if (nrow(subs) == 0) return(NULL)
      rich <- rowSums(rrarefy(subs, sample = d) > 0)
      data.frame(
        sample   = names(rich),
        depth    = d,
        richness = as.numeric(rich),
        iter     = i
      )
    })
    bind_rows(iter_list)
  })

  rare_df <- bind_rows(rare_list)

  # Average across iterations
  rare_summary <- rare_df %>%
    group_by(sample, depth) %>%
    summarise(
      mean_richness = mean(richness),
      se_richness   = sd(richness) / sqrt(n()),
      .groups = "drop"
    )

  # Add metadata
  if (!is.null(group_var) && group_var %in% sample_variables(ps)) {
    meta_df <- data.frame(sample_data(ps)) %>%
      select(all_of(group_var)) %>%
      rownames_to_column("sample")
    rare_summary <- left_join(rare_summary, meta_df, by = "sample")
  }

  # Mark samples that did not reach maximum depth
  max_per_sample <- rare_summary %>%
    group_by(sample) %>%
    summarise(max_depth_reached = max(depth), .groups = "drop")

  plateau_check <- max_per_sample %>%
    filter(max_depth_reached < max_depth * 0.9) %>%
    pull(sample)

  if (length(plateau_check) > 0) {
    cat("  Samples that may not have reached plateau (< 90% of max depth):\n")
    cat("   ", paste(plateau_check, collapse = ", "), "\n")
  }

  # --- Plot -----------------------------------------------------------------
  colour_var <- if (!is.null(group_var) && group_var %in% colnames(rare_summary)) group_var else NULL

  p <- ggplot(rare_summary,
              aes(x = depth, y = mean_richness,
                  group = sample,
                  colour = if (!is.null(colour_var)) .data[[colour_var]] else sample)) +
    {
      if (show_se) {
        geom_ribbon(
          aes(ymin = mean_richness - se_richness,
              ymax = mean_richness + se_richness,
              fill = if (!is.null(colour_var)) .data[[colour_var]] else sample),
          alpha = 0.08, colour = NA
        )
      }
    } +
    geom_line(alpha = 0.6, linewidth = 0.5) +
    geom_vline(xintercept = min(rowSums(otu_t)),
               linetype = "dashed", colour = "#e74c3c",
               linewidth = 0.8) +
    annotate("text",
             x = min(rowSums(otu_t)) * 1.02,
             y = max(rare_summary$mean_richness) * 0.95,
             label = paste0("Min depth\n",
                            format(min(rowSums(otu_t)), big.mark = ",")),
             hjust = 0, size = 3, colour = "#e74c3c") +
    scale_x_continuous(labels = label_comma()) +
    scale_colour_brewer(palette = "Set2", name = colour_var) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(
      title    = "Rarefaction curves — observed richness",
      subtitle = paste0(nrow(otu_t), " samples | ",
                        n_iter, " iterations averaged"),
      x        = "Sequencing depth (reads)",
      y        = "Observed ASV richness",
      caption  = if (length(plateau_check) > 0)
        paste0("Samples not reaching plateau: ",
               paste(plateau_check, collapse = ", "))
      else "All samples appear to have reached plateau"
    ) +
    theme_microbiome()

  return(list(plot = p, data = rare_summary))
}


# =============================================================================
# SECTION 4 — DIVERSITY OVER METADATA GRADIENTS
# =============================================================================

#' Plot alpha diversity as a function of a continuous or ordinal metadata variable.
#'
#' @param diversity_df  Data frame from calculate_alpha_diversity().
#' @param metric        Diversity metric to plot.
#' @param x_var         Continuous or ordinal metadata variable for x axis.
#' @param group_var     Optional grouping variable for faceting or colouring.
#' @param add_smooth    Whether to add a loess/lm smoothing line. Default = TRUE.
#' @param smooth_method Smoothing method: "loess" or "lm". Default = "loess".
#' @return A ggplot object.

plot_diversity_gradient <- function(diversity_df,
                                     metric       = "shannon",
                                     x_var        = "age",
                                     group_var    = NULL,
                                     add_smooth   = TRUE,
                                     smooth_method = "loess") {

  cat("=== Diversity gradient plot:", metric, "~", x_var, "===\n")

  if (!x_var %in% colnames(diversity_df)) {
    stop("Variable '", x_var, "' not found in diversity data frame.")
  }

  metric_labels <- c(
    observed   = "Observed richness",
    shannon    = "Shannon entropy (H')",
    simpson    = "Simpson diversity (1-D)",
    pielou     = "Pielou's evenness (J')",
    chao1      = "Chao1 richness",
    faith_pd   = "Faith's phylogenetic diversity"
  )

  y_label <- if (metric %in% names(metric_labels)) metric_labels[[metric]] else metric

  # Correlation coefficient
  cor_val  <- tryCatch(
    cor.test(diversity_df[[x_var]], diversity_df[[metric]],
             method = "spearman"),
    error = function(e) NULL
  )
  cor_label <- if (!is.null(cor_val)) {
    paste0("Spearman rho = ", round(cor_val$estimate, 3),
           ", p = ", signif(cor_val$p.value, 3))
  } else ""

  p <- ggplot(diversity_df,
              aes(x = .data[[x_var]], y = .data[[metric]])) +
    {
      if (!is.null(group_var) && group_var %in% colnames(diversity_df)) {
        list(
          geom_point(aes(colour = .data[[group_var]]), alpha = 0.7, size = 2),
          scale_colour_brewer(palette = "Set2", name = group_var)
        )
      } else {
        list(geom_point(colour = "#3498db", alpha = 0.7, size = 2))
      }
    } +
    {
      if (add_smooth) {
        geom_smooth(method  = smooth_method,
                    se      = TRUE,
                    colour  = "#e74c3c",
                    fill    = "#e74c3c",
                    alpha   = 0.12,
                    linewidth = 0.9)
      }
    } +
    {
      if (!is.null(group_var) && group_var %in% colnames(diversity_df)) {
        facet_wrap(~ .data[[group_var]])
      }
    } +
    labs(
      title    = paste0(y_label, " ~ ", x_var),
      subtitle = cor_label,
      x        = x_var,
      y        = y_label
    ) +
    theme_microbiome()

  return(p)
}


# =============================================================================
# SECTION 5 — LONGITUDINAL ALPHA DIVERSITY
# =============================================================================

#' Plot alpha diversity trajectories over time for paired or longitudinal data.
#'
#' @param diversity_df  Data frame from calculate_alpha_diversity().
#' @param metric        Diversity metric.
#' @param time_var      Metadata variable for time points.
#' @param subject_var   Metadata variable for subject/individual ID.
#' @param group_var     Optional grouping variable.
#' @return A ggplot object.

plot_longitudinal_diversity <- function(diversity_df,
                                         metric      = "shannon",
                                         time_var    = "timepoint",
                                         subject_var = "subject_id",
                                         group_var   = NULL) {

  cat("=== Longitudinal diversity plot ===\n")

  required_vars <- c(time_var, subject_var)
  missing_vars  <- setdiff(required_vars, colnames(diversity_df))
  if (length(missing_vars) > 0) {
    stop("Missing variables: ", paste(missing_vars, collapse = ", "))
  }

  metric_labels <- c(
    observed = "Observed richness", shannon = "Shannon entropy (H')",
    simpson  = "Simpson diversity (1-D)", pielou = "Pielou's evenness (J')"
  )
  y_label <- if (metric %in% names(metric_labels)) metric_labels[[metric]] else metric

  # Order time variable
  diversity_df[[time_var]] <- factor(diversity_df[[time_var]])

  # Summary per time point and group
  group_summary <- diversity_df %>%
    group_by(.data[[time_var]],
             if (!is.null(group_var)) .data[[group_var]] else NULL) %>%
    summarise(
      mean_div = mean(.data[[metric]], na.rm = TRUE),
      se_div   = sd(.data[[metric]], na.rm = TRUE) / sqrt(n()),
      .groups  = "drop"
    )

  if (!is.null(group_var) && group_var %in% colnames(group_summary)) {
    colnames(group_summary)[2] <- group_var
  }

  p <- ggplot(diversity_df,
              aes(x = .data[[time_var]],
                  y = .data[[metric]],
                  group = .data[[subject_var]])) +
    # Individual trajectories (thin, faint)
    geom_line(alpha = 0.2, linewidth = 0.4, colour = "grey60") +
    geom_point(alpha = 0.25, size = 1.2, colour = "grey60") +
    # Group mean ± SE on top
    {
      if (!is.null(group_var) && group_var %in% colnames(group_summary)) {
        list(
          geom_line(
            data    = group_summary,
            aes(x = .data[[time_var]], y = mean_div,
                group = .data[[group_var]],
                colour = .data[[group_var]]),
            linewidth = 1.2, inherit.aes = FALSE
          ),
          geom_ribbon(
            data  = group_summary,
            aes(x = .data[[time_var]],
                ymin  = mean_div - se_div,
                ymax  = mean_div + se_div,
                group = .data[[group_var]],
                fill  = .data[[group_var]]),
            alpha = 0.2, inherit.aes = FALSE
          ),
          scale_colour_brewer(palette = "Set2", name = group_var),
          scale_fill_brewer(palette = "Set2", guide = "none")
        )
      } else {
        list(
          geom_line(
            data      = group_summary,
            aes(x = .data[[time_var]], y = mean_div, group = 1),
            colour    = "#e74c3c", linewidth = 1.2,
            inherit.aes = FALSE
          ),
          geom_ribbon(
            data    = group_summary,
            aes(x = .data[[time_var]],
                ymin  = mean_div - se_div,
                ymax  = mean_div + se_div,
                group = 1),
            fill    = "#e74c3c", alpha = 0.2, inherit.aes = FALSE
          )
        )
      }
    } +
    labs(
      title    = paste0(y_label, " over time"),
      subtitle = paste0("Individual trajectories (grey) + group mean ± SE"),
      x        = time_var,
      y        = y_label,
      caption  = paste0("n = ",
                        n_distinct(diversity_df[[subject_var]]),
                        " individuals | ",
                        nlevels(diversity_df[[time_var]]),
                        " time points")
    ) +
    theme_microbiome()

  return(p)
}


# =============================================================================
# SECTION 6 — STATISTICAL SUMMARY TABLE
# =============================================================================

#' Produce a publication-ready statistical summary of all alpha diversity metrics.
#'
#' @param diversity_df  Data frame from calculate_alpha_diversity().
#' @param group_var     Metadata variable defining groups.
#' @param test          "wilcox" or "kruskal". Default = "wilcox".
#' @param p_adjust      P-value adjustment method. Default = "BH".
#' @return A tidy data frame with group means, SDs, and test results.

summarise_alpha_stats <- function(diversity_df,
                                   group_var = "group",
                                   test      = "wilcox",
                                   p_adjust  = "BH") {

  cat("=== Alpha diversity statistical summary ===\n")

  metrics <- c("observed", "chao1", "shannon", "simpson",
               "invsimpson", "pielou", "fisher", "faith_pd")
  metrics <- intersect(metrics, colnames(diversity_df))

  # Descriptive statistics per group
  desc_stats <- diversity_df %>%
    select(all_of(c(group_var, metrics))) %>%
    pivot_longer(-all_of(group_var), names_to = "metric", values_to = "value") %>%
    filter(!is.na(value)) %>%
    group_by(.data[[group_var]], metric) %>%
    summarise(
      n      = n(),
      mean   = round(mean(value), 4),
      sd     = round(sd(value), 4),
      median = round(median(value), 4),
      iqr    = round(IQR(value), 4),
      .groups = "drop"
    )

  # Statistical tests
  n_groups <- n_distinct(diversity_df[[group_var]])

  stat_list <- lapply(metrics, function(m) {
    formula_str <- as.formula(paste(m, "~", group_var))

    if (n_groups == 2 && test == "wilcox") {
      result <- tryCatch({
        wt <- wilcox.test(formula_str, data = diversity_df, exact = FALSE)
        data.frame(
          metric   = m,
          test     = "Wilcoxon rank-sum",
          statistic = round(wt$statistic, 3),
          p_value  = round(wt$p.value, 5),
          stringsAsFactors = FALSE
        )
      }, error = function(e) NULL)
    } else {
      result <- tryCatch({
        kt <- kruskal.test(formula_str, data = diversity_df)
        data.frame(
          metric    = m,
          test      = "Kruskal-Wallis",
          statistic = round(kt$statistic, 3),
          p_value   = round(kt$p.value, 5),
          stringsAsFactors = FALSE
        )
      }, error = function(e) NULL)
    }
    result
  })

  stat_df <- bind_rows(stat_list)

  # Apply p-value correction
  stat_df$p_adjusted <- round(p.adjust(stat_df$p_value, method = p_adjust), 5)
  stat_df$significance <- case_when(
    stat_df$p_adjusted < 0.001 ~ "***",
    stat_df$p_adjusted < 0.01  ~ "**",
    stat_df$p_adjusted < 0.05  ~ "*",
    stat_df$p_adjusted < 0.1   ~ ".",
    TRUE                       ~ "ns"
  )

  cat("  Significant metrics (p.adj < 0.05):",
      sum(stat_df$p_adjusted < 0.05, na.rm = TRUE), "\n\n")

  print(stat_df %>% select(metric, test, statistic, p_value, p_adjusted, significance))
  cat("\n")

  return(list(
    descriptive = desc_stats,
    tests       = stat_df
  ))
}


# =============================================================================
# SECTION 7 — DIVERSITY CORRELATION MATRIX
# =============================================================================

#' Plot pairwise correlations between all alpha diversity metrics.
#'
#' @param diversity_df  Data frame from calculate_alpha_diversity().
#' @param method        Correlation method: "spearman" or "pearson". Default = "spearman".
#' @return A ggplot correlation matrix.

plot_diversity_correlations <- function(diversity_df,
                                         method = "spearman") {

  cat("=== Alpha diversity correlation matrix ===\n")

  metrics <- c("observed", "chao1", "ace", "shannon",
               "simpson", "invsimpson", "pielou", "fisher", "faith_pd")
  metrics <- intersect(metrics, colnames(diversity_df))

  metric_data <- diversity_df[, metrics, drop = FALSE]
  
  # Remove metrics that are completely NA or have too few values
  metric_data <- metric_data[, colSums(!is.na(metric_data)) > 1, drop = FALSE]
  
  # Remove metrics with zero variance
  metric_data <- metric_data[, sapply(metric_data, function(x) {
    sd(x, na.rm = TRUE) > 0
  }), drop = FALSE]
  
  # Update metric list after filtering
  metrics <- colnames(metric_data)
  
  # Keep only samples with complete data for the remaining metrics
  metric_data <- metric_data[complete.cases(metric_data), ]

  # Compute correlation matrix
  cor_mat <- cor(metric_data, method = method, use = "pairwise.complete.obs")

  # Melt to long format
  cor_df <- as.data.frame(cor_mat) %>%
    rownames_to_column("metric1") %>%
    pivot_longer(-metric1, names_to = "metric2", values_to = "correlation") %>%
    mutate(
      metric1 = factor(metric1, levels = metrics),
      metric2 = factor(metric2, levels = rev(metrics))
    )

  p <- ggplot(cor_df, aes(x = metric1, y = metric2, fill = correlation)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = round(correlation, 2)),
              size = 3, fontface = "bold",
              colour = ifelse(abs(cor_df$correlation) > 0.7, "white", "black")) +
    scale_fill_gradient2(
      low      = "#3498db",
      mid      = "white",
      high     = "#e74c3c",
      midpoint = 0,
      limits   = c(-1, 1),
      name     = paste0(str_to_title(method), "\ncorrelation")
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      title    = "Alpha diversity metric correlations",
      subtitle = paste0(str_to_title(method),
                        " correlation | n = ", nrow(metric_data), " samples"),
      x = NULL, y = NULL
    ) +
    theme_microbiome() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      panel.border = element_rect(colour = "grey80", fill = NA)
    )

  return(p)
}


# =============================================================================
# SECTION 8 — COMPLETE ALPHA DIVERSITY WORKFLOW WRAPPER
# =============================================================================

#' Run the complete alpha diversity pipeline.
#'
#' @param ps          A filtered phyloseq object (raw counts).
#' @param group_var   Metadata variable for group comparisons.
#' @param time_var    Optional time variable for longitudinal analysis.
#' @param subject_var Optional subject ID for longitudinal analysis.
#' @param rare_depth  Rarefaction depth. Default = minimum sample depth.
#' @param output_dir  Directory to save all outputs.
#' @return A named list of all results.

run_alpha_diversity <- function(ps,
                                 group_var   = NULL,
                                 time_var    = NULL,
                                 subject_var = NULL,
                                 rare_depth  = NULL,
                                 output_dir  = "alpha_output") {

  cat("\n", strrep("=", 60), "\n")
  cat("  ALPHA DIVERSITY PIPELINE\n")
  cat(strrep("=", 60), "\n\n")

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  results <- list()

  # --- Step 1: Calculate all metrics ----------------------------------------
  results$diversity_df <- calculate_alpha_diversity(
    ps          = ps,
    rarefaction = TRUE,
    rare_depth  = rare_depth
  )

  write.csv(results$diversity_df,
            file.path(output_dir, "alpha_diversity_metrics.csv"),
            row.names = FALSE)

  # --- Step 2: Main diversity plots -----------------------------------------
  cat("--- Plot 1: Alpha diversity by group ---\n")
  alpha_plots <- plot_alpha_diversity(
    diversity_df = results$diversity_df,
    metrics      = c("observed", "shannon", "simpson", "pielou", "chao1"),
    group_var    = group_var,
    test         = ifelse(
      n_distinct(results$diversity_df[[group_var]]) == 2,
      "wilcox", "kruskal"
    )
  )
  results$p_alpha      <- alpha_plots$plot
  results$alpha_stats  <- alpha_plots$stats
  ggsave(file.path(output_dir, "01_alpha_diversity_by_group.pdf"),
         results$p_alpha, width = 16, height = 8)

  # --- Step 3: Rarefaction curves -------------------------------------------
  cat("--- Plot 2: Rarefaction curves ---\n")
  rare_res <- plot_rarefaction_curves(ps, group_var = group_var)
  results$p_rarefaction <- rare_res$plot
  ggsave(file.path(output_dir, "02_rarefaction_curves.pdf"),
         results$p_rarefaction, width = 10, height = 7)

  # --- Step 4: Correlation matrix -------------------------------------------
  cat("--- Plot 3: Metric correlation matrix ---\n")
  results$p_correlation <- plot_diversity_correlations(results$diversity_df)
  ggsave(file.path(output_dir, "03_diversity_correlations.pdf"),
         results$p_correlation, width = 9, height = 8)

  # --- Step 5: Statistical summary ------------------------------------------
  cat("--- Table 1: Statistical summary ---\n")
  if (!is.null(group_var)) {
    stat_summary <- summarise_alpha_stats(
      results$diversity_df,
      group_var = group_var
    )
    results$stat_summary <- stat_summary
    write.csv(stat_summary$tests,
              file.path(output_dir, "alpha_statistical_tests.csv"),
              row.names = FALSE)
    write.csv(stat_summary$descriptive,
              file.path(output_dir, "alpha_descriptive_stats.csv"),
              row.names = FALSE)
  }

  # --- Step 6: Longitudinal (optional) --------------------------------------
  if (!is.null(time_var) && !is.null(subject_var)) {
    cat("--- Plot 4: Longitudinal diversity ---\n")
    results$p_longitudinal <- plot_longitudinal_diversity(
      results$diversity_df,
      metric      = "shannon",
      time_var    = time_var,
      subject_var = subject_var,
      group_var   = group_var
    )
    ggsave(file.path(output_dir, "04_longitudinal_diversity.pdf"),
           results$p_longitudinal, width = 10, height = 6)
  }

  cat("\n", strrep("=", 60), "\n")
  cat("  ALPHA DIVERSITY PIPELINE COMPLETE\n")
  cat(strrep("=", 60), "\n")
  cat("  Output saved to:", output_dir, "\n")
  cat("  Metrics CSV    : alpha_diversity_metrics.csv\n")
  cat("  Stats CSV      : alpha_statistical_tests.csv\n")
  cat("  Plots saved    :", ifelse(!is.null(time_var), 4, 3), "PDF files\n\n")

  return(invisible(results))
}


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# Load filtered phyloseq from QC module (must be raw counts)
# ps <- readRDS("qc_output/phyloseq_qc_filtered.rds")

# --- Option A: Full pipeline -------------------------------------------------
# results <- run_alpha_diversity(
#   ps          = ps,
#   group_var   = "disease_status",
#   rare_depth  = 10000,
#   output_dir  = "results/alpha"
# )

# --- Option B: Longitudinal study --------------------------------------------
# results <- run_alpha_diversity(
#   ps          = ps,
#   group_var   = "treatment",
#   time_var    = "week",
#   subject_var = "patient_id",
#   output_dir  = "results/alpha_longitudinal"
# )

# --- Option C: Individual steps ----------------------------------------------
# div_df <- calculate_alpha_diversity(ps, rare_depth = 10000)
# div_df %>% select(sample, observed, shannon, pielou) %>% head()
#
# plot_alpha_diversity(div_df, group_var = "disease_status")$plot
# plot_diversity_gradient(div_df, metric = "shannon", x_var = "age",
#                         group_var = "disease_status")
# plot_diversity_correlations(div_df)
