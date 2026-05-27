# =============================================================================
# MICROBIOME DIFFERENTIAL ABUNDANCE ANALYSIS
# =============================================================================
# Description : Comprehensive differential abundance pipeline using multiple
#               methods with consensus scoring and visualisation
# Input       : Filtered phyloseq object (output from microbiome_qc.R)
# Output      : Differential taxa tables, volcano plots, effect size plots,
#               heatmaps, multi-method consensus results
# Author      : Patricia
# Dependencies: phyloseq, ANCOMBC, DESeq2, vegan, ggplot2, dplyr, tidyr,
#               patchwork, scales, tibble, RColorBrewer, microbiome,
#               ALDEx2 (optional), Maaslin2 (optional)
# =============================================================================

# --- 1. LOAD LIBRARIES -------------------------------------------------------

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(tibble)
  library(RColorBrewer)
  library(microbiome)
  library(stringr)
  library(forcats)
})

# Conditionally load method-specific packages
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
# SECTION 1 — DATA PREPARATION
# =============================================================================

#' Prepare a phyloseq object for differential abundance analysis.
#'
#' Agglomerates taxa, removes zero-variance taxa, and applies minimum
#' prevalence filtering appropriate for DA methods.
#'
#' @param ps              A phyloseq object (raw counts).
#' @param rank            Taxonomic rank. Default = "Genus".
#' @param min_prevalence  Minimum fraction of samples a taxon must appear in.
#' @param min_count       Minimum total read count across all samples.
#' @return A filtered, agglomerated phyloseq object.

prepare_da_data <- function(ps,
                             rank           = "Genus",
                             min_prevalence = 0.10,
                             min_count      = 10) {

  cat("=== Preparing data for differential abundance ===\n")
  cat("  Rank:", rank, "\n")

  # Agglomerate to chosen rank
  ps_agg <- tax_glom(ps, taxrank = rank, NArm = FALSE)

  # Rename taxa using rank label
  new_names <- as.character(tax_table(ps_agg)[, rank])
  new_names[is.na(new_names)] <- paste0("Unknown_", seq_len(sum(is.na(new_names))))
  # Handle duplicate names
  dup_idx    <- duplicated(new_names)
  new_names[dup_idx] <- paste0(new_names[dup_idx], "_dup",
                                seq_len(sum(dup_idx)))
  taxa_names(ps_agg) <- new_names

  n_before <- ntaxa(ps_agg)

  # Prevalence filter
  otu_mat   <- as.matrix(otu_table(ps_agg))
  if (!taxa_are_rows(ps_agg)) otu_mat <- t(otu_mat)

  prevalence <- rowSums(otu_mat > 0) / ncol(otu_mat)
  total_count <- rowSums(otu_mat)

  keep       <- prevalence >= min_prevalence & total_count >= min_count
  ps_filt    <- prune_taxa(keep, ps_agg)

  cat("  Taxa before filtering:", n_before, "\n")
  cat("  Taxa after filtering: ", ntaxa(ps_filt), "\n")
  cat("  Removed:              ", n_before - ntaxa(ps_filt), "\n\n")

  return(ps_filt)
}


# =============================================================================
# SECTION 2 — ANCOM-BC
# =============================================================================

#' Run ANCOM-BC differential abundance analysis.
#'
#' ANCOM-BC accounts for compositional bias and provides bias-corrected
#' log fold changes. Recommended as the primary DA method.
#'
#' @param ps          A phyloseq object (raw counts).
#' @param formula     Model formula as string. E.g. "group" or "group + age + sex".
#' @param group_var   The primary grouping variable (used for result extraction).
#' @param reference   Reference group level. Default = first alphabetically.
#' @param p_adj_method P-value adjustment method. Default = "BH".
#' @param alpha       Significance threshold. Default = 0.05.
#' @return A tidy data frame of ANCOM-BC results.

run_ancombc <- function(ps,
                         formula    = "group",
                         group_var  = "group",
                         reference  = NULL,
                         p_adj_method = "BH",
                         alpha      = 0.05) {

  cat("=== ANCOM-BC ===\n")

  if (!pkg_available("ANCOMBC")) {
    stop("ANCOMBC package not installed. Run: BiocManager::install('ANCOMBC')")
  }
  library(ANCOMBC)

  # Set reference level
  if (!is.null(reference)) {
    sample_data(ps)[[group_var]] <- relevel(
      factor(sample_data(ps)[[group_var]]),
      ref = reference
    )
  }

  cat("  Running ANCOM-BC2...\n")
  set.seed(42)

  ancombc_res <- tryCatch({
    withCallingHandlers(
      ancombc2(
        data         = ps,
        tax_level    = NULL,
        fix_formula  = formula,
        rand_formula = NULL,
        p_adj_method = p_adj_method,
        pseudo_sens  = TRUE,
        prv_cut      = 0,
        lib_cut      = 0,
        s0_perc      = 0.05,
        group        = group_var,
        struc_zero   = TRUE,
        neg_lb       = TRUE,
        alpha        = alpha,
        n_cl         = 1,
        verbose      = FALSE
      ),
      warning = function(w) invokeRestart("muffleWarning")
    )
  }, error = function(e) {
    cat("  ANCOM-BC2 failed:", conditionMessage(e), "\n")
    cat("  Trying ANCOM-BC1...\n")
    ancombc(
      phyloseq     = ps,
      formula      = formula,
      p_adj_method = p_adj_method,
      zero_cut     = 0.90,
      lib_cut      = 0,
      group        = group_var,
      struc_zero   = TRUE,
      neg_lb       = TRUE,
      tol          = 1e-5,
      max_iter     = 100,
      conserve     = TRUE,
      alpha        = alpha,
      global       = FALSE
    )
  })

  # --- Extract results -------------------------------------------------------
  if (inherits(ancombc_res, "ANCOMBC2") || (is.list(ancombc_res) && "res" %in% names(ancombc_res) && "taxon" %in% colnames(ancombc_res$res))) {
    res_df <- ancombc_res$res
    lfc_cols  <- grep("^lfc_(?!\\(Intercept\\))", colnames(res_df), value = TRUE, perl = TRUE)
    q_cols    <- grep("^q_(?!\\(Intercept\\))",   colnames(res_df), value = TRUE, perl = TRUE)
    diff_cols <- grep("^diff_(?!\\(Intercept\\)|robust)", colnames(res_df), value = TRUE, perl = TRUE)
    
    n <- min(length(lfc_cols), length(q_cols), length(diff_cols))
    if (n == 0) stop("ANCOMBC2 result has no comparison columns to extract.")
    
    tidy_list <- lapply(seq_len(n), function(i) {
      data.frame(
        taxon      = res_df$taxon,
        comparison = sub("^lfc_", "", lfc_cols[i]),
        lfc        = res_df[[lfc_cols[i]]],
        q_value    = res_df[[q_cols[i]]],
        diff_abund = res_df[[diff_cols[i]]],
        method     = "ANCOM-BC2",
        stringsAsFactors = FALSE
      )
    })
    result_df <- dplyr::bind_rows(tidy_list)
  } else {
    # ANCOM-BC1 format
    res        <- ancombc_res$res
    result_df  <- data.frame(
      taxon      = rownames(res$beta),
      comparison = group_var,
      lfc        = as.numeric(res$beta[, 1]),
      se         = as.numeric(res$se[, 1]),
      q_value    = as.numeric(res$q_val[, 1]),
      diff_abund = as.logical(res$diff_abn[, 1]),
      method     = "ANCOM-BC",
      stringsAsFactors = FALSE
    )
  }

  n_sig <- sum(result_df$diff_abund, na.rm = TRUE)
  cat("  Differentially abundant taxa:", n_sig, "\n\n")

  return(result_df)
}


# =============================================================================
# SECTION 3 — DESeq2
# =============================================================================

#' Run DESeq2 differential abundance analysis.
#'
#' DESeq2 was designed for RNA-seq but is widely used for microbiome data.
#' Includes a zero-inflation pseudocount strategy appropriate for sparse data.
#'
#' @param ps          A phyloseq object (raw counts).
#' @param group_var   Grouping variable.
#' @param reference   Reference group. Default = first alphabetically.
#' @param alpha       Significance threshold. Default = 0.05.
#' @param lfc_threshold Minimum absolute log2 fold change. Default = 1.
#' @return A tidy data frame of DESeq2 results.

run_deseq2 <- function(ps,
                       group_var     = "group",
                       reference     = NULL,
                       alpha         = 0.05,
                       lfc_threshold = 1) {
  
  cat("=== DESeq2 ===\n")
  
  if (!pkg_available("DESeq2")) {
    stop("DESeq2 not installed. Run: BiocManager::install('DESeq2')")
  }
  
  library(DESeq2)
  
  # Ensure grouping variable is a factor
  sample_data(ps)[[group_var]] <- factor(sample_data(ps)[[group_var]])
  
  # Set reference level
  if (!is.null(reference)) {
    sample_data(ps)[[group_var]] <- relevel(
      sample_data(ps)[[group_var]],
      ref = reference
    )
  }
  
  # Add pseudocount of 1 to handle zeros
  ps_counts <- ps
  otu_table(ps_counts) <- otu_table(
    round(as.matrix(otu_table(ps_counts)) + 1),
    taxa_are_rows = taxa_are_rows(ps_counts)
  )
  
  cat("  Converting to DESeqDataSet...\n")
  
  dds <- phyloseq_to_deseq2(
    ps_counts,
    as.formula(paste("~", group_var))
  )
  
  dds <- estimateSizeFactors(dds, type = "poscounts")
  
  cat("  Running DESeq2...\n")
  
  set.seed(42)
  dds <- DESeq(
    dds,
    test = "Wald",
    fitType = "local",
    quiet = TRUE
  )
  
  group_levels <- levels(sample_data(ps_counts)[[group_var]])
  ref_level <- if (!is.null(reference)) reference else group_levels[1]
  other_levels <- setdiff(group_levels, ref_level)
  
  result_list <- lapply(other_levels, function(lvl) {
    
    res <- results(
      dds,
      contrast = c(group_var, lvl, ref_level),
      alpha = alpha,
      lfcThreshold = lfc_threshold,
      pAdjustMethod = "BH"
    )
    
    data.frame(
      taxon = rownames(res),
      baseMean = res$baseMean,
      lfc = res$log2FoldChange,
      lfcSE = res$lfcSE,
      stat = res$stat,
      pvalue = res$pvalue,
      q_value = res$padj,
      comparison = paste0(lvl, "_vs_", ref_level),
      diff_abund = !is.na(res$padj) &
        res$padj < alpha &
        abs(res$log2FoldChange) >= lfc_threshold,
      method = "DESeq2",
      stringsAsFactors = FALSE
    )
  })
  
  result_df <- dplyr::bind_rows(result_list)
  
  n_sig <- sum(result_df$diff_abund, na.rm = TRUE)
  cat("  Differentially abundant taxa:", n_sig, "\n\n")
  
  return(result_df)
}

# =============================================================================
# SECTION 4 — ALDEx2
# =============================================================================

#' Run ALDEx2 differential abundance analysis.
#'
#' ALDEx2 uses a Dirichlet-multinomial model with CLR transformation,
#' making it robust to compositional effects.
#'
#' @param ps          A phyloseq object (raw counts).
#' @param group_var   Grouping variable (currently supports 2 groups only).
#' @param mc_samples  Number of Monte Carlo samples. Default = 128.
#' @param alpha       Significance threshold. Default = 0.05.
#' @return A tidy data frame of ALDEx2 results.

run_aldex2 <- function(ps,
                        group_var  = "group",
                        mc_samples = 128,
                        alpha      = 0.05) {

  cat("=== ALDEx2 ===\n")

  if (!pkg_available("ALDEx2")) {
    cat("  ALDEx2 not installed. Skipping. ",
        "Run: BiocManager::install('ALDEx2')\n\n")
    return(NULL)
  }
  library(ALDEx2)

  groups      <- as.character(sample_data(ps)[[group_var]])
  group_levels <- unique(groups)

  if (length(group_levels) != 2) {
    cat("  ALDEx2 currently supports 2-group comparisons only.",
        "Skipping for", length(group_levels), "groups.\n\n")
    return(NULL)
  }

  otu_mat <- as.matrix(otu_table(ps))
  if (taxa_are_rows(ps)) otu_mat <- otu_mat else otu_mat <- t(otu_mat)

  cat("  Running ALDEx2 with", mc_samples, "MC samples...\n")
  set.seed(42)

  aldex_res <- tryCatch({
    clr_data <- aldex.clr(otu_mat, conds = groups,
                           mc.samples = mc_samples,
                           denom = "all", verbose = FALSE)
    t_test   <- aldex.ttest(clr_data, paired.test = FALSE, verbose = FALSE)
    effect   <- aldex.effect(clr_data, CI = TRUE, verbose = FALSE)
#    aldex.plot(clr_data, type = "MA", test = "welch")   # suppress but run
    cbind(t_test, effect)
  }, error = function(e) {
    cat("  ALDEx2 error:", e$message, "\n")
    return(NULL)
  })

  if (is.null(aldex_res)) return(NULL)

  result_df <- aldex_res %>%
    rownames_to_column("taxon") %>%
    mutate(
      comparison = paste0(group_levels[2], "_vs_", group_levels[1]),
      lfc        = effect,
      q_value    = wi.eBH,        # Welch's t-test BH-adjusted p
      diff_abund = wi.eBH < alpha & abs(effect) > 1,
      method     = "ALDEx2"
    ) %>%
    dplyr::select(taxon, comparison, lfc, q_value, diff_abund, method,
                  dplyr::everything())

  n_sig <- sum(result_df$diff_abund, na.rm = TRUE)
  cat("  Differentially abundant taxa:", n_sig, "\n\n")

  return(result_df)
}


# =============================================================================
# SECTION 5 — CONSENSUS SCORING
# =============================================================================

#' Build a consensus differential abundance table across multiple methods.
#'
#' Taxa are scored by how many methods identify them as significantly
#' differentially abundant. Consensus taxa are more robust than
#' single-method results.
#'
#' @param result_list Named list of result data frames from individual methods.
#' @param alpha       Significance threshold. Default = 0.05.
#' @param min_methods Minimum methods that must agree. Default = 2.
#' @return A tidy data frame with consensus scores and direction agreement.

build_consensus <- function(result_list,
                            alpha       = 0.05,
                            min_methods = 2) {
  
  cat("=== Building multi-method consensus ===\n")
  
  result_list <- Filter(Negate(is.null), result_list)
  
  if (length(result_list) == 0) {
    cat("  No differential abundance methods available.\n\n")
    
    empty_df <- data.frame(
      taxon = character(),
      comparison = character(),
      lfc = numeric(),
      q_value = numeric(),
      diff_abund = logical(),
      method = character(),
      stringsAsFactors = FALSE
    )
    
    return(list(
      consensus = empty_df,
      all_results = empty_df
    ))
  }
  
  if (length(result_list) == 1) {
    cat("  Only 1 method available. Using single-method results.\n\n")
    
    single <- dplyr::bind_rows(result_list) %>%
      dplyr::mutate(
        n_methods_total = 1,
        n_methods_sig = as.integer(diff_abund),
        consensus_da = diff_abund,
        direction = dplyr::case_when(
          diff_abund & lfc > 0 ~ "increased",
          diff_abund & lfc < 0 ~ "decreased",
          TRUE ~ "not_significant"
        ),
        confidence = ifelse(diff_abund, "single_method", "low"),
        mean_lfc = lfc,
        median_lfc = lfc,
        min_q = q_value,
        methods_detected = ifelse(diff_abund, method, "")
      )
    
    return(list(
      consensus = single,
      all_results = single
    ))
  }
  
  std_list <- lapply(result_list, function(df) {
    df %>%
      dplyr::select(taxon, comparison, lfc, q_value, diff_abund, method) %>%
      dplyr::mutate(
        direction = dplyr::case_when(
          diff_abund & lfc > 0 ~ "increased",
          diff_abund & lfc < 0 ~ "decreased",
          TRUE ~ "not_significant"
        )
      )
  })
  
  combined <- dplyr::bind_rows(std_list)
  
  consensus <- combined %>%
    dplyr::group_by(taxon, comparison) %>%
    dplyr::summarise(
      n_methods_total  = dplyr::n(),
      n_methods_sig    = sum(diff_abund, na.rm = TRUE),
      n_increased      = sum(direction == "increased", na.rm = TRUE),
      n_decreased      = sum(direction == "decreased", na.rm = TRUE),
      mean_lfc         = round(mean(lfc, na.rm = TRUE), 4),
      median_lfc       = round(median(lfc, na.rm = TRUE), 4),
      min_q            = round(min(q_value, na.rm = TRUE), 6),
      methods_detected = paste(method[diff_abund], collapse = "; "),
      .groups          = "drop"
    ) %>%
    dplyr::mutate(
      consensus_da = n_methods_sig >= min_methods,
      direction = dplyr::case_when(
        n_increased > n_decreased ~ "increased",
        n_decreased > n_increased ~ "decreased",
        TRUE ~ "conflicting"
      ),
      confidence = dplyr::case_when(
        n_methods_sig == n_methods_total ~ "high",
        n_methods_sig >= min_methods ~ "moderate",
        TRUE ~ "low"
      )
    ) %>%
    dplyr::arrange(dplyr::desc(n_methods_sig), min_q)
  
  n_consensus <- sum(consensus$consensus_da)
  
  cat("  Methods used:", length(result_list), "\n")
  cat("  Consensus DA taxa (>=", min_methods, "methods):", n_consensus, "\n")
  cat("    Increased:", sum(consensus$consensus_da & consensus$direction == "increased"), "\n")
  cat("    Decreased:", sum(consensus$consensus_da & consensus$direction == "decreased"), "\n\n")
  
  return(list(
    consensus = consensus,
    all_results = combined
  ))
}

# =============================================================================
# SECTION 6 — VOLCANO PLOT
# =============================================================================

#' Generate a volcano plot from differential abundance results.
#'
#' @param result_df   Data frame with columns: taxon, lfc, q_value, diff_abund.
#' @param lfc_col     Column name for log fold change. Default = "lfc".
#' @param q_col       Column name for adjusted p-value. Default = "q_value".
#' @param alpha       Significance threshold. Default = 0.05.
#' @param lfc_threshold Fold change threshold for labelling. Default = 1.
#' @param top_n_label Number of top taxa to label. Default = 15.
#' @param title       Plot title.
#' @return A ggplot object.

plot_volcano <- function(result_df,
                          lfc_col       = "lfc",
                          q_col         = "q_value",
                          alpha         = 0.05,
                          lfc_threshold = 1,
                          top_n_label   = 15,
                          title         = "Differential Abundance") {

  cat("=== Volcano plot ===\n")

  df <- result_df %>%
    filter(!is.na(.data[[lfc_col]]), !is.na(.data[[q_col]])) %>%
    mutate(
      neg_log_q   = -log10(.data[[q_col]] + 1e-10),
      status      = case_when(
        .data[[q_col]] < alpha & .data[[lfc_col]] >  lfc_threshold ~ "Increased",
        .data[[q_col]] < alpha & .data[[lfc_col]] < -lfc_threshold ~ "Decreased",
        TRUE                                                         ~ "Not significant"
      )
    )

  # Top taxa to label (by q-value, among significant)
  top_label <- df %>%
    filter(status != "Not significant") %>%
    arrange(.data[[q_col]]) %>%
    slice_head(n = top_n_label)

  colour_map <- c(
    "Increased"       = "#e74c3c",
    "Decreased"       = "#3498db",
    "Not significant" = "#bdc3c7"
  )

  n_up   <- sum(df$status == "Increased")
  n_down <- sum(df$status == "Decreased")

  p <- ggplot(df, aes(x = .data[[lfc_col]], y = neg_log_q,
                       colour = status, size = status)) +
    # Background points
    geom_point(alpha = 0.6) +
    # Significance thresholds
    geom_hline(yintercept = -log10(alpha),
               linetype = "dashed", colour = "grey50", linewidth = 0.7) +
    geom_vline(xintercept =  lfc_threshold,
               linetype = "dashed", colour = "grey50", linewidth = 0.7) +
    geom_vline(xintercept = -lfc_threshold,
               linetype = "dashed", colour = "grey50", linewidth = 0.7) +
    # Labels for top taxa
    ggrepel::geom_text_repel(
      data          = top_label,
      aes(label     = taxon),
      size          = 2.8,
      fontface      = "italic",
      max.overlaps  = 20,
      segment.size  = 0.3,
      segment.alpha = 0.5,
      box.padding   = 0.3
    ) +
    # Annotations
    annotate("text", x =  max(df[[lfc_col]], na.rm=TRUE) * 0.85,
             y = max(df$neg_log_q) * 0.95,
             label = paste0("Up ", n_up, " increased"),
             colour = "#e74c3c", size = 3.5, fontface = "bold") +
    annotate("text", x = min(df[[lfc_col]], na.rm=TRUE) * 0.85,
             y = max(df$neg_log_q) * 0.95,
             label = paste0("Down ", n_down, " decreased"),
             colour = "#3498db", size = 3.5, fontface = "bold") +
    scale_colour_manual(values = colour_map, name = "Status") +
    scale_size_manual(values = c("Increased" = 2.5,
                                  "Decreased" = 2.5,
                                  "Not significant" = 1.2),
                      guide = "none") +
    labs(
      title    = title,
      subtitle = paste0("Significance: q < ", alpha,
                        " | LFC threshold: ±", lfc_threshold),
      x        = "Log fold change",
      y        = expression(-log[10](q-value))
    ) +
    theme_microbiome()

  return(p)
}


# =============================================================================
# SECTION 7 — EFFECT SIZE PLOT
# =============================================================================

#' Plot ranked effect sizes (log fold changes) for significant taxa.
#'
#' @param result_df       Data frame with taxon, lfc, q_value, diff_abund.
#' @param lfc_col         Log fold change column. Default = "lfc".
#' @param q_col           Adjusted p-value column. Default = "q_value".
#' @param alpha           Significance threshold. Default = 0.05.
#' @param se_col          Standard error column for error bars (optional).
#' @param top_n           Maximum number of taxa to display. Default = 30.
#' @param group_comparison Title annotation for the comparison shown.
#' @return A ggplot object.

plot_effect_sizes <- function(result_df,
                               lfc_threshold    = 1,
                               lfc_col          = "lfc",
                               q_col            = "q_value",
                               alpha            = 0.05,
                               se_col           = NULL,
                               top_n            = 30,
                               group_comparison = "") {

  cat("=== Effect size plot ===\n")

  df_sig <- result_df %>%
    filter(!is.na(.data[[q_col]]),
           .data[[q_col]] < alpha) %>%
    arrange(desc(abs(.data[[lfc_col]]))) %>%
    slice_head(n = top_n) %>%
    mutate(
      taxon     = fct_reorder(taxon, .data[[lfc_col]]),
      direction = ifelse(.data[[lfc_col]] > 0, "Increased", "Decreased")
    )

  if (nrow(df_sig) == 0) {
    cat("  No significant taxa to plot.\n\n")
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = "No significant taxa",
                      size = 6, colour = "grey60") +
             theme_void())
  }

  cat("  Plotting", nrow(df_sig), "significant taxa\n\n")

  colour_map <- c("Increased" = "#e74c3c", "Decreased" = "#3498db")

  p <- ggplot(df_sig, aes(x = .data[[lfc_col]], y = taxon,
                            colour = direction, fill = direction)) +
    geom_col(width = 0.7, alpha = 0.85) +
    {
      if (!is.null(se_col) && se_col %in% colnames(df_sig)) {
        geom_errorbarh(
          aes(xmin = .data[[lfc_col]] - 1.96 * .data[[se_col]],
              xmax = .data[[lfc_col]] + 1.96 * .data[[se_col]]),
          height = 0.3, linewidth = 0.6, colour = "grey30"
        )
      }
    } +
    geom_vline(xintercept = 0, linewidth = 0.8, colour = "grey30") +
    # q-value labels
    geom_text(
      aes(x     = ifelse(.data[[lfc_col]] > 0,
                         .data[[lfc_col]] + 0.02 * diff(range(.data[[lfc_col]])),
                         .data[[lfc_col]] - 0.02 * diff(range(.data[[lfc_col]]))),
          label = paste0("q=", formatC(.data[[q_col]], format = "e", digits = 1))),
      size  = 2.2, colour = "grey40",
      hjust = ifelse(df_sig[[lfc_col]] > 0, 0, 1)
    ) +
    scale_colour_manual(values = colour_map, guide = "none") +
    scale_fill_manual(values = colour_map, name = "Direction") +
    labs(
      title    = paste0("Differential abundance - ", group_comparison),
      subtitle = paste0(nrow(df_sig), " significant taxa (q < ", alpha, ")"),
      x        = "Log fold change",
      y        = NULL
    ) +
    theme_microbiome() +
    theme(axis.text.y = element_text(face = "italic", size = 8))

  return(p)
}


# =============================================================================
# SECTION 8 — DA HEATMAP
# =============================================================================

#' Heatmap of significant taxa abundance across all samples.
#'
#' @param ps          A phyloseq object.
#' @param da_taxa     Character vector of significant taxon names.
#' @param group_var   Grouping variable for column annotation.
#' @param transform   Transformation: "clr", "log10", or "zscore". Default = "zscore".
#' @param rank        Taxonomic rank (must match taxa names in ps). Default = "Genus".
#' @return A ggplot heatmap.

plot_da_heatmap <- function(ps,
                             da_taxa,
                             group_var   = "group",
                             transform   = "zscore",
                             rank        = "Genus",
                             sort_by_group   = TRUE,
                             cluster_cols    = TRUE) {

  # Agglomerate and filter to DA taxa
  ps_agg  <- tax_glom(ps, taxrank = rank, NArm = FALSE)
  taxa_names(ps_agg) <- as.character(tax_table(ps_agg)[, rank])
  shared  <- intersect(da_taxa, taxa_names(ps_agg))

  if (length(shared) == 0) {
    cat("  No matching DA taxa found in phyloseq object.\n\n")
    return(NULL)
  }

  ps_da   <- prune_taxa(shared, ps_agg)
  cat("  Plotting", ntaxa(ps_da), "DA taxa\n")

  otu_mat <- as.matrix(otu_table(ps_da))
  if (!taxa_are_rows(ps_da)) otu_mat <- t(otu_mat)

  # Transformation
  if (transform == "clr") {
    otu_plot <- log(otu_mat + 0.5) - rowMeans(log(otu_mat + 0.5))
  } else if (transform == "log10") {
    otu_plot <- log10(otu_mat + 1)
  } else {   # z-score per taxon
    otu_plot <- t(scale(t(otu_mat)))
  }

  # Hierarchical clustering of taxa (rows)
  row_clust <- hclust(dist(otu_plot), method = "ward.D2")
  row_order <- rownames(otu_plot)[row_clust$order]
  # Column ordering
  cat("  sort_by_group:", sort_by_group, "| cluster_cols:", cluster_cols, "\n")
  if (sort_by_group && !is.null(group_var)) {
    meta_df_ord <- data.frame(sample_data(ps_da))
    meta_df_ord <- data.frame(sample_data(ps_da))
    groups      <- unique(as.character(meta_df_ord[[group_var]]))
    col_order   <- unlist(lapply(groups, function(g) {
      samps <- rownames(meta_df_ord)[meta_df_ord[[group_var]] == g]
      samps <- intersect(samps, colnames(otu_plot))
      if (cluster_cols && length(samps) > 2) {
        samps[hclust(dist(t(otu_plot[, samps, drop=FALSE])), method="ward.D2")$order]
      } else { samps }
    }))
  } else {
    col_order <- colnames(otu_plot)[hclust(dist(t(otu_plot)), method="ward.D2")$order]
  }

  # Melt
  heat_df <- as.data.frame(otu_plot) %>%
    rownames_to_column("taxon") %>%
    pivot_longer(-taxon, names_to = "sample", values_to = "value") %>%
    mutate(
      taxon  = factor(taxon, levels = rev(row_order)),
      sample = factor(sample, levels = col_order)
    )

  # Add group annotation
  meta_df <- data.frame(sample_data(ps)) %>%
    dplyr::select(all_of(group_var)) %>%
    rownames_to_column("sample")
  heat_df <- left_join(heat_df, meta_df, by = "sample")

  n_groups <- n_distinct(heat_df[[group_var]])
  colours  <- brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)]
  group_levels <- unique(heat_df[[group_var]])
  names(colours) <- group_levels

  # Annotation strip
  annot_df <- heat_df %>% dplyr::select(sample, all_of(group_var)) %>% distinct()

  p_annot <- ggplot(annot_df, aes(x = sample, y = 1,
                                   fill = .data[[group_var]])) +
    geom_tile() +
    scale_fill_manual(values = colours, name = group_var) +
    scale_x_discrete(expand = c(0, 0)) +
    theme_void() +
    theme(legend.position = "right")

  # Heatmap
  fill_label <- switch(transform,
                        clr = "CLR",
                        log10 = "Log10",
                        zscore = "Z-score")

  p_heat <- ggplot(heat_df, aes(x = sample, y = taxon, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.1) +
    scale_fill_gradient2(
      low      = "#3498db",
      mid      = "white",
      high     = "#e74c3c",
      midpoint = 0,
      name     = fill_label
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "Samples", y = NULL) +
    theme_microbiome() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(face = "italic", size = 8),
      panel.border = element_rect(colour = "grey80", fill = NA)
    )

  combined <- p_annot / p_heat +
    plot_layout(heights = c(0.06, 1), guides = "collect") +
    plot_annotation(
      title    = paste0("Differentially abundant taxa - ", rank, " level"),
      subtitle = paste0(ntaxa(ps_da), " taxa | Transform: ", transform),
      theme    = theme(plot.title = element_text(size = 13, face = "bold"))
    )

  cat("\n")
  return(combined)
}


# =============================================================================
# SECTION 9 — BOXPLOTS FOR TOP DA TAXA
# =============================================================================

#' Boxplots of individual taxon abundance for the top DA taxa.
#'
#' @param ps          A phyloseq object.
#' @param da_taxa     Named character vector of significant taxa (name = taxon, value = direction).
#' @param group_var   Grouping variable.
#' @param rank        Taxonomic rank. Default = "Genus".
#' @param top_n       Number of taxa to plot. Default = 12.
#' @param transform   "relative" or "log10". Default = "relative".
#' @return A patchwork of individual taxon boxplots.

plot_da_boxplots <- function(ps,
                              da_taxa,
                              group_var = "group",
                              rank      = "Genus",
                              top_n     = 12,
                              transform = "relative") {

  cat("=== DA boxplots ===\n")

  # Prepare data
  ps_rel  <- transform_sample_counts(ps, function(x) x / sum(x))
  ps_agg  <- tax_glom(ps_rel, taxrank = rank, NArm = FALSE)
  taxa_names(ps_agg) <- as.character(tax_table(ps_agg)[, rank])

  # Keep top DA taxa
  top_taxa <- head(da_taxa, top_n)
  shared   <- intersect(top_taxa, taxa_names(ps_agg))
  ps_sub   <- prune_taxa(shared, ps_agg)

  ps_melt <- psmelt(ps_sub) %>%
    rename(taxon = OTU) %>%
    mutate(
      Abundance = if (transform == "log10") log10(Abundance + 1e-4) else Abundance,
      taxon     = factor(taxon, levels = shared)
    )

  groups   <- unique(as.character(ps_melt[[group_var]]))
  n_groups <- length(groups)
  colours  <- brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)]
  names(colours) <- groups

  y_label <- if (transform == "log10") "Log10(relative abundance)" else "Relative abundance"

  # Individual taxon plots
  plots <- lapply(shared, function(tx) {
    df_sub <- ps_melt %>% filter(taxon == tx)

    # Wilcoxon or Kruskal
    stat_res <- tryCatch({
      if (n_groups == 2) {
        wt <- wilcox.test(Abundance ~ .data[[group_var]], data = df_sub)
        p_lab <- paste0("p=", signif(wt$p.value, 2))
      } else {
        kt <- kruskal.test(Abundance ~ .data[[group_var]], data = df_sub)
        p_lab <- paste0("KW p=", signif(kt$p.value, 2))
      }
      p_lab
    }, error = function(e) "")

    ggplot(df_sub, aes(x = .data[[group_var]], y = Abundance,
                       fill = .data[[group_var]])) +
      geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.55,
                   linewidth = 0.4) +
      geom_jitter(aes(colour = .data[[group_var]]),
                  width = 0.15, alpha = 0.5, size = 0.8) +
      scale_fill_manual(values = colours, guide = "none") +
      scale_colour_manual(values = colours, guide = "none") +
      annotate("text", x = Inf, y = Inf, vjust = 1.5, hjust = 1.2,
               label = stat_res, size = 2.5, colour = "grey40") +
      labs(
        title = tx,
        x = NULL,
        y = if (tx == shared[1]) y_label else NULL
      ) +
      theme_microbiome() +
      theme(
        plot.title   = element_text(size = 8, face = "italic", hjust = 0.5),
        axis.text.x  = element_text(angle = 30, hjust = 1, size = 7),
        axis.title.y = element_text(size = 8)
      )
  })

  n_cols <- min(4, length(plots))
  combined <- wrap_plots(plots, ncol = n_cols) +
    plot_annotation(
      title    = paste0("Top ", length(plots), " differential taxa - ", rank),
      subtitle = paste0(y_label, " | Grouped by ", group_var),
      theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                       plot.subtitle = element_text(size = 10, colour = "grey40"))
    )

  cat("  Boxplots produced for", length(plots), "taxa\n\n")
  return(combined)
}


# =============================================================================
# SECTION 10 — COMPLETE DIFFERENTIAL ABUNDANCE WORKFLOW WRAPPER
# =============================================================================

#' Run the complete differential abundance pipeline with multiple methods.
#'
#' @param ps          A filtered phyloseq object (raw counts).
#' @param group_var   Primary grouping variable.
#' @param formula     Full model formula for ANCOM-BC. Default = group_var.
#' @param reference   Reference group. Default = first level.
#' @param rank        Taxonomic rank. Default = "Genus".
#' @param alpha       Significance threshold. Default = 0.05.
#' @param lfc_threshold Log fold change threshold. Default = 1.
#' @param min_methods Minimum methods for consensus. Default = 2.
#' @param output_dir  Directory to save all outputs.
#' @return A named list of all results.

run_differential_abundance <- function(ps,
                                        group_var     = NULL,
                                        formula       = NULL,
                                        reference     = NULL,
                                        rank          = "Genus",
                                        alpha         = 0.05,
                                        lfc_threshold = 1,
                                        min_methods   = 2,
                                        output_dir    = "da_output") {

  cat("\n", strrep("=", 60), "\n")
  cat("  DIFFERENTIAL ABUNDANCE PIPELINE\n")
  cat(strrep("=", 60), "\n\n")

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (is.null(formula)) formula <- group_var

  results <- list()

  # --- Prepare data ---------------------------------------------------------
  ps_da <- prepare_da_data(ps, rank = rank,
                            min_prevalence = 0.10, min_count = 10)

  # --- Run each method ------------------------------------------------------
#  results$ancombc <- tryCatch(
#    run_ancombc(ps_da, formula = formula, group_var = group_var,
#                reference = reference, alpha = alpha),
#    error = function(e) { cat("ANCOM-BC failed:", e$message, "\n"); NULL }
#  )

  results$ancombc <- NULL
  cat("ANCOM-BC skipped: unavailable on this system.\n")
  
  results$deseq2 <- tryCatch(
    run_deseq2(ps_da, group_var = group_var, reference = reference,
               alpha = alpha, lfc_threshold = lfc_threshold),
    error = function(e) { cat("DESeq2 failed:", e$message, "\n"); NULL }
  )

  results$aldex2 <- tryCatch(
    run_aldex2(ps_da, group_var = group_var, alpha = alpha),
    error = function(e) { cat("ALDEx2 failed:", e$message, "\n"); NULL }
  )

  # --- Consensus -----------------------------------------------------------
  cat("--- Building consensus ---\n")
  consensus_res <- build_consensus(
    result_list = list(
      ancombc = results$ancombc,
      deseq2  = results$deseq2,
      aldex2  = results$aldex2
    ),
    alpha       = alpha,
    min_methods = min_methods
  )
  results$consensus <- consensus_res

  # Save all results tables
  write.csv(consensus_res$consensus,
            file.path(output_dir, "da_consensus_results.csv"),
            row.names = FALSE)
  write.csv(consensus_res$all_results,
            file.path(output_dir, "da_all_methods_results.csv"),
            row.names = FALSE)

  # --- Significant taxa from consensus ------------------------------------
  sig_taxa <- consensus_res$consensus %>%
    filter(consensus_da) %>%
    arrange(desc(n_methods_sig), min_q) %>%
    pull(taxon)

  cat("  Consensus DA taxa:", length(sig_taxa), "\n\n")

  # --- Plots ----------------------------------------------------------------
  # Use ANCOM-BC results for volcano/effect size (most recommended)
  primary_res <- if (!is.null(results$ancombc)) results$ancombc else
                 if (!is.null(results$deseq2))  results$deseq2 else
                 consensus_res$all_results

  cat("--- Plot 1: Volcano plot ---\n")
  results$p_volcano <- plot_volcano(
    primary_res,
    alpha         = alpha,
    lfc_threshold = lfc_threshold,
    top_n_label   = 20,
    title         = paste0("Differential Abundance - ", group_var)
  )
  ggsave(file.path(output_dir, "01_volcano.pdf"),
         results$p_volcano, width = 10, height = 8)

  cat("--- Plot 2: Effect size plot ---\n")
  results$p_effects <- plot_effect_sizes(
    primary_res,
    alpha            = alpha,
    top_n            = 30,
    se_col           = if ("se" %in% colnames(primary_res)) "se" else NULL,
    group_comparison = group_var
  )
  
  ggsave(file.path(output_dir, "02_effect_sizes.pdf"),
         results$p_effects, width = 10, height = max(6, length(sig_taxa) * 0.35))

  if (length(sig_taxa) > 0) {
    cat("--- Plot 3: DA heatmap ---\n")
    results$p_heatmap <- plot_da_heatmap(ps, sig_taxa,
                                          group_var = group_var,
                                          rank = rank, transform = "zscore")
    if (!is.null(results$p_heatmap)) {
      ggsave(file.path(output_dir, "03_da_heatmap.pdf"),
             results$p_heatmap,
             width = max(10, nsamples(ps) * 0.18),
             height = max(6, length(sig_taxa) * 0.35))
    }

    cat("--- Plot 4: DA boxplots ---\n")
    results$p_boxplots <- plot_da_boxplots(
      ps, sig_taxa, group_var = group_var,
      rank = rank, top_n = 12
    )
    ggsave(file.path(output_dir, "04_da_boxplots.pdf"),
           results$p_boxplots, width = 16, height = 12)
  }

  # --- Multi-method comparison Venn -----------------------------------------
  cat("--- Table 2: Method comparison ---\n")
  method_comparison <- consensus_res$all_results %>%
    filter(diff_abund) %>%
    group_by(method) %>%
    summarise(n_significant = n(), .groups = "drop")
  print(method_comparison)
  write.csv(method_comparison,
            file.path(output_dir, "method_comparison.csv"),
            row.names = FALSE)

  cat("\n", strrep("=", 60), "\n")
  cat("  DIFFERENTIAL ABUNDANCE PIPELINE COMPLETE\n")
  cat(strrep("=", 60), "\n")
  cat("  Output saved to:", output_dir, "\n")
  cat("  Tables: consensus results, all methods, method comparison\n")
  cat("  Plots:  up to 4 PDF files\n")
  cat("  Consensus DA taxa:", length(sig_taxa), "\n\n")

  return(invisible(results))
}


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# ps <- readRDS("qc_output/phyloseq_qc_filtered.rds")

# --- Option A: Full pipeline with all methods --------------------------------
# results <- run_differential_abundance(
#   ps            = ps,
#   group_var     = "disease_status",
#   formula       = "disease_status + age + sex",
#   reference     = "Healthy",
#   rank          = "Genus",
#   alpha         = 0.05,
#   lfc_threshold = 1,
#   min_methods   = 2,
#   output_dir    = "results/differential_abundance"
# )

# --- Option B: Single method -------------------------------------------------
# ps_da     <- prepare_da_data(ps, rank = "Genus")
# ancom_res <- run_ancombc(ps_da, formula = "disease_status",
#                           group_var = "disease_status", reference = "Healthy")
# plot_volcano(ancom_res)
# plot_effect_sizes(ancom_res, top_n = 25)

# --- Option C: Access consensus results -------------------------------------
# results$consensus$consensus %>%
#   filter(consensus_da) %>%
#   select(taxon, direction, confidence, n_methods_sig, mean_lfc, min_q)
