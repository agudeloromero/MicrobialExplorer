# =============================================================================
# MICROBIOME DATA IMPORT AND QUALITY CONTROL
# =============================================================================
# Description : Comprehensive QC pipeline for 16S microbiome data
# Input       : OTU/ASV table, taxonomy table, metadata, optional phylogenetic tree
# Output      : Filtered phyloseq object, QC report, diagnostic plots
# Author      : Patricia
# Dependencies: phyloseq, vegan, ggplot2, dplyr, tidyr, scales, patchwork,
#               microbiome, decontam, knitr
# =============================================================================

# --- 1. LOAD LIBRARIES -------------------------------------------------------

suppressPackageStartupMessages({
  library(phyloseq)       # Core microbiome data structure
  library(vegan)          # Ecological diversity and ordination
  library(ggplot2)        # Visualisation
  library(dplyr)          # Data manipulation
  library(tidyr)          # Data tidying
  library(scales)         # Axis formatting
  library(patchwork)      # Combining plots
  library(microbiome)     # Additional microbiome utilities
  library(decontam)       # Contamination detection
  library(RColorBrewer)   # Colour palettes
})

# Set global ggplot2 theme
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
# SECTION 1 — DATA IMPORT
# =============================================================================

#' Import microbiome data from multiple formats and build a phyloseq object.
#'
#' @param otu_file    Path to OTU/ASV count table (CSV or TSV). Rows = taxa, Cols = samples.
#' @param tax_file    Path to taxonomy table (CSV or TSV). Rows = taxa, Cols = ranks.
#' @param meta_file   Path to sample metadata (CSV or TSV). Rows = samples.
#' @param tree_file   Optional path to Newick phylogenetic tree.
#' @param sep         Field separator. Default is auto-detected.
#' @param otu_rownames Column index of row names in OTU table. Default = 1.
#' @return A phyloseq object.

import_microbiome_data <- function(otu_file,
                                   tax_file,
                                   meta_file,
                                   tree_file   = NULL,
                                   sep         = NULL,
                                   otu_rownames = 1) {

  cat("=== Importing microbiome data ===\n")

  # --- Auto-detect separator ------------------------------------------------
  detect_sep <- function(file) {
    first_line <- readLines(file, n = 1)
    if (grepl("\t", first_line)) return("\t")
    return(",")
  }

  sep_otu  <- if (is.null(sep)) detect_sep(otu_file)  else sep
  sep_tax  <- if (is.null(sep)) detect_sep(tax_file)  else sep
  sep_meta <- if (is.null(sep)) detect_sep(meta_file) else sep

  # --- Load OTU table -------------------------------------------------------
  cat("  Loading OTU/ASV table:", otu_file, "\n")
  otu_raw <- read.table(otu_file, sep = sep_otu, header = TRUE,
                        row.names = otu_rownames, check.names = FALSE,
                        stringsAsFactors = FALSE)

  # Ensure all values are numeric
  otu_mat <- as.matrix(apply(otu_raw, 2, as.numeric))
  rownames(otu_mat) <- rownames(otu_raw)

  # Remove any columns that are entirely zero
  otu_mat <- otu_mat[, colSums(otu_mat) > 0]

  cat("    OTU table dimensions:", nrow(otu_mat), "taxa x", ncol(otu_mat), "samples\n")

  # --- Load taxonomy table --------------------------------------------------
  cat("  Loading taxonomy table:", tax_file, "\n")
  tax_raw <- read.table(tax_file, sep = sep_tax, header = TRUE,
                        row.names = 1, check.names = FALSE,
                        stringsAsFactors = FALSE)

  # Standardise rank names
  standard_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  colnames(tax_raw) <- standard_ranks[seq_len(ncol(tax_raw))]

  # Remove ambiguous taxonomy strings
  ambiguous_patterns <- c("^k__$", "^p__$", "^c__$", "^o__$", "^f__$",
                          "uncultured", "unknown", "unidentified",
                          "\\[", "\\]", "metagenome")
  for (pattern in ambiguous_patterns) {
    tax_raw[grep(pattern, tax_raw, ignore.case = TRUE)] <- NA
  }

  tax_mat <- as.matrix(tax_raw)
  cat("    Taxonomy table dimensions:", nrow(tax_mat), "taxa x", ncol(tax_mat), "ranks\n")

  # --- Load metadata --------------------------------------------------------
  cat("  Loading metadata:", meta_file, "\n")
  meta_raw <- read.table(meta_file, sep = sep_meta, header = TRUE,
                         row.names = 1, check.names = FALSE,
                         stringsAsFactors = FALSE)

  cat("    Metadata dimensions:", nrow(meta_raw), "samples x", ncol(meta_raw), "variables\n")

  # --- Check sample consistency ---------------------------------------------
  otu_samples  <- colnames(otu_mat)
  meta_samples <- rownames(meta_raw)

  in_both      <- intersect(otu_samples, meta_samples)
  otu_only     <- setdiff(otu_samples, meta_samples)
  meta_only    <- setdiff(meta_samples, otu_samples)

  if (length(otu_only) > 0) {
    warning("Samples in OTU table but not metadata: ",
            paste(otu_only, collapse = ", "))
  }
  if (length(meta_only) > 0) {
    warning("Samples in metadata but not OTU table: ",
            paste(meta_only, collapse = ", "))
  }

  cat("    Matched samples:", length(in_both), "\n")

  # Subset to shared samples
  otu_mat  <- otu_mat[, in_both]
  meta_raw <- meta_raw[in_both, , drop = FALSE]

  # --- Build phyloseq object ------------------------------------------------
  OTU      <- otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX      <- tax_table(tax_mat)
  META     <- sample_data(meta_raw)
  taxa_names(TAX) <- rownames(tax_mat)

  if (!is.null(tree_file) && file.exists(tree_file)) {
    cat("  Loading phylogenetic tree:", tree_file, "\n")
    library(ape)
    tree <- read.tree(tree_file)
    ps   <- phyloseq(OTU, TAX, META, phy_tree(tree))
  } else {
    ps <- phyloseq(OTU, TAX, META)
  }

  cat("  Phyloseq object created successfully.\n")
  cat("    Samples:", nsamples(ps), "\n")
  cat("    Taxa:   ", ntaxa(ps), "\n\n")

  return(ps)
}


# =============================================================================
# SECTION 2 — SEQUENCING DEPTH QC
# =============================================================================

#' Assess sequencing depth across all samples and flag outliers.
#'
#' @param ps          A phyloseq object.
#' @param min_reads   Minimum acceptable reads per sample. Default = 1000.
#' @param group_var   Optional metadata variable for colouring plots.
#' @return A list: plot, summary table, and flagged samples.

qc_sequencing_depth <- function(ps,
                                min_reads = 1000,
                                group_var = NULL) {

  cat("=== Sequencing depth QC ===\n")

  # Calculate read counts per sample
  depth_df <- data.frame(
    sample     = sample_names(ps),
    reads      = sample_sums(ps),
    stringsAsFactors = FALSE
  )

  # Add metadata group variable if provided
  if (!is.null(group_var) && group_var %in% sample_variables(ps)) {
    depth_df[[group_var]] <- as.character(sample_data(ps)[[group_var]])
  }

  depth_df <- depth_df %>%
    arrange(reads) %>%
    mutate(
      sample_order = row_number(),
      flag         = reads < min_reads
    )

  # Summary statistics
  depth_summary <- depth_df %>%
    summarise(
      n_samples     = n(),
      min_reads     = min(reads),
      max_reads     = max(reads),
      median_reads  = median(reads),
      mean_reads    = round(mean(reads)),
      sd_reads      = round(sd(reads)),
      n_below_threshold = sum(flag),
      pct_below     = round(100 * mean(flag), 1)
    )

  cat("  Total samples:", depth_summary$n_samples, "\n")
  cat("  Read count range:", depth_summary$min_reads, "-", depth_summary$max_reads, "\n")
  cat("  Median reads per sample:", depth_summary$median_reads, "\n")
  cat("  Samples below threshold (", min_reads, "):", depth_summary$n_below_threshold, "\n\n")

  # Flagged samples
  flagged <- depth_df %>% filter(flag) %>% pull(sample)

  # --- Plot 1: Bar plot of read depth per sample ----------------------------
  colour_var <- if (!is.null(group_var) && group_var %in% colnames(depth_df)) group_var else NULL

  p1 <- ggplot(depth_df, aes(x = sample_order, y = reads,
                              fill = if (!is.null(colour_var)) .data[[colour_var]] else flag)) +
    geom_col(width = 0.85) +
    geom_hline(yintercept = min_reads, linetype = "dashed",
               colour = "#e74c3c", linewidth = 0.8) +
    annotate("text", x = 1, y = min_reads * 1.05,
             label = paste0("Min threshold: ", format(min_reads, big.mark = ",")),
             hjust = 0, size = 3, colour = "#e74c3c") +
    scale_y_continuous(labels = label_comma()) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title    = "Sequencing depth per sample",
      subtitle = paste0("n = ", nrow(depth_df), " samples | Threshold = ",
                        format(min_reads, big.mark = ","), " reads"),
      x        = "Samples (sorted by read count)",
      y        = "Total reads",
      fill     = if (!is.null(colour_var)) colour_var else "Below threshold",
      caption  = if (length(flagged) > 0) paste("Flagged:", paste(flagged, collapse = ", ")) else "No samples flagged"
    ) +
    theme_microbiome()

  # --- Plot 2: Histogram of read depth distribution -------------------------
  p2 <- ggplot(depth_df, aes(x = reads)) +
    geom_histogram(bins = 30, fill = "#3498db", colour = "white", alpha = 0.85) +
    geom_vline(xintercept = min_reads, linetype = "dashed",
               colour = "#e74c3c", linewidth = 0.8) +
    geom_vline(xintercept = median(depth_df$reads), linetype = "dashed",
               colour = "#27ae60", linewidth = 0.8) +
    annotate("text", x = median(depth_df$reads) * 1.02,
             y = Inf, vjust = 1.5, hjust = 0,
             label = paste0("Median: ", format(round(median(depth_df$reads)), big.mark = ",")),
             size = 3, colour = "#27ae60") +
    scale_x_continuous(labels = label_comma()) +
    labs(
      title    = "Distribution of read depth",
      x        = "Total reads per sample",
      y        = "Number of samples"
    ) +
    theme_microbiome()

  combined_plot <- p1 / p2 +
    plot_annotation(
      title   = "QC Module 1: Sequencing Depth",
      theme   = theme(plot.title = element_text(size = 14, face = "bold"))
    )

  return(list(
    plot     = combined_plot,
    summary  = depth_summary,
    flagged  = flagged,
    data     = depth_df
  ))
}


# =============================================================================
# SECTION 3 — RAREFACTION CURVE ANALYSIS
# =============================================================================

#' Generate rarefaction curves to assess sampling effort sufficiency.
#'
#' @param ps          A phyloseq object.
#' @param step        Step size for rarefaction. Default = 500.
#' @param group_var   Optional metadata variable for colouring curves.
#' @param n_samples   Number of samples to subsample if dataset is large.
#' @return A list: plot and rarefaction data frame.

qc_rarefaction_curves <- function(ps,
                                  step      = 500,
                                  group_var = NULL,
                                  n_samples = NULL) {
  
  cat("=== Rarefaction curve analysis ===\n")
  
  otu_mat <- as(phyloseq::otu_table(ps), "matrix")
  
  if (phyloseq::taxa_are_rows(ps)) {
    otu_t <- t(otu_mat)   # samples x taxa
  } else {
    otu_t <- otu_mat      # already samples x taxa
  }
  
  otu_t <- as.matrix(otu_t)
  storage.mode(otu_t) <- "numeric"
  
  # Subsample if requested
  if (!is.null(n_samples) && n_samples < nrow(otu_t)) {
    set.seed(42)
    idx <- sample(seq_len(nrow(otu_t)), n_samples)
    otu_t <- otu_t[idx, , drop = FALSE]
    cat("  Subsampling to", n_samples, "samples for rarefaction curves.\n")
  }
  
  sample_depths <- rowSums(otu_t)
  max_depth <- min(sample_depths)
  
  rare_steps <- seq(step, max_depth, by = step)
  if (length(rare_steps) == 0 || tail(rare_steps, 1) < max_depth) {
    rare_steps <- c(rare_steps, max_depth)
  }
  
  cat("  Max rarefaction depth:", max_depth, "\n")
  cat("  Steps:", length(rare_steps), "\n")
  
  rare_list <- lapply(rare_steps, function(depth) {
    keep <- sample_depths >= depth
    
    richness <- suppressWarnings(
      vegan::rarefy(otu_t[keep, , drop = FALSE], sample = depth)
    )
    
    data.frame(
      sample   = names(richness),
      depth    = depth,
      richness = as.numeric(richness),
      stringsAsFactors = FALSE
    )
  })
  
  rare_df <- dplyr::bind_rows(rare_list)
  
  if (!is.null(group_var) && group_var %in% phyloseq::sample_variables(ps)) {
    meta_df <- data.frame(phyloseq::sample_data(ps)) %>%
      dplyr::select(dplyr::all_of(group_var)) %>%
      tibble::rownames_to_column("sample")
    
    rare_df <- dplyr::left_join(rare_df, meta_df, by = "sample")
  }
  
  knee_df <- rare_df %>%
    dplyr::group_by(sample) %>%
    dplyr::arrange(depth) %>%
    dplyr::mutate(
      delta_richness = c(NA, diff(richness)),
      slope = delta_richness / step
    ) %>%
    dplyr::filter(!is.na(slope)) %>%
    dplyr::summarise(
      knee_depth = depth[which.min(abs(slope - 0.01))],
      max_richness = max(richness),
      .groups = "drop"
    )
  
  cat("  Median knee point:", median(knee_df$knee_depth), "reads\n\n")
  
  p <- ggplot(rare_df, aes(x = depth, y = richness, group = sample)) +
    {
      if (!is.null(group_var) && group_var %in% colnames(rare_df)) {
        list(
          geom_line(aes(colour = .data[[group_var]]), alpha = 0.6, linewidth = 0.5),
          scale_colour_brewer(palette = "Set2")
        )
      } else {
        list(
          geom_line(colour = "#3498db", alpha = 0.5, linewidth = 0.5)
        )
      }
    } +
    geom_vline(
      xintercept = median(knee_df$knee_depth),
      linetype = "dotted",
      colour = "#e67e22",
      linewidth = 0.9
    ) +
    annotate(
      "text",
      x = median(knee_df$knee_depth) * 1.02,
      y = max(rare_df$richness) * 0.95,
      label = paste0(
        "Median plateau:\n",
        format(round(median(knee_df$knee_depth)), big.mark = ",")
      ),
      hjust = 0,
      size = 3,
      colour = "#e67e22"
    ) +
    scale_x_continuous(labels = label_comma()) +
    labs(
      title = "Rarefaction curves",
      subtitle = "Sampling effort assessment - curves should plateau before the minimum depth",
      x = "Sequencing depth (reads)",
      y = "Observed ASV richness",
      colour = group_var,
      caption = paste0(
        "Min sample depth: ",
        format(min(sample_depths), big.mark = ",")
      )
    ) +
    theme_microbiome()
  
  return(list(
    plot = p,
    data = rare_df,
    knee = knee_df
  ))
}


# =============================================================================
# SECTION 4 — TAXA PREVALENCE AND ABUNDANCE FILTERING
# =============================================================================

#' Filter low-prevalence and low-abundance taxa.
#'
#' @param ps              A phyloseq object.
#' @param min_prevalence  Minimum fraction of samples a taxon must appear in. Default = 0.05.
#' @param min_abundance   Minimum total read count across all samples. Default = 10.
#' @param plot_filter     Whether to produce a prevalence-abundance scatter plot. Default = TRUE.
#' @return A list: filtered phyloseq, plot, and filtering summary.

qc_filter_taxa <- function(ps,
                           min_prevalence = 0.05,
                           min_abundance  = 10,
                           plot_filter    = TRUE) {

  cat("=== Taxa prevalence and abundance filtering ===\n")

  otu_mat  <- as.matrix(otu_table(ps))
  if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)

  n_samples <- ncol(otu_mat)

  # Compute per-taxon metrics
  taxa_df <- data.frame(
    taxon       = rownames(otu_mat),
    prevalence  = rowSums(otu_mat > 0) / n_samples,
    total_abund = rowSums(otu_mat),
    mean_abund  = rowMeans(otu_mat),
    max_abund   = apply(otu_mat, 1, max),
    stringsAsFactors = FALSE
  )

  # Apply filters
  pass <- taxa_df$prevalence >= min_prevalence &
          taxa_df$total_abund >= min_abundance

  taxa_df$pass <- pass

  n_before <- ntaxa(ps)
  n_after  <- sum(pass)
  n_removed <- n_before - n_after
  pct_removed <- round(100 * n_removed / n_before, 1)

  cat("  Taxa before filtering:", n_before, "\n")
  cat("  Taxa passing filters: ", n_after, "\n")
  cat("  Taxa removed:         ", n_removed, "(", pct_removed, "%)\n\n")

  # Filter the phyloseq object
  keep_taxa <- taxa_df$taxon[pass]
  ps_filtered <- prune_taxa(keep_taxa, ps)

  # --- Prevalence-abundance scatter plot -----------------------------------
  p <- NULL
  if (plot_filter) {
    # Add phylum information for colouring
    if (!is.null(tax_table(ps))) {
      tax_df <- as.data.frame(tax_table(ps)) %>%
        tibble::rownames_to_column("taxon") %>%
        dplyr::select(taxon, Phylum) %>%
        mutate(Phylum = ifelse(is.na(Phylum) | Phylum == "", "Unknown", Phylum))
      taxa_df <- left_join(taxa_df, tax_df, by = "taxon")

      # Keep only top phyla for colouring
      top_phyla <- taxa_df %>%
        group_by(Phylum) %>%
        summarise(total = sum(total_abund)) %>%
        arrange(desc(total)) %>%
        slice_head(n = 8) %>%
        pull(Phylum)

      taxa_df <- taxa_df %>%
        mutate(Phylum_plot = ifelse(Phylum %in% top_phyla, Phylum, "Other"))
    }

    p <- ggplot(taxa_df, aes(x = prevalence, y = log10(total_abund + 1))) +
      {
        if ("Phylum_plot" %in% colnames(taxa_df)) {
          list(
            geom_point(aes(colour = Phylum_plot, shape = pass),
                       alpha = 0.6, size = 1.5),
            scale_colour_brewer(palette = "Paired", name = "Phylum")
          )
        } else {
          list(
            geom_point(aes(shape = pass, colour = pass),
                       alpha = 0.6, size = 1.5)
          )
        }
      } +
      geom_vline(xintercept = min_prevalence,
                 linetype = "dashed", colour = "#e74c3c", linewidth = 0.8) +
      geom_hline(yintercept = log10(min_abundance + 1),
                 linetype = "dashed", colour = "#e74c3c", linewidth = 0.8) +
      scale_shape_manual(values = c(`FALSE` = 4, `TRUE` = 16),
                         labels = c("Removed", "Kept"),
                         name   = "Filter") +
      annotate("text", x = min_prevalence + 0.01, y = 0.2,
               label = paste0("Min prevalence: ", min_prevalence * 100, "%"),
               hjust = 0, size = 3, colour = "#e74c3c") +
      annotate("text", x = 0.01, y = log10(min_abundance + 1) + 0.1,
               label = paste0("Min abundance: ", min_abundance),
               hjust = 0, size = 3, colour = "#e74c3c") +
      labs(
        title    = "Taxa prevalence vs. abundance",
        subtitle = paste0("Removed ", n_removed, " of ", n_before,
                          " taxa (", pct_removed, "%)"),
        x        = "Prevalence (fraction of samples)",
        y        = "Log10 total abundance"
      ) +
      theme_microbiome()
  }

  return(list(
    ps_filtered = ps_filtered,
    plot        = p,
    summary     = list(
      n_before    = n_before,
      n_after     = n_after,
      n_removed   = n_removed,
      pct_removed = pct_removed
    ),
    taxa_stats  = taxa_df
  ))
}


# =============================================================================
# SECTION 5 — CONTAMINATION DETECTION
# =============================================================================

#' Detect potential contaminants using decontam.
#'
#' @param ps              A phyloseq object.
#' @param method          "frequency" (needs DNA concentration), "prevalence"
#'                        (needs negative control column), or "combined".
#' @param neg_control_var Metadata column name for negative control indicator.
#' @param conc_var        Metadata column name for DNA concentration.
#' @param threshold       Probability threshold for contamination. Default = 0.1.
#' @return A list: cleaned phyloseq, contaminant taxa, and plot.

qc_decontam <- function(ps,
                        method          = "prevalence",
                        neg_control_var = "is_negative_control",
                        conc_var        = "dna_concentration",
                        threshold       = 0.1) {

  cat("=== Contamination detection (decontam) ===\n")

  if (!requireNamespace("decontam", quietly = TRUE)) {
    cat("  decontam not available. Skipping contamination detection.\n")
    return(list(ps_clean = ps, contaminants = character(0), plot = NULL))
  }

  # Check required metadata columns
  available_vars <- sample_variables(ps)

  if (method == "prevalence" || method == "combined") {
    if (!neg_control_var %in% available_vars) {
      cat("  Column '", neg_control_var, "' not found in metadata.\n")
      cat("  Skipping contamination detection.\n")
      return(list(ps_clean = ps, contaminants = character(0), plot = NULL))
    }
    negative_controls <- sample_data(ps)[[neg_control_var]]
  }

  if (method == "frequency" || method == "combined") {
    if (!conc_var %in% available_vars) {
      cat("  Column '", conc_var, "' not found in metadata.\n")
      cat("  Skipping contamination detection.\n")
      return(list(ps_clean = ps, contaminants = character(0), plot = NULL))
    }
    dna_conc <- as.numeric(sample_data(ps)[[conc_var]])
  }

  # Run decontam
  contam_df <- tryCatch({
    if (method == "prevalence") {
      isContaminant(ps, method = "prevalence",
                    neg      = negative_controls,
                    threshold = threshold)
    } else if (method == "frequency") {
      isContaminant(ps, method = "frequency",
                    conc     = dna_conc,
                    threshold = threshold)
    } else {
      isContaminant(ps, method = "combined",
                    neg      = negative_controls,
                    conc     = dna_conc,
                    threshold = threshold)
    }
  }, error = function(e) {
    cat("  decontam error:", conditionMessage(e), "\n")
    return(NULL)
  })

  if (is.null(contam_df)) {
    return(list(ps_clean = ps, contaminants = character(0), plot = NULL))
  }

  # Identify contaminants
  contaminants <- rownames(contam_df)[contam_df$contaminant == TRUE]
  n_contam     <- length(contaminants)

  cat("  Contaminants identified:", n_contam, "\n")
  if (n_contam > 0) {
    cat("  Contaminant taxa:\n")
    cat(paste0("    - ", contaminants, "\n"))
  }

  # Remove contaminants
  ps_clean <- prune_taxa(!contam_df$contaminant, ps)
  cat("  Taxa after decontam:", ntaxa(ps_clean), "\n\n")

  # --- Plot: Decontam score distribution ------------------------------------
  contam_df$taxon      <- rownames(contam_df)
  contam_df$status     <- ifelse(contam_df$contaminant, "Contaminant", "Retained")
  contam_df$prev_score <- contam_df$p

  p <- ggplot(contam_df, aes(x = prev_score, fill = status)) +
    geom_histogram(bins = 40, colour = "white", alpha = 0.85) +
    geom_vline(xintercept = threshold, linetype = "dashed",
               colour = "#e74c3c", linewidth = 0.8) +
    scale_fill_manual(values = c("Contaminant" = "#e74c3c", "Retained" = "#27ae60")) +
    labs(
      title    = "Contamination probability scores (decontam)",
      subtitle = paste0("Method: ", method, " | Threshold: ", threshold,
                        " | Contaminants removed: ", n_contam),
      x        = "Contamination probability",
      y        = "Number of taxa",
      fill     = "Status"
    ) +
    theme_microbiome()

  return(list(
    ps_clean     = ps_clean,
    contaminants = contaminants,
    scores       = contam_df,
    plot         = p
  ))
}


# =============================================================================
# SECTION 6 — SAMPLE FILTERING
# =============================================================================

#' Remove samples that fail QC thresholds.
#'
#' @param ps          A phyloseq object.
#' @param min_reads   Minimum reads to retain a sample. Default = 1000.
#' @param min_taxa    Minimum number of taxa to retain a sample. Default = 10.
#' @return A list: filtered phyloseq and removed samples.

qc_filter_samples <- function(ps,
                               min_reads = 1000,
                               min_taxa  = 10) {

  cat("=== Sample filtering ===\n")

  reads_per_sample <- sample_sums(ps)
  taxa_per_sample  <- colSums(as.matrix(otu_table(ps)) > 0)
  if (!taxa_are_rows(ps)) taxa_per_sample <- rowSums(as.matrix(otu_table(ps)) > 0)

  pass_reads <- reads_per_sample >= min_reads
  pass_taxa  <- taxa_per_sample  >= min_taxa
  pass_both  <- pass_reads & pass_taxa

  removed_reads <- names(pass_reads)[!pass_reads]
  removed_taxa  <- names(pass_taxa)[!pass_taxa]
  removed_all   <- names(pass_both)[!pass_both]

  cat("  Samples before filtering:", nsamples(ps), "\n")
  cat("  Removed (low reads <", min_reads, "):", length(removed_reads), "\n")
  cat("  Removed (low taxa  <", min_taxa,  "):", length(removed_taxa), "\n")
  cat("  Total removed:           ", length(removed_all), "\n")

  ps_filtered <- prune_samples(pass_both, ps)
  cat("  Samples after filtering: ", nsamples(ps_filtered), "\n\n")

  return(list(
    ps_filtered     = ps_filtered,
    removed_samples = removed_all,
    removed_reads   = removed_reads,
    removed_taxa    = removed_taxa
  ))
}


# =============================================================================
# SECTION 7 — QC SUMMARY REPORT
# =============================================================================

#' Generate a comprehensive QC summary table and final diagnostic plot.
#'
#' @param ps_raw      The original unfiltered phyloseq object.
#' @param ps_filtered The filtered phyloseq object.
#' @param depth_result Output from qc_sequencing_depth().
#' @param filter_result Output from qc_filter_taxa().
#' @return A list: summary data frame and panel plot.

qc_summary_report <- function(ps_raw,
                              ps_filtered,
                              depth_result  = NULL,
                              filter_result = NULL) {

  cat("=== QC Summary Report ===\n")

  summary_df <- data.frame(
    Step                  = c("Raw data",
                               "After sample filtering",
                               "After taxa filtering",
                               "Final dataset"),
    N_samples             = c(nsamples(ps_raw),
                               NA,
                               NA,
                               nsamples(ps_filtered)),
    N_taxa                = c(ntaxa(ps_raw),
                               NA,
                               NA,
                               ntaxa(ps_filtered)),
    Total_reads           = c(sum(sample_sums(ps_raw)),
                               NA,
                               NA,
                               sum(sample_sums(ps_filtered))),
    Median_reads_sample   = c(median(sample_sums(ps_raw)),
                               NA,
                               NA,
                               median(sample_sums(ps_filtered)))
  )

  print(summary_df)

  # --- Final sample summary plot: reads vs taxa coloured by QC pass --------
  otu_mat <- as(phyloseq::otu_table(ps_filtered), "matrix")
  
  if (phyloseq::taxa_are_rows(ps_filtered)) {
    n_taxa <- colSums(otu_mat > 0)
  } else {
    n_taxa <- rowSums(otu_mat > 0)
  }
  
  final_df <- data.frame(
    sample = phyloseq::sample_names(ps_filtered),
    reads  = phyloseq::sample_sums(ps_filtered),
    n_taxa = n_taxa,
    stringsAsFactors = FALSE
  )
  

  p_final <- ggplot(final_df, aes(x = reads, y = n_taxa)) +
    geom_point(colour = "#27ae60", alpha = 0.75, size = 3) +
    geom_smooth(method = "lm", se = TRUE, colour = "#2980b9", linetype = "dashed",
                linewidth = 0.8, alpha = 0.15) +
    scale_x_continuous(labels = label_comma()) +
    labs(
      title    = "Final dataset: reads vs taxon richness per sample",
      subtitle = paste0(nsamples(ps_filtered), " samples | ",
                        ntaxa(ps_filtered), " taxa | ",
                        format(sum(sample_sums(ps_filtered)), big.mark = ","),
                        " total reads"),
      x        = "Total reads",
      y        = "Observed taxa"
    ) +
    theme_microbiome()

  return(list(
    summary = summary_df,
    plot    = p_final
  ))
}


# =============================================================================
# SECTION 8 — COMPLETE QC WORKFLOW WRAPPER
# =============================================================================

#' Run the complete QC pipeline from import to final filtered phyloseq.
#'
#' @param otu_file       Path to OTU/ASV count table.
#' @param tax_file       Path to taxonomy table.
#' @param meta_file      Path to metadata.
#' @param tree_file      Optional path to phylogenetic tree.
#' @param min_reads      Minimum reads per sample. Default = 1000.
#' @param min_prevalence Minimum taxon prevalence. Default = 0.05.
#' @param min_abundance  Minimum taxon total abundance. Default = 10.
#' @param group_var      Metadata variable for grouping in plots.
#' @param output_dir     Directory to save plots and reports.
#' @return Final filtered phyloseq object, with all QC plots saved.

run_microbiome_qc <- function(otu_file,
                               tax_file,
                               meta_file,
                               tree_file      = NULL,
                               min_reads      = 1000,
                               min_prevalence = 0.05,
                               min_abundance  = 10,
                               group_var      = NULL,
                               output_dir     = "qc_output") {

  cat("\n", strrep("=", 60), "\n")
  cat("  MICROBIOME QC PIPELINE\n")
  cat(strrep("=", 60), "\n\n")

  # Create output directory
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # --- Step 1: Import -------------------------------------------------------
  ps_raw <- import_microbiome_data(
    otu_file  = otu_file,
    tax_file  = tax_file,
    meta_file = meta_file,
    tree_file = tree_file
  )

  # --- Step 2: Sequencing depth QC -----------------------------------------
  depth_res <- qc_sequencing_depth(ps_raw, min_reads = min_reads,
                                   group_var = group_var)
  ggsave(file.path(output_dir, "01_sequencing_depth.pdf"),
         depth_res$plot, width = 12, height = 10)

  # --- Step 3: Sample filtering ---------------------------------------------
  sample_res <- qc_filter_samples(ps_raw, min_reads = min_reads)
  ps_sample_filtered <- sample_res$ps_filtered

  # --- Step 4: Rarefaction curves -------------------------------------------
  rare_res <- qc_rarefaction_curves(ps_sample_filtered,
                                    group_var = group_var)
  ggsave(file.path(output_dir, "02_rarefaction_curves.pdf"),
         rare_res$plot, width = 10, height = 7)

  # --- Step 5: Taxa filtering -----------------------------------------------
  taxa_res <- qc_filter_taxa(ps_sample_filtered,
                              min_prevalence = min_prevalence,
                              min_abundance  = min_abundance)
  ggsave(file.path(output_dir, "03_taxa_prevalence_abundance.pdf"),
         taxa_res$plot, width = 10, height = 7)
  ps_taxa_filtered <- taxa_res$ps_filtered

  # --- Step 6: Contamination detection (if applicable) ---------------------
  decontam_res <- qc_decontam(ps_taxa_filtered)
  ps_clean <- decontam_res$ps_clean
  if (!is.null(decontam_res$plot)) {
    ggsave(file.path(output_dir, "04_decontam.pdf"),
           decontam_res$plot, width = 9, height = 6)
  }

  # --- Step 7: Summary report -----------------------------------------------
  summary_res <- qc_summary_report(
    ps_raw      = ps_raw,
    ps_filtered = ps_clean,
    depth_result  = depth_res,
    filter_result = taxa_res
  )
  ggsave(file.path(output_dir, "05_final_qc_summary.pdf"),
         summary_res$plot, width = 8, height = 6)

  # Save filtered phyloseq object
  saveRDS(ps_clean, file.path(output_dir, "phyloseq_qc_filtered.rds"))

  cat("\n", strrep("=", 60), "\n")
  cat("  QC PIPELINE COMPLETE\n")
  cat(strrep("=", 60), "\n")
  cat("  Output saved to:", output_dir, "\n")
  cat("  Final phyloseq : phyloseq_qc_filtered.rds\n")
  cat("  Plots saved    : 5 PDF files\n")
  cat("  Final dataset  :", nsamples(ps_clean), "samples,",
      ntaxa(ps_clean), "taxa\n\n")

  return(ps_clean)
}


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# --- Option A: Run the complete pipeline -------------------------------------
# ps_filtered <- run_microbiome_qc(
#   otu_file       = "data/otu_table.csv",
#   tax_file       = "data/taxonomy.csv",
#   meta_file      = "data/metadata.csv",
#   tree_file      = "data/tree.nwk",       # optional
#   min_reads      = 5000,
#   min_prevalence = 0.05,
#   min_abundance  = 10,
#   group_var      = "disease_status",
#   output_dir     = "results/qc"
# )

# --- Option B: Run steps individually ----------------------------------------
# ps          <- import_microbiome_data("otu.csv", "tax.csv", "meta.csv")
# depth_res   <- qc_sequencing_depth(ps, min_reads = 5000, group_var = "group")
# depth_res$plot
# depth_res$flagged
# rare_res    <- qc_rarefaction_curves(ps, group_var = "group")
# rare_res$plot
# filter_res  <- qc_filter_taxa(ps, min_prevalence = 0.05, min_abundance = 10)
# filter_res$plot
# ps_clean    <- filter_res$ps_filtered
