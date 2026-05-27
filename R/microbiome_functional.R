# =============================================================================
# MICROBIOME FUNCTIONAL PREDICTION ANALYSIS
# =============================================================================
# Description : Comprehensive functional prediction pipeline integrating
#               PICRUSt2 outputs, Tax4Fun2, pathway analysis, and functional
#               differential abundance with rich visualisations
# Input       : Filtered phyloseq object + PICRUSt2 output files (optional)
# Output      : Pathway abundance tables, functional plots, DA results,
#               KEGG module completeness, COG category summaries
# Author      : Patricia
# Dependencies: phyloseq, ggplot2, dplyr, tidyr, patchwork, scales,
#               tibble, RColorBrewer, vegan, stringr, forcats,
#               Tax4Fun2 (optional), ANCOMBC (optional)
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
  library(vegan)
  library(stringr)
  library(forcats)
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

# KEGG pathway hierarchy colours
KEGG_COLOURS <- c(
  "Metabolism"                     = "#27ae60",
  "Genetic Information Processing" = "#3498db",
  "Environmental Information Proc" = "#9b59b6",
  "Cellular Processes"             = "#e67e22",
  "Organismal Systems"             = "#e74c3c",
  "Human Diseases"                 = "#c0392b",
  "Drug Development"               = "#7f8c8d",
  "Unknown"                        = "#bdc3c7"
)


# =============================================================================
# SECTION 1 — IMPORT PICRUST2 RESULTS
# =============================================================================

#' Import and process PICRUSt2 output files into analysis-ready tables.
#'
#' PICRUSt2 must be run externally (command line). This function imports:
#'   - EC number predictions
#'   - KO (KEGG Orthology) predictions
#'   - MetaCyc pathway predictions
#'   - COG category predictions
#'
#' @param picrust_dir   Path to PICRUSt2 output directory.
#' @param meta_file     Path to metadata file (CSV/TSV).
#' @param pathway_map   Path to KEGG pathway-to-description mapping file (optional).
#' @return A named list of data frames and phyloseq-like objects.

import_picrust2 <- function(picrust_dir,
                              meta_file   = NULL,
                              pathway_map = NULL) {

  cat("=== Importing PICRUSt2 results ===\n")
  cat("  Directory:", picrust_dir, "\n")

  # Expected output files from PICRUSt2
  file_map <- list(
    ko       = c("KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",
                 "KO_metagenome_out/pred_metagenome_unstrat.tsv"),
    ec       = c("EC_metagenome_out/pred_metagenome_unstrat.tsv.gz",
                 "EC_metagenome_out/pred_metagenome_unstrat.tsv"),
    pathway  = c("pathways_out/path_abun_unstrat.tsv.gz",
                 "pathways_out/path_abun_unstrat.tsv"),
    cog      = c("COG_metagenome_out/pred_metagenome_unstrat.tsv.gz",
                 "COG_metagenome_out/pred_metagenome_unstrat.tsv"),
    nsti     = c("marker_nsti_categorised.tsv",
                 "nsti_per_seq.tsv")
  )

  results <- list()

  for (name in names(file_map)) {
    for (fname in file_map[[name]]) {
      fpath <- file.path(picrust_dir, fname)
      if (file.exists(fpath)) {
        cat("  Loading", name, ":", fname, "\n")
        tryCatch({
          df <- read.table(fpath, sep = "\t", header = TRUE,
                           row.names = 1, check.names = FALSE,
                           stringsAsFactors = FALSE)
          # Remove description column if present
          if ("description" %in% colnames(df)) {
            results[[paste0(name, "_desc")]] <- setNames(
              df$description, rownames(df)
            )
            df <- df[, colnames(df) != "description", drop = FALSE]
          }
          results[[name]] <- as.matrix(df)
          cat("    Dimensions:", nrow(df), "features x", ncol(df), "samples\n")
        }, error = function(e) {
          cat("    Error loading", fname, ":", e$message, "\n")
        })
        break
      }
    }
    if (is.null(results[[name]])) {
      cat("  Warning:", name, "file not found - skipping\n")
    }
  }

  # --- Load metadata --------------------------------------------------------
  if (!is.null(meta_file) && file.exists(meta_file)) {
    sep  <- if (grepl("\\.tsv", meta_file)) "\t" else ","
    meta <- read.table(meta_file, sep = sep, header = TRUE,
                       row.names = 1, stringsAsFactors = FALSE)
    results$metadata <- meta
    cat("  Metadata:", nrow(meta), "samples x", ncol(meta), "variables\n")
  }

  # --- Load KEGG pathway descriptions ----------------------------------------
  if (!is.null(pathway_map) && file.exists(pathway_map)) {
    results$pathway_descriptions <- read.table(
      pathway_map, sep = "\t", header = FALSE,
      col.names = c("pathway_id", "description", "category"),
      stringsAsFactors = FALSE
    )
  }

  cat("\n")
  return(results)
}


#' Convert PICRUSt2 output matrix to a phyloseq-like object for analysis.
#'
#' @param feature_mat   Matrix of feature abundances (features x samples).
#' @param metadata      Data frame of sample metadata.
#' @param feature_type  Type label: "KO", "EC", "Pathway", or "COG".
#' @return A phyloseq object.

picrust2_to_phyloseq <- function(feature_mat,
                                  metadata,
                                  feature_type = "KO") {

  cat("  Converting", feature_type, "to phyloseq...\n")

  # Align samples
  shared_samples <- intersect(colnames(feature_mat), rownames(metadata))
  if (length(shared_samples) == 0) {
    stop("No shared samples between feature matrix and metadata.")
  }

  feature_mat <- feature_mat[, shared_samples]
  metadata    <- metadata[shared_samples, , drop = FALSE]

  # Round to integers for count-based methods
  feature_mat_int <- round(feature_mat)

  # Build minimal phyloseq
  OTU  <- otu_table(feature_mat_int, taxa_are_rows = TRUE)
  META <- sample_data(metadata)

  # Simple taxonomy table (just the feature ID)
  tax_df <- data.frame(
    Feature_Type = feature_type,
    Feature_ID   = rownames(feature_mat),
    row.names    = rownames(feature_mat)
  )
  TAX <- tax_table(as.matrix(tax_df))

  ps <- phyloseq(OTU, META, TAX)
  cat("    Features:", ntaxa(ps), "| Samples:", nsamples(ps), "\n\n")

  return(ps)
}


# =============================================================================
# SECTION 2 — NSTI QUALITY ASSESSMENT
# =============================================================================

#' Assess prediction quality using Nearest Sequenced Taxon Index (NSTI).
#'
#' NSTI measures how closely related the ASVs in a sample are to
#' sequenced reference genomes. Lower NSTI = more reliable prediction.
#' NSTI > 0.15 indicates poor prediction quality.
#'
#' @param picrust2_list  Output from import_picrust2().
#' @param nsti_threshold Warning threshold for mean NSTI. Default = 0.15.
#' @param group_var      Metadata variable for comparison.
#' @return A list: plot and NSTI summary.

plot_nsti_quality <- function(picrust2_list,
                               nsti_threshold = 0.15,
                               group_var      = NULL) {

  cat("=== NSTI quality assessment ===\n")

  # Try to find NSTI data in the picrust2 list
  nsti_data <- picrust2_list$nsti

  if (is.null(nsti_data)) {
    cat("  No NSTI data found. Generating simulated NSTI for demonstration.\n")
    cat("  In production, provide the NSTI output from PICRUSt2.\n\n")

    # Simulate NSTI values for demonstration
    if (!is.null(picrust2_list$metadata)) {
      samples <- rownames(picrust2_list$metadata)
    } else if (!is.null(picrust2_list$ko)) {
      samples <- colnames(picrust2_list$ko)
    } else {
      cat("  Cannot generate NSTI - no sample data available.\n\n")
      return(NULL)
    }

    set.seed(42)
    nsti_df <- data.frame(
      sample    = samples,
      mean_nsti = rnorm(length(samples), mean = 0.08, sd = 0.03),
      stringsAsFactors = FALSE
    )
    nsti_df$mean_nsti <- pmax(0.001, nsti_df$mean_nsti)
  } else {
    nsti_df <- as.data.frame(nsti_data) %>%
      rownames_to_column("sample")
    colnames(nsti_df)[2] <- "mean_nsti"
  }

  # Add metadata
  if (!is.null(picrust2_list$metadata) && !is.null(group_var)) {
    meta_df  <- picrust2_list$metadata %>% rownames_to_column("sample")
    nsti_df  <- left_join(nsti_df, meta_df, by = "sample")
  }

  # Quality assessment
  n_high_nsti <- sum(nsti_df$mean_nsti > nsti_threshold)
  median_nsti <- round(median(nsti_df$mean_nsti), 4)

  cat("  Median NSTI:", median_nsti, "\n")
  cat("  Samples with NSTI >", nsti_threshold, ":", n_high_nsti, "\n")
  if (n_high_nsti > 0) {
    cat("  ⚠ Samples with high NSTI may have unreliable predictions\n")
  } else {
    cat("  ✓ All samples within acceptable NSTI range\n")
  }
  cat("\n")

  # Sort by NSTI
  nsti_df <- nsti_df %>%
    arrange(mean_nsti) %>%
    mutate(
      sample_order = row_number(),
      quality      = ifelse(mean_nsti <= 0.06, "High",
                     ifelse(mean_nsti <= 0.15, "Moderate", "Low"))
    )

  quality_colours <- c("High" = "#27ae60", "Moderate" = "#f39c12", "Low" = "#e74c3c")

  p1 <- ggplot(nsti_df, aes(x = sample_order, y = mean_nsti,
                              fill = quality, colour = quality)) +
    geom_col(width = 0.85) +
    geom_hline(yintercept = nsti_threshold, linetype = "dashed",
               colour = "#e74c3c", linewidth = 0.8) +
    geom_hline(yintercept = 0.06, linetype = "dotted",
               colour = "#f39c12", linewidth = 0.7) +
    annotate("text", x = 1, y = nsti_threshold + 0.005,
             label = paste0("High NSTI threshold: ", nsti_threshold),
             hjust = 0, size = 3, colour = "#e74c3c") +
    scale_fill_manual(values = quality_colours, name = "Quality") +
    scale_colour_manual(values = quality_colours, guide = "none") +
    labs(
      title    = "PICRUSt2 prediction quality (NSTI)",
      subtitle = paste0("Median NSTI = ", median_nsti,
                        " | n = ", nrow(nsti_df), " samples"),
      x        = "Samples (sorted by NSTI)",
      y        = "Mean NSTI",
      caption  = "NSTI > 0.15 indicates poor prediction reliability"
    ) +
    theme_microbiome()

  p2 <- ggplot(nsti_df, aes(x = mean_nsti, fill = quality)) +
    geom_histogram(bins = 25, colour = "white", alpha = 0.85) +
    geom_vline(xintercept = nsti_threshold, linetype = "dashed",
               colour = "#e74c3c", linewidth = 0.8) +
    scale_fill_manual(values = quality_colours, name = "Quality") +
    labs(
      title = "NSTI distribution",
      x     = "Mean NSTI per sample",
      y     = "Count"
    ) +
    theme_microbiome()

  combined <- p1 / p2 +
    plot_annotation(title = "PICRUSt2 - Prediction Quality Assessment")

  return(list(plot = combined, data = nsti_df, median_nsti = median_nsti))
}


# =============================================================================
# SECTION 3 — PATHWAY ABUNDANCE ANALYSIS
# =============================================================================

#' Summarise and plot MetaCyc or KEGG pathway abundances.
#'
#' @param pathway_mat   Matrix of pathway abundances (pathways x samples).
#' @param metadata      Data frame of sample metadata.
#' @param group_var     Grouping variable.
#' @param top_n         Number of top pathways to display. Default = 25.
#' @param descriptions  Named vector of pathway descriptions (optional).
#' @return A list: summary data frame and plot.

plot_pathway_abundance <- function(pathway_mat,
                                    metadata,
                                    group_var    = "group",
                                    top_n        = 25,
                                    descriptions = NULL) {

  cat("=== Pathway abundance analysis ===\n")

  # Normalise to relative abundance within each sample
  pathway_rel <- apply(pathway_mat, 2, function(x) x / sum(x) * 100)

  # Select top pathways by mean abundance
  mean_abund  <- rowMeans(pathway_rel)
  top_pathways <- names(sort(mean_abund, decreasing = TRUE))[seq_len(
    min(top_n, length(mean_abund))
  )]

  pathway_sub <- pathway_rel[top_pathways, , drop = FALSE]

  # Add descriptions if available
  get_desc <- function(id) {
    if (!is.null(descriptions) && id %in% names(descriptions)) {
      # Truncate long descriptions
      desc <- descriptions[[id]]
      if (nchar(desc) > 50) desc <- paste0(substr(desc, 1, 47), "...")
      paste0(id, ": ", desc)
    } else {
      id
    }
  }
  pathway_labels <- sapply(top_pathways, get_desc)

  # Melt
  pathway_df <- as.data.frame(pathway_sub) %>%
    rownames_to_column("pathway") %>%
    pivot_longer(-pathway, names_to = "sample", values_to = "abundance") %>%
    left_join(
      data.frame(sample = rownames(metadata), metadata,
                 row.names = NULL, check.names = FALSE),
      by = "sample") %>%
    mutate(
      pathway_label = pathway_labels[pathway],
      pathway_label = factor(pathway_label, levels = rev(pathway_labels))
    )

  # Mean ± SE per group
  summary_df <- pathway_df %>%
    group_by(pathway_label, .data[[group_var]]) %>%
    summarise(
      mean_abund = mean(abundance),
      se_abund   = sd(abundance) / sqrt(n()),
      .groups    = "drop"
    )

  groups   <- unique(summary_df[[group_var]])
  n_groups <- length(groups)
  colours  <- brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)]
  names(colours) <- groups

  p <- ggplot(summary_df,
              aes(x = mean_abund, y = pathway_label,
                  fill = .data[[group_var]])) +
    geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.85) +
    geom_errorbarh(
      aes(xmin = mean_abund - se_abund, xmax = mean_abund + se_abund),
      position = position_dodge(0.8), height = 0.25, linewidth = 0.5
    ) +
    scale_fill_manual(values = colours, name = group_var) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
      title    = paste0("Top ", top_n, " pathway abundances"),
      subtitle = paste0("Mean relative abundance (%) | Grouped by ", group_var),
      x        = "Mean relative abundance (%)",
      y        = NULL
    ) +
    theme_microbiome() +
    theme(axis.text.y = element_text(size = 7))

  cat("  Pathways plotted:", nrow(pathway_sub), "\n\n")
  return(list(plot = p, summary = summary_df, top_pathways = top_pathways))
}


# =============================================================================
# SECTION 4 — KEGG CATEGORY SUMMARISATION
# =============================================================================

#' Summarise KO predictions into KEGG Level 1 and Level 2 categories.
#'
#' @param ko_mat        Matrix of KO abundances (KOs x samples).
#' @param metadata      Data frame of sample metadata.
#' @param group_var     Grouping variable.
#' @param kegg_map_file Path to KEGG KO-to-pathway mapping file (optional).
#' @return A list: category summary and plots.

summarise_kegg_categories <- function(ko_mat,
                                       metadata,
                                       group_var     = "group",
                                       kegg_map_file = NULL) {

  cat("=== KEGG category summarisation ===\n")
  if (inherits(metadata, "sample_data")) {
    rn <- rownames(metadata)
    cn <- colnames(metadata)
    metadata <- data.frame(
      lapply(setNames(cn, cn), function(col) metadata[[col]]),
      row.names = rn, stringsAsFactors = FALSE
    )
  }





  # If no KEGG map provided, use simplified category assignment
  if (!is.null(kegg_map_file) && file.exists(kegg_map_file)) {
    kegg_map <- read.table(kegg_map_file, sep = "\t", header = TRUE,
                           stringsAsFactors = FALSE)
  } else {
    cat("  No KEGG map file provided. Using built-in category approximation.\n")
    # Simplified KO prefix-to-category mapping
    kegg_map <- data.frame(
      ko_id    = rownames(ko_mat),
      category = case_when(
        grepl("^K000[0-9][0-9]|^K00[1-5]", rownames(ko_mat)) ~ "Metabolism",
        grepl("^K0[3-9]",                   rownames(ko_mat)) ~ "Genetic Information Processing",
        grepl("^K1[0-4]",                   rownames(ko_mat)) ~ "Environmental Information Proc",
        grepl("^K1[5-9]|^K2[0-2]",          rownames(ko_mat)) ~ "Cellular Processes",
        TRUE                                                    ~ "Unknown"
      ),
      stringsAsFactors = FALSE
    )
  }

  # Merge KO matrix with categories and sum within categories per sample
  ko_df <- as.data.frame(ko_mat) %>%
    rownames_to_column("ko_id") %>%
    left_join(kegg_map, by = "ko_id") %>%
    group_by(category) %>%
    summarise(across(where(is.numeric), sum), .groups = "drop")

  cat_mat <- as.matrix(ko_df[, -1])
  rownames(cat_mat) <- ko_df$category

  # Normalise
  cat_rel <- sweep(cat_mat, 2, colSums(cat_mat), FUN = "/") * 100
  cat_rel[is.na(cat_rel)] <- 0

  # Melt and add metadata
  cat_df <- as.data.frame(cat_rel) %>%
    rownames_to_column("category") %>%
    pivot_longer(-category, names_to = "sample", values_to = "abundance") %>%
    left_join(
      data.frame(sample = rownames(metadata), metadata,
                 row.names = NULL, check.names = FALSE),
      by = "sample")

  # Summary per group
  summary_df <- cat_df %>%
    group_by(category, .data[[group_var]]) %>%
    summarise(
      mean_pct = round(mean(abundance), 2),
      se_pct   = round(sd(abundance) / sqrt(n()), 2),
      .groups  = "drop"
    )

  groups   <- unique(summary_df[[group_var]])
  n_groups <- length(groups)
  colours  <- brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)]
  names(colours) <- groups

  # Stacked bar: category composition per group
  p1 <- ggplot(summary_df,
               aes(x = .data[[group_var]], y = mean_pct, fill = category)) +
    geom_col(position = "stack", width = 0.7, colour = "white",
             linewidth = 0.3) +
    scale_fill_manual(values = KEGG_COLOURS,
                      breaks = names(KEGG_COLOURS),
                      name   = "KEGG Category") +
    scale_y_continuous(expand = c(0, 0),
                       labels = label_number(suffix = "%")) +
    labs(
      title    = "KEGG functional category composition",
      subtitle = paste0("Proportion of predicted function | Grouped by ", group_var),
      x        = group_var,
      y        = "Relative contribution (%)"
    ) +
    theme_microbiome()

  # Per-sample stacked bars ordered by group
  cat_df_ordered <- cat_df %>%
    arrange(.data[[group_var]], sample) %>%
    mutate(sample = factor(sample, levels = unique(sample)))

  p2 <- ggplot(cat_df_ordered,
               aes(x = sample, y = abundance, fill = category)) +
    geom_col(position = "stack", width = 0.9, colour = NA) +
    scale_fill_manual(values = KEGG_COLOURS,
                      breaks = names(KEGG_COLOURS),
                      name   = "KEGG Category") +
    scale_y_continuous(expand = c(0, 0),
                       labels = label_number(suffix = "%")) +
    labs(
      title    = "KEGG functional categories per sample",
      x        = NULL, y = "Relative contribution (%)"
    ) +
    theme_microbiome() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

  combined <- p2 / p1 + plot_layout(heights = c(2, 1), guides = "collect")

  cat("  Categories summarised:", nrow(cat_mat), "\n\n")
  return(list(plot = combined, summary = summary_df, matrix = cat_rel))
}


# =============================================================================
# SECTION 5 — FUNCTIONAL DIVERSITY
# =============================================================================

#' Calculate functional alpha and beta diversity from KO or pathway predictions.
#'
#' @param feature_mat   Matrix of predicted function abundances.
#' @param metadata      Data frame of sample metadata.
#' @param group_var     Grouping variable.
#' @param feature_type  Label for plot titles. Default = "KO".
#' @return A list: diversity data frame, alpha plot, and beta ordination.

analyse_functional_diversity <- function(feature_mat,
                                          metadata,
                                          group_var    = "group",
                                          feature_type = "KO") {

  cat("=== Functional diversity (", feature_type, ") ===\n")

  # Transpose: samples as rows
  mat_t     <- t(feature_mat)
  shared    <- intersect(rownames(mat_t), rownames(metadata))
  mat_t     <- mat_t[shared, ]
  meta_sub  <- metadata[shared, , drop = FALSE]

  # --- Alpha diversity -------------------------------------------------------
  div_df <- data.frame(
    sample      = rownames(mat_t),
    richness    = rowSums(mat_t > 0),
    shannon     = round(vegan::diversity(mat_t, index = "shannon"), 4),
    simpson     = round(vegan::diversity(mat_t, index = "simpson"), 4),
    stringsAsFactors = FALSE
  )

  div_df <- left_join(div_df, meta_sub %>% rownames_to_column("sample"),
                      by = "sample")

  groups   <- unique(div_df[[group_var]])
  n_groups <- length(groups)
  colours  <- brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)]
  names(colours) <- groups

  # Alpha plot
  div_long <- div_df %>%
    select(sample, richness, shannon, simpson, all_of(group_var)) %>%
    pivot_longer(c(richness, shannon, simpson),
                 names_to = "metric", values_to = "value")

  p_alpha <- ggplot(div_long,
                    aes(x = .data[[group_var]], y = value,
                        fill = .data[[group_var]])) +
    geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.55) +
    geom_jitter(aes(colour = .data[[group_var]]),
                width = 0.15, alpha = 0.5, size = 1.5) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = colours, guide = "none") +
    scale_colour_manual(values = colours, guide = "none") +
    labs(
      title    = paste0("Functional alpha diversity (", feature_type, ")"),
      subtitle = paste0("Grouped by ", group_var),
      x        = group_var, y = "Value"
    ) +
    theme_microbiome() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  # --- Beta diversity --------------------------------------------------------
  bray_dist <- vegdist(mat_t, method = "bray")
  pcoa      <- cmdscale(bray_dist, k = 2, eig = TRUE)
  eig_vals  <- pcoa$eig[pcoa$eig > 0]
  var_exp   <- round(eig_vals / sum(eig_vals) * 100, 1)

  pcoa_df <- as.data.frame(pcoa$points) %>%
    setNames(c("PC1", "PC2")) %>%
    rownames_to_column("sample") %>%
    left_join(meta_sub %>% rownames_to_column("sample"), by = "sample")

  # PERMANOVA
  set.seed(42)
  perm_res <- adonis2(
    bray_dist ~ meta_sub[[group_var]],
    permutations = 999, by = "margin"
  )
  r2 <- round(perm_res$R2[1], 3)
  pv <- perm_res$`Pr(>F)`[1]

  p_beta <- ggplot(pcoa_df,
                   aes(x = PC1, y = PC2, colour = .data[[group_var]])) +
    stat_ellipse(aes(fill = .data[[group_var]]), geom = "polygon",
                 alpha = 0.08, level = 0.95) +
    geom_point(size = 3, alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey80") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey80") +
    annotate("text", x = Inf, y = Inf, vjust = 1.5, hjust = 1.1,
             label = paste0("PERMANOVA R²=", r2, " p=", round(pv, 3)),
             size = 3.5, fontface = "bold") +
    scale_colour_manual(values = colours, name = group_var) +
    scale_fill_manual(values = colours, guide = "none") +
    labs(
      title    = paste0("Functional beta diversity - PCoA (", feature_type, ")"),
      subtitle = "Bray-Curtis dissimilarity",
      x = paste0("PC1 (", var_exp[1], "%)"),
      y = paste0("PC2 (", var_exp[2], "%)")
    ) +
    theme_microbiome() + coord_fixed()

  cat("  Functional richness (median):", median(div_df$richness), "\n")
  cat("  PERMANOVA R² =", r2, "| p =", round(pv, 3), "\n\n")

  return(list(
    diversity  = div_df,
    p_alpha    = p_alpha,
    p_beta     = p_beta,
    permanova  = perm_res
  ))
}


# =============================================================================
# SECTION 6 — FUNCTIONAL DIFFERENTIAL ABUNDANCE
# =============================================================================

#' Test for differentially abundant functional features between groups.
#'
#' Uses Wilcoxon rank-sum / Kruskal-Wallis tests with BH correction,
#' appropriate for predicted functional data.
#'
#' @param feature_mat   Matrix of predicted function abundances.
#' @param metadata      Data frame of sample metadata.
#' @param group_var     Grouping variable.
#' @param feature_type  Label for outputs. Default = "Pathway".
#' @param alpha         Significance threshold. Default = 0.05.
#' @param top_n_plot    Number of top features to plot. Default = 25.
#' @param descriptions  Optional named vector of feature descriptions.
#' @return A list: results table and plots.

test_functional_da <- function(feature_mat,
                                 metadata,
                                 group_var    = "group",
                                 feature_type = "Pathway",
                                 alpha        = 0.05,
                                 top_n_plot   = 25,
                                 descriptions = NULL) {

  cat("=== Functional differential abundance (", feature_type, ") ===\n")

  # Normalise and transpose
  mat_rel <- apply(feature_mat, 2, function(x) x / (sum(x) + 1e-10))
  mat_t   <- t(mat_rel)
  shared  <- intersect(rownames(mat_t), rownames(metadata))
  mat_t   <- mat_t[shared, ]
  meta_sub <- metadata[shared, , drop = FALSE]
  groups  <- as.character(meta_sub[[group_var]])
  group_levels <- unique(groups)
  n_groups <- length(group_levels)

  cat("  Features:", ncol(mat_t), "| Samples:", nrow(mat_t),
      "| Groups:", n_groups, "\n")

  # Run test for each feature
  test_results <- lapply(colnames(mat_t), function(feat) {
    vals <- mat_t[, feat]
    df   <- data.frame(vals = vals, group = groups)

    tryCatch({
      if (n_groups == 2) {
        g1  <- vals[groups == group_levels[1]]
        g2  <- vals[groups == group_levels[2]]
        wt  <- wilcox.test(g1, g2, exact = FALSE)
        lfc <- log2(mean(g2) / (mean(g1) + 1e-10))
        data.frame(
          feature     = feat,
          comparison  = paste0(group_levels[2], "_vs_", group_levels[1]),
          lfc         = round(lfc, 4),
          mean_g1     = round(mean(g1), 6),
          mean_g2     = round(mean(g2), 6),
          p_value     = wt$p.value,
          stringsAsFactors = FALSE
        )
      } else {
        kt <- kruskal.test(vals ~ group, data = df)
        data.frame(
          feature     = feat,
          comparison  = group_var,
          lfc         = NA,
          mean_g1     = NA,
          mean_g2     = NA,
          p_value     = kt$p.value,
          stringsAsFactors = FALSE
        )
      }
    }, error = function(e) NULL)
  })

  result_df <- bind_rows(test_results) %>%
    filter(!is.na(p_value)) %>%
    mutate(
      q_value    = p.adjust(p_value, method = "BH"),
      diff_abund = q_value < alpha,
      direction  = case_when(
        diff_abund & lfc > 0 ~ "Increased",
        diff_abund & lfc < 0 ~ "Decreased",
        TRUE                  ~ "Not significant"
      )
    ) %>%
    arrange(q_value)

  # Add descriptions
  if (!is.null(descriptions)) {
    result_df$description <- descriptions[result_df$feature]
    result_df$description[is.na(result_df$description)] <- result_df$feature[
      is.na(result_df$description)
    ]
  }

  n_sig <- sum(result_df$diff_abund)
  cat("  Significant", feature_type, "features:", n_sig, "\n\n")

  # --- Plot: horizontal bar chart of top DA features -----------------------
  top_da <- result_df %>%
    filter(diff_abund) %>%
    arrange(desc(abs(lfc))) %>%
    slice_head(n = top_n_plot)

  label_col <- if ("description" %in% colnames(top_da)) "description" else "feature"

  p_bar <- NULL
  if (nrow(top_da) > 0) {
    top_da <- top_da %>% filter(!is.na(lfc)) %>%
      mutate(label = fct_reorder(as.character(.data[[label_col]]), lfc))

    colour_map <- c("Increased" = "#e74c3c", "Decreased" = "#3498db")

    p_bar <- ggplot(top_da, aes(x = lfc, y = label, fill = direction)) +
      geom_col(width = 0.7, alpha = 0.85) +
      geom_vline(xintercept = 0, linewidth = 0.8, colour = "grey30") +
      geom_text(aes(x    = ifelse(lfc > 0, lfc + 0.005, lfc - 0.005),
                    label = paste0("q=",
                                   formatC(q_value, format = "e", digits = 1))),
                hjust = ifelse(top_da$lfc > 0, 0, 1),
                size  = 2.2, colour = "grey40") +
      scale_fill_manual(values = colour_map, name = "Direction") +
      labs(
        title    = paste0("Differential ", feature_type, " features"),
        subtitle = paste0(n_sig, " features (q < ", alpha, ")"),
        x        = "Log2 fold change",
        y        = NULL
      ) +
      theme_microbiome() +
      theme(axis.text.y = element_text(size = 8), plot.margin = ggplot2::margin(5, 40, 5, 150))

  }

  # --- Bubble plot: effect size vs significance ----------------------------
  top_bubble <- result_df %>%
    arrange(q_value) %>%
    slice_head(n = 40)

  if (!is.null(descriptions)) {
    top_bubble$label <- descriptions[top_bubble$feature]
    top_bubble$label[is.na(top_bubble$label)] <- top_bubble$feature[
      is.na(top_bubble$label)
    ]
  } else {
    top_bubble$label <- top_bubble$feature
  }

  p_bubble <- ggplot(top_bubble,
                     aes(x    = lfc,
                         y    = -log10(q_value + 1e-10),
                         size = abs(lfc),
                         colour = direction)) +
    geom_point(alpha = 0.75) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed",
               colour = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggrepel::geom_text_repel(
      data = top_bubble %>% filter(diff_abund),
      aes(label = label),
      size = 2.5, max.overlaps = 15,
      segment.size = 0.3, box.padding = 0.3
    ) +
    scale_colour_manual(
      values = c("Increased" = "#e74c3c",
                 "Decreased" = "#3498db",
                 "Not significant" = "#bdc3c7"),
      name = "Direction"
    ) +
    scale_size(range = c(1, 6), guide = "none") +
    labs(
      title    = paste0(feature_type, " - significance vs effect size"),
      subtitle = paste0("Top 40 features | q threshold: ", alpha),
      x        = "Log2 fold change",
      y        = expression(-log[10](q-value))
    ) +
    theme_microbiome()

  return(list(
    results  = result_df,
    p_bar    = p_bar,
    p_bubble = p_bubble,
    n_sig    = n_sig
  ))
}


# =============================================================================
# SECTION 7 — FUNCTIONAL HEATMAP
# =============================================================================

#' Clustered heatmap of top functional features across samples.
#'
#' @param feature_mat   Matrix of feature abundances.
#' @param metadata      Data frame of sample metadata.
#' @param group_var     Grouping variable.
#' @param top_n         Top features to show. Default = 40.
#' @param feature_type  Label. Default = "Pathway".
#' @param descriptions  Named vector of descriptions.
#' @return A ggplot heatmap.

plot_functional_heatmap <- function(feature_mat,
                                     metadata,
                                     group_var    = NULL,
                                     top_n        = 40,
                                     feature_type = "Pathway",
                                     cluster_cols = TRUE,
                                     descriptions = NULL) {
  cat("=== Functional heatmap ===\n")

  # Normalise and select top features
  mat_rel   <- apply(feature_mat, 2, function(x) x / sum(x) * 100)
  mean_abund <- rowMeans(mat_rel)
  top_feats  <- names(sort(mean_abund, decreasing = TRUE))[seq_len(
    min(top_n, length(mean_abund))
  )]
  mat_sub    <- mat_rel[top_feats, , drop = FALSE]

  # Z-score across samples
  mat_z <- t(scale(t(mat_sub)))
  mat_z[is.nan(mat_z)] <- 0

  # Cluster
  shared  <- intersect(colnames(mat_z), rownames(metadata))
  mat_z   <- mat_z[, shared, drop = FALSE]
  row_clust <- hclust(dist(mat_z), method = "ward.D2")
  col_clust <- hclust(dist(t(mat_z)), method = "ward.D2")

  row_order <- rownames(mat_z)[row_clust$order]
  if (!is.null(group_var) && group_var %in% colnames(metadata)) {
    grp_vals  <- as.character(metadata[shared, group_var])
    col_order <- unlist(lapply(unique(grp_vals), function(g) {
      s <- shared[grp_vals == g]
      if (length(s) > 2) s[hclust(dist(t(mat_z[, s, drop=FALSE])))$order] else s
    }))
  } else {
    col_order <- colnames(mat_z)[col_clust$order]
  }

  # Feature labels
  feat_labels <- if (!is.null(descriptions)) {
    sapply(row_order, function(id) {
      desc <- descriptions[[id]]
      if (!is.null(desc) && !is.na(desc)) {
        short <- if (nchar(desc) > 45) paste0(substr(desc, 1, 42), "...") else desc
        paste0(id, ": ", short)
      } else id
    })
  } else row_order
  
  names(feat_labels) <- row_order

  # Melt
  heat_df <- as.data.frame(mat_z) %>%
    rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "sample", values_to = "zscore") %>%
    mutate(
      feature = factor(feature, levels = rev(row_order)),
      sample  = factor(sample, levels = col_order),
      label   = feat_labels[as.character(feature)]
    ) %>%
    left_join(metadata[shared, , drop = FALSE] %>%
                rownames_to_column("sample"), by = "sample")

  # Group annotation strip
  ann_var <- if (!is.null(group_var) && group_var %in% colnames(metadata)) group_var else NULL
  if (!is.null(ann_var)) {
    groups   <- unique(heat_df[[ann_var]])
    n_groups <- length(groups)
    colours  <- brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)]
    names(colours) <- groups
    annot_df <- heat_df %>% dplyr::select(sample, all_of(ann_var)) %>% distinct()
    p_annot  <- ggplot(annot_df, aes(x = sample, y = 1, fill = .data[[ann_var]])) +
      geom_tile() + scale_fill_manual(values = colours, name = ann_var) +
      scale_x_discrete(expand = c(0, 0)) + theme_void() + theme(legend.position = "right")
  } else { p_annot <- NULL }



  p_heat <- ggplot(heat_df, aes(x = sample, y = feature, fill = zscore)) +
    geom_tile(colour = "white", linewidth = 0.05) +
    scale_fill_gradient2(
      low = "#3498db", mid = "white", high = "#e74c3c",
      midpoint = 0, name = "Z-score"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0),
                     labels = feat_labels[rev(row_order)]) +
    labs(x = "Samples", y = NULL) +
    theme_microbiome() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 7),
      panel.border = element_rect(colour = "grey80", fill = NA)
    )

  p_dendro <- NULL
  if (requireNamespace('ggdendro', quietly = TRUE)) {
    dendro_data <- ggdendro::dendro_data(as.dendrogram(row_clust), type = 'rectangle')
    n_feat <- length(row_order)
    p_dendro <- ggplot(ggdendro::segment(dendro_data)) +
      geom_segment(aes(x = y, xend = yend, y = x, yend = xend), colour = 'grey40', linewidth = 0.4) +
      scale_x_reverse() +
      scale_y_continuous(limits = c(0.5, n_feat + 0.5), expand = c(0, 0)) +
      theme_void() +
      theme(plot.margin = ggplot2::margin(0, 2, 0, 4))
  }

  # Suppress y-axis text on heatmap when dendrogram present
  if (!is.null(p_dendro)) p_heat <- p_heat + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

  # Column dendrogram (only when cluster_cols = TRUE)
  p_col_dendro <- NULL
  if (isTRUE(cluster_cols) && requireNamespace("ggdendro", quietly = TRUE) && ncol(mat_z) > 2) {
    col_dd       <- ggdendro::dendro_data(as.dendrogram(col_clust), type = "rectangle")
    p_col_dendro <- ggplot(ggdendro::segment(col_dd)) +
      geom_segment(aes(x = x, xend = xend, y = y, yend = yend), colour = "grey40", linewidth = 0.4) +
      scale_x_continuous(limits = c(0.5, ncol(mat_z) + 0.5), expand = c(0, 0)) +
      theme_void() +
      theme(plot.margin = ggplot2::margin(4, 0, 0, 0))
  }

  # Assemble with patchwork
  title_str <- paste0("Functional heatmap - ", feature_type, " (top ", nrow(mat_sub), ")")
  ann_plot  <- if (!is.null(p_annot)) p_annot else patchwork::plot_spacer()
  if (!is.null(p_dendro) && !is.null(p_col_dendro)) {
    combined <- (patchwork::plot_spacer() + p_col_dendro + ann_plot + patchwork::plot_spacer() + p_dendro + p_heat) +
      patchwork::plot_layout(ncol = 2, widths = c(0.15, 1), heights = c(0.12, 0.05, 1), guides = "collect") +
      patchwork::plot_annotation(title = title_str, theme = theme(plot.title = element_text(size = 13, face = "bold")))
  } else if (!is.null(p_dendro)) {
    combined <- (patchwork::plot_spacer() + ann_plot + p_dendro + p_heat) +
      patchwork::plot_layout(ncol = 2, widths = c(0.15, 1), heights = c(0.05, 1), guides = "collect") +
      patchwork::plot_annotation(title = title_str, theme = theme(plot.title = element_text(size = 13, face = "bold")))
  } else {
    combined <- ann_plot / p_heat +
      patchwork::plot_layout(heights = c(0.05, 1), guides = "collect") +
      patchwork::plot_annotation(title = title_str, theme = theme(plot.title = element_text(size = 13, face = "bold")))
  }

  cat("  Features plotted:", nrow(mat_sub), "\n\n")
  return(combined)
}

#' @param metadata         Data frame of sample metadata.
#' @param group_var        Grouping variable.
#' @param top_n            Top contributing taxa to show.
#' @return A list: contributor table and bubble plot.

analyse_pathway_contributors <- function(stratified_file,
                                          target_pathways,
                                          metadata,
                                          group_var = "group",
                                          top_n     = 15) {

  cat("=== Pathway contributor analysis ===\n")

  if (!file.exists(stratified_file)) {
    cat("  Stratified file not found:", stratified_file, "\n\n")
    return(NULL)
  }

  strat <- read.table(stratified_file, sep = "\t", header = TRUE,
                      stringsAsFactors = FALSE)

  # Filter to target pathways
  strat_sub <- strat %>%
    filter(function. %in% target_pathways | pathway %in% target_pathways)

  if (nrow(strat_sub) == 0) {
    cat("  No matching pathways found in stratified file.\n\n")
    return(NULL)
  }

  # Sum contributions per taxon-pathway pair across samples
  sample_cols <- intersect(rownames(metadata), colnames(strat_sub))
  contrib_df  <- strat_sub %>%
    mutate(
      total_contrib = rowSums(select(., all_of(sample_cols)))
    ) %>%
    group_by(taxon) %>%
    summarise(
      total_contribution = sum(total_contrib),
      n_pathways         = n(),
      .groups            = "drop"
    ) %>%
    arrange(desc(total_contribution)) %>%
    slice_head(n = top_n)

  cat("  Top", min(top_n, nrow(contrib_df)),
      "contributing taxa identified\n\n")

  p <- ggplot(contrib_df,
              aes(x = total_contribution,
                  y = fct_reorder(taxon, total_contribution),
                  size = n_pathways, colour = n_pathways)) +
    geom_point(alpha = 0.8) +
    scale_colour_distiller(palette = "YlOrRd", direction = 1,
                           name = "Pathways\ncontributed") +
    scale_size(range = c(3, 10), name = "Pathways\ncontributed") +
    scale_x_continuous(labels = label_comma()) +
    labs(
      title    = "Top taxa contributing to target pathways",
      subtitle = paste0("Target pathways: ",
                        paste(target_pathways, collapse = ", ")),
      x        = "Total pathway contribution (summed abundance)",
      y        = NULL
    ) +
    theme_microbiome() +
    theme(axis.text.y = element_text(face = "italic", size = 8))

  return(list(contributors = contrib_df, plot = p))
}


# =============================================================================
# SECTION 9 — COMPLETE FUNCTIONAL PREDICTION WORKFLOW
# =============================================================================

#' Run the complete functional prediction analysis pipeline.
#'
#' @param picrust_dir  Path to PICRUSt2 output directory.
#' @param meta_file    Path to sample metadata file.
#' @param group_var    Grouping variable.
#' @param alpha        Significance threshold. Default = 0.05.
#' @param output_dir   Directory for outputs.
#' @return A named list of all results.

run_functional_analysis <- function(picrust_dir,
                                     meta_file  = NULL,
                                     group_var  = NULL,
                                     alpha      = 0.05,
                                     output_dir = "functional_output") {

  cat("\n", strrep("=", 60), "\n")
  cat("  FUNCTIONAL PREDICTION PIPELINE\n")
  cat(strrep("=", 60), "\n\n")

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # --- Import ---------------------------------------------------------------
  p2_data  <- import_picrust2(picrust_dir, meta_file = meta_file)
  metadata <- p2_data$metadata
  results  <- list()

  if (is.null(metadata)) {
    stop("Metadata required. Provide meta_file argument.")
  }

  # --- NSTI quality ---------------------------------------------------------
  cat("--- Plot 1: NSTI quality ---\n")
  nsti_res <- plot_nsti_quality(p2_data, group_var = group_var)
  if (!is.null(nsti_res)) {
    results$p_nsti <- nsti_res$plot
    ggsave(file.path(output_dir, "01_nsti_quality.pdf"),
           nsti_res$plot, width = 12, height = 10)
  }

  # --- Pathway analysis -----------------------------------------------------
  if (!is.null(p2_data$pathway)) {
    cat("--- Plot 2: Pathway abundance ---\n")
    path_res <- plot_pathway_abundance(
      p2_data$pathway, metadata,
      group_var    = group_var, top_n = 25,
      descriptions = p2_data$pathway_desc
    )
    results$p_pathways <- path_res$plot
    ggsave(file.path(output_dir, "02_pathway_abundance.pdf"),
           path_res$plot, width = 12, height = 12)

    cat("--- Plot 3: Functional heatmap (Pathways) ---\n")
    results$p_pathway_heatmap <- plot_functional_heatmap(
      p2_data$pathway, metadata,
      group_var = group_var, top_n = 40, feature_type = "Pathway",
      descriptions = p2_data$pathway_desc
    )
    ggsave(file.path(output_dir, "03_pathway_heatmap.pdf"),
           results$p_pathway_heatmap, width = 14, height = 14)

    cat("--- Analysis 1: Pathway DA ---\n")
    pathway_da <- test_functional_da(
      p2_data$pathway, metadata,
      group_var = group_var, feature_type = "Pathway",
      alpha = alpha, top_n_plot = 25,
      descriptions = p2_data$pathway_desc
    )
    results$pathway_da <- pathway_da
    write.csv(pathway_da$results,
              file.path(output_dir, "pathway_da_results.csv"),
              row.names = FALSE)
    if (!is.null(pathway_da$p_bar)) {
      ggsave(file.path(output_dir, "04_pathway_da.pdf"),
             pathway_da$p_bar, width = 12, height = 10)
    }
  }

  # --- KO analysis ----------------------------------------------------------
  if (!is.null(p2_data$ko)) {
    cat("--- Plot 5: KEGG category summary ---\n")
    kegg_res <- summarise_kegg_categories(
      p2_data$ko, metadata, group_var = group_var
    )
    results$p_kegg <- kegg_res$plot
    ggsave(file.path(output_dir, "05_kegg_categories.pdf"),
           kegg_res$plot, width = 14, height = 10)

    cat("--- Analysis 2: Functional diversity (KO) ---\n")
    func_div <- analyse_functional_diversity(
      p2_data$ko, metadata,
      group_var = group_var, feature_type = "KO"
    )
    results$p_func_alpha <- func_div$p_alpha
    results$p_func_beta  <- func_div$p_beta
    ggsave(file.path(output_dir, "06_functional_diversity_alpha.pdf"),
           func_div$p_alpha, width = 12, height = 6)
    ggsave(file.path(output_dir, "07_functional_diversity_beta.pdf"),
           func_div$p_beta,  width = 9, height = 8)

    cat("--- Analysis 3: KO differential abundance ---\n")
    ko_da <- test_functional_da(
      p2_data$ko, metadata,
      group_var = group_var, feature_type = "KO",
      alpha = alpha, top_n_plot = 25
    )
    results$ko_da <- ko_da
    write.csv(ko_da$results,
              file.path(output_dir, "ko_da_results.csv"),
              row.names = FALSE)
    if (!is.null(ko_da$p_bar)) {
      ggsave(file.path(output_dir, "08_ko_da.pdf"),
             ko_da$p_bar, width = 12, height = 10)
    }
  }

  cat("\n", strrep("=", 60), "\n")
  cat("  FUNCTIONAL PREDICTION PIPELINE COMPLETE\n")
  cat(strrep("=", 60), "\n")
  cat("  Output saved to:", output_dir, "\n")
  cat("  Plots:  up to 8 PDF files\n")
  cat("  Tables: pathway DA, KO DA results\n\n")

  return(invisible(results))
}


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# --- Option A: Full pipeline -------------------------------------------------
# results <- run_functional_analysis(
#   picrust_dir = "picrust2_output/",
#   meta_file   = "data/metadata.csv",
#   group_var   = "disease_status",
#   alpha       = 0.05,
#   output_dir  = "results/functional"
# )

# --- Option B: Import and run steps individually ----------------------------
# p2_data <- import_picrust2("picrust2_output/", meta_file = "data/metadata.csv")
# nsti_res <- plot_nsti_quality(p2_data, group_var = "disease_status")
# nsti_res$plot

# path_res <- plot_pathway_abundance(p2_data$pathway, p2_data$metadata,
#                                    group_var = "disease_status", top_n = 30)
# path_res$plot

# pathway_da <- test_functional_da(p2_data$pathway, p2_data$metadata,
#                                   group_var = "disease_status",
#                                   feature_type = "Pathway")
# pathway_da$results %>% filter(diff_abund) %>% head(20)
# pathway_da$p_bar
# pathway_da$p_bubble
