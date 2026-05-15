# Suppress harmless phylo/tidytree cache warning
suppressMessages(suppressWarnings({
  if (requireNamespace("tidytree", quietly = TRUE)) library(tidytree)
}))

# =============================================================================
# utils.R — Shared utilities for all ViromeAnalyst modules
# =============================================================================
# Source this file at the TOP of every module script:
#   source("R/utils.R")
#
# Contains: shared theme, colour palettes, helper functions, operators
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
})

# =============================================================================
# SHARED GGPLOT2 THEME
# =============================================================================

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

theme_network <- function() {
  theme_void() +
    theme(
      legend.title    = element_text(size = 10, face = "bold"),
      legend.text     = element_text(size = 9),
      plot.title      = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle   = element_text(size = 10, colour = "grey40", hjust = 0.5),
      plot.caption    = element_text(size = 8, colour = "grey60", hjust = 1),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

# =============================================================================
# COLOUR PALETTES
# =============================================================================

# 21 colours for taxonomic bar plots ("Other" always last)
TAXA_PALETTE <- c(
  "#3498db", "#e74c3c", "#2ecc71", "#f39c12", "#9b59b6",
  "#1abc9c", "#e67e22", "#34495e", "#e91e63", "#00bcd4",
  "#8bc34a", "#ff5722", "#607d8b", "#795548", "#ffc107",
  "#673ab7", "#009688", "#ff9800", "#4caf50", "#f44336",
  "#b0bec5"   # "Other" — always last
)

# KEGG Level 1 category colours
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
# HELPER FUNCTIONS
# =============================================================================

#' Check if a package is available without loading it
pkg_available <- function(pkg) requireNamespace(pkg, quietly = TRUE)

#' Null-coalescing operator: return a if not null, else b
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Run an expression safely with timing and error catching
#'
#' @param name   Label shown in console output
#' @param expr   Expression to evaluate (use quote())
#' @return Result of expression, or NULL if an error occurred
safe_run <- function(name, expr) {
  t_start <- proc.time()["elapsed"]
  cat("[START]", name, "\n")

  result <- tryCatch({
    eval(expr)
  }, error = function(e) {
    cat("[FAIL] ", name, ":", conditionMessage(e), "\n")
    return(NULL)
  }, warning = function(w) {
    cat("[WARN] ", name, ":", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  })

  t_secs <- round(proc.time()["elapsed"] - t_start, 1)
  if (!is.null(result)) cat("[OK]   ", name, "(", t_secs, "s)\n")

  invisible(result)
}

#' Print a diagnostic summary of a phyloseq object
#'
#' @param ps    A phyloseq object
#' @param label A descriptive label for the output header
check_phyloseq <- function(ps, label = "phyloseq object") {
  cat("\n===", label, "===\n")
  cat("  Samples:", nsamples(ps), "\n")
  cat("  Taxa:   ", ntaxa(ps), "\n")
  cat("  Ranks:  ", paste(rank_names(ps), collapse = ", "), "\n")
  cat("  Total reads:", format(sum(sample_sums(ps)), big.mark = ","), "\n")
  cat("  Min reads:  ", format(min(sample_sums(ps)), big.mark = ","), "\n")
  cat("  taxa_are_rows:", taxa_are_rows(ps), "\n")

  meta_df   <- data.frame(sample_data(ps))
  num_vars  <- names(meta_df)[sapply(meta_df, is.numeric)]
  cat_vars  <- names(meta_df)[sapply(meta_df, function(x) is.character(x) | is.factor(x))]
  cat("  Numeric metadata: ", paste(num_vars, collapse = ", "), "\n")
  cat("  Group/factor vars:", paste(cat_vars,  collapse = ", "), "\n")

  # Check for any issues
  issues <- c()
  if (any(is.na(as.matrix(otu_table(ps))))) issues <- c(issues, "NA values in OTU table")
  if (any(as.matrix(otu_table(ps)) < 0))   issues <- c(issues, "Negative values in OTU table")
  if (min(sample_sums(ps)) == 0)            issues <- c(issues, "Sample(s) with 0 reads")

  if (length(issues) > 0) {
    cat("  ⚠ Issues:\n")
    for (i in issues) cat("    -", i, "\n")
  } else {
    cat("  ✓ No issues detected\n")
  }
  cat("\n")

  invisible(ps)
}

#' Standardise sample names across OTU table and metadata
#'
#' @param ps A phyloseq object
#' @return A phyloseq object with cleaned sample names
clean_sample_names <- function(ps) {
  clean <- function(x) trimws(x)
  new_names <- clean(sample_names(ps))
  sample_names(ps) <- new_names
  rownames(sample_data(ps)) <- new_names
  ps
}

#' Get the top N taxa by mean relative abundance
#'
#' @param ps    A phyloseq object
#' @param top_n Number of top taxa
#' @return Character vector of taxon names
get_top_taxa <- function(ps, top_n = 20) {
  ps_rel     <- transform_sample_counts(ps, function(x) x / sum(x))
  otu_mat    <- as.matrix(otu_table(ps_rel))
  if (!taxa_are_rows(ps_rel)) otu_mat <- t(otu_mat)
  mean_abund <- rowMeans(otu_mat)
  names(sort(mean_abund, decreasing = TRUE))[seq_len(min(top_n, length(mean_abund)))]
}

#' Apply CLR transformation safely (handles zeros with pseudocount)
#'
#' @param mat   A matrix (taxa × samples)
#' @param pseudo Pseudocount to add before log transform. Default = 0.5
#' @return CLR-transformed matrix (samples × taxa)
clr_transform <- function(mat, pseudo = 0.5) {
  mat_pseudo <- mat + pseudo
  t(apply(mat_pseudo, 2, function(x) log(x) - mean(log(x))))
}

#' Log memory usage to console
log_memory <- function(label = "") {
  if (pkg_available("pryr")) {
    mem_mb <- round(as.numeric(pryr::mem_used()) / 1e6, 1)
    cat(sprintf("[MEM] %s: %.1f MB in use\n", label, mem_mb))
  }
}

cat("[utils.R] Shared utilities loaded.\n")
# ANCOMBC availability check — graceful fallback to DESeq2 + ALDEx2
ANCOMBC_AVAILABLE <- tryCatch({
  suppressPackageStartupMessages(library(ANCOMBC))
  TRUE
}, error = function(e) {
  FALSE
})

if (!ANCOMBC_AVAILABLE) {
  message("[utils.R] ANCOMBC not available — DA analysis will use DESeq2 + ALDEx2")
}


# Safe metadata extractor — handles newer dplyr/phyloseq compatibility
get_meta <- function(ps) {
  as.data.frame(as.matrix(sample_data(ps)))
}

# Safe sample_data conversion
meta_df_safe <- function(ps) {
  df <- data.frame(sample_data(ps))
  return(df)
}


# Safe metadata extractor — handles newer dplyr/phyloseq compatibility
get_meta <- function(ps) {
  as.data.frame(as.matrix(sample_data(ps)))
}

# Safe sample_data conversion
meta_df_safe <- function(ps) {
  df <- data.frame(sample_data(ps))
  return(df)
}
