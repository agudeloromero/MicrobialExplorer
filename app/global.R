# =============================================================================
# global.R — MicrobialExplorer Shiny App
# =============================================================================
# Loaded once when the app starts. Runs in the global environment, so
# everything defined here is shared across all user sessions.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  # Shiny
  library(shiny)
  library(shinyjs)
  library(bslib)
  library(DT)
  library(bsplus)
  library(shinycssloaders)

  # Phyloseq / biology
  library(phyloseq)
  library(ape)

  # Tidyverse
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(forcats)

  # Plotting
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(patchwork)
  library(ggdendro)

  # Stats / ecology
  library(vegan)
  library(scales)

  # Suppress noisy tidytree conflict BEFORE loading anything that triggers it
  suppressMessages(suppressWarnings({
    if (requireNamespace("tidytree", quietly = TRUE)) library(tidytree)
  }))
})

# ── Source R pipeline modules (relative to app/ directory) ───────────────────
r_dir <- file.path(dirname(getwd()), "R")
if (!dir.exists(r_dir)) r_dir <- "../R"   # fallback when running from app/

for (f in list.files(r_dir, pattern = "\\.R$", full.names = TRUE)) {
  source(f)
}

# ── Source Shiny modules from app/modules/ ────────────────────────────────────
mod_dir <- file.path(getwd(), "modules")
if (dir.exists(mod_dir)) {
  for (f in list.files(mod_dir, pattern = "\\.R$", full.names = TRUE)) {
    source(f)
  }
  message("[global.R] Sourced modules from: ", mod_dir)
} else {
  warning("[global.R] Could not find modules directory: ", mod_dir)
}

# ── Optional package availability flags ──────────────────────────────────────
ANCOMBC_AVAILABLE  <- requireNamespace("ANCOMBC",     quietly = TRUE)
PICANTE_AVAILABLE  <- requireNamespace("picante",     quietly = TRUE)
SPIECEASI_AVAILABLE <- requireNamespace("SpiecEasi",  quietly = TRUE)

if (!ANCOMBC_AVAILABLE)
  message("[global.R] ANCOMBC not available — DA will use DESeq2 + ALDEx2 consensus")
if (!PICANTE_AVAILABLE)
  message("[global.R] picante not available — Faith's PD will be skipped")

# ── App-wide constants ────────────────────────────────────────────────────────
APP_TITLE   <- "MicrobialExplorer"
APP_VERSION <- "v1.0"
MAX_UPLOAD_MB <- 200

options(shiny.maxRequestSize = MAX_UPLOAD_MB * 1024^2)

# Accepted file extensions for the four upload slots
ACCEPTED_OTU      <- c(".csv", ".tsv", ".txt", ".biom")
ACCEPTED_TAX      <- c(".csv", ".tsv", ".txt")
ACCEPTED_META     <- c(".csv", ".tsv", ".txt")
ACCEPTED_TREE     <- c(".nwk", ".newick", ".tre", ".tree", ".nex")
ACCEPTED_PICRUST  <- c(".tsv", ".tsv.gz")

# ── Data type detection ───────────────────────────────────────────────────────

#' Detect whether OTU values are counts, relative abundances, or normalised
#'
#' @param otu_mat  A numeric matrix (taxa × samples)
#' @return One of "counts", "relative", or "normalised"
detect_data_type <- function(otu_mat) {
  # Force to plain base R matrix — phyloseq otu_table objects
  # do not support logical matrix subsetting
  m    <- matrix(as.numeric(otu_mat), nrow = nrow(otu_mat), ncol = ncol(otu_mat))
  vals <- m[m > 0]
  if (length(vals) == 0) return("unknown")
  if (all(vals == floor(vals))) return("counts")
  if (max(vals) <= 1 && all(colSums(m) <= 1.01)) return("relative")
  return("normalised")
}

#' Detect features of a phyloseq object used to enable / disable modules
#'
#' @param ps  A phyloseq object
#' @return A named list of dataset features
detect_features <- function(ps) {
  raw     <- otu_table(ps)
  otu_mat <- matrix(as.numeric(raw), nrow = nrow(raw), ncol = ncol(raw))
  if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)

  meta_df     <- tryCatch(data.frame(sample_data(ps)), error = function(e) NULL)
  numeric_vars <- if (!is.null(meta_df))
    names(which(sapply(meta_df, is.numeric))) else character(0)
  cat_vars     <- if (!is.null(meta_df))
    names(which(sapply(meta_df, function(x) is.character(x) || is.factor(x)))) else character(0)
  sample_vars  <- if (!is.null(meta_df)) colnames(meta_df) else character(0)

  avail_ranks <- tryCatch(rank_names(ps), error = function(e) character(0))

  has_tree <- tryCatch({
    !is.null(phy_tree(ps, errorIfNULL = FALSE))
  }, error = function(e) FALSE)

  has_taxonomy <- tryCatch({
    !is.null(tax_table(ps, errorIfNULL = FALSE))
  }, error = function(e) FALSE)

  data_type <- detect_data_type(otu_mat)

  # Determine which group variables have >= 2 levels
  group_vars <- if (!is.null(meta_df)) {
    Filter(function(v) {
      n_distinct(meta_df[[v]], na.rm = TRUE) >= 2
    }, cat_vars)
  } else character(0)

  binary_group_vars <- if (!is.null(meta_df)) {
    Filter(function(v) {
      n_distinct(meta_df[[v]], na.rm = TRUE) == 2
    }, cat_vars)
  } else character(0)

  list(
    data_type         = data_type,
    has_taxonomy      = has_taxonomy,
    has_tree          = has_tree,
    avail_ranks       = avail_ranks,
    n_samples         = nsamples(ps),
    n_taxa            = ntaxa(ps),
    sample_vars       = sample_vars,
    numeric_vars      = numeric_vars,
    cat_vars          = cat_vars,
    group_vars        = group_vars,
    binary_group_vars = binary_group_vars
  )
}

# ── Module availability logic ─────────────────────────────────────────────────

#' Return a named logical vector of which analysis modules are available
#'
#' @param features  Output of detect_features()
#' @return Named logical vector
get_module_availability <- function(features) {

  is_counts    <- features$data_type == "counts"
  has_tax      <- features$has_taxonomy
  has_groups   <- length(features$group_vars) >= 1
  has_binary   <- length(features$binary_group_vars) >= 1
  has_longit   <- all(c("timepoint", "subject_id") %in% features$sample_vars) ||
    any(grepl("time|visit|week|day|month", features$sample_vars, ignore.case = TRUE))

  c(
    qc           = TRUE,
    composition  = has_tax,
    alpha        = is_counts,
    beta         = TRUE,
    da           = is_counts && has_groups,
    functional   = TRUE,           # gracefully informs user if no PICRUSt2 data
    network      = TRUE,
    ml           = is_counts && has_binary,
    correlation  = TRUE,
    longitudinal = has_longit
  )
}

# ── Colour helpers ────────────────────────────────────────────────────────────
# These mirror TAXA_PALETTE in utils.R for use in Shiny UI elements
ME_BLUE  <- "#2c3e50"
ME_TEAL  <- "#1abc9c"
ME_RED   <- "#e74c3c"
ME_GREY  <- "#bdc3c7"

# ── Read helper for the four supported delimited formats ─────────────────────

#' Read a delimited OTU/taxonomy/metadata file, handling .csv / .tsv / .txt
#'
#' @param path    File path
#' @param header  Logical — first row is column names
#' @return data.frame

read_table_file <- function(path, header = TRUE) {
  raw_lines <- readLines(path, n = 2, warn = FALSE)
  # Detect separator from content, not extension (Shiny temp files have no ext)
  n_comma <- lengths(regmatches(raw_lines[2], gregexpr(",",  raw_lines[2])))
  n_tab   <- lengths(regmatches(raw_lines[2], gregexpr("\t", raw_lines[2])))
  sep <- if (n_comma >= n_tab) "," else "\t"
  # Detect row-names column: blank/quoted-blank first header field
  first_field <- gsub('"', "", trimws(strsplit(raw_lines[1], sep, fixed=TRUE)[[1]][1]))
  header_n <- length(strsplit(raw_lines[1], sep, fixed=TRUE)[[1]])
  data_n   <- length(strsplit(raw_lines[2], sep, fixed=TRUE)[[1]])
  has_rownames_col <- (first_field == "") || (data_n == header_n + 1)
  df <- suppressWarnings(
    read.table(path, sep=sep, header=header, row.names=NULL,
               check.names=FALSE, comment.char="", quote="")
  )
  if (has_rownames_col) {
    rownames(df) <- make.unique(as.character(df[[1]]), sep="_dup")
  df[[1]] <- NULL
  # Strip residual quote characters from names and values (artifact of quote="")
  colnames(df) <- gsub('"', "", colnames(df))
  rownames(df) <- gsub('"', "", rownames(df))
  df[] <- lapply(df, function(x) if(is.character(x)) gsub('"', "", x) else x)
  # Strip residual quote characters from values (artifact of quote="")
  df[] <- lapply(df, function(x) if(is.character(x)) gsub('"', "", x) else x)
  }
  df
}

# ── Demo / example data loader ────────────────────────────────────────────────

#' Load a built-in demo dataset from phyloseq
#'
#' @param name  One of "GlobalPatterns", "enterotype", "soilrep"
#' @return phyloseq object
load_demo_data <- function(name = "GlobalPatterns") {
  if (!requireNamespace("phyloseq", quietly = TRUE))
    stop("phyloseq required for demo data")
  switch(name,
    GlobalPatterns = {
      data("GlobalPatterns", package = "phyloseq", envir = environment())
      get("GlobalPatterns", envir = environment())
    },
    enterotype = {
      data("enterotype", package = "phyloseq", envir = environment())
      get("enterotype", envir = environment())
    },
    soilrep = {
      data("soilrep", package = "phyloseq", envir = environment())
      get("soilrep", envir = environment())
    },
    network_enriched = {
      d <- file.path(dirname(getwd()), "data", "simulated")
      otu  <- read_table_file(file.path(d, "network_enriched_otu.csv"))
      tax  <- read_table_file(file.path(d, "network_enriched_taxonomy.csv"))
      meta <- read_table_file(file.path(d, "network_enriched_metadata.csv"))
      tree <- ape::read.tree(file.path(d, "network_enriched_tree.nwk"))
      OTU  <- otu_table(as.matrix(otu), taxa_are_rows = TRUE)
      TAX  <- tax_table(as.matrix(tax))
      SAM  <- sample_data(meta)
      phyloseq(OTU, TAX, SAM, phy_tree(tree))
    },
        stop("Unknown demo dataset: ", name)
  )
}

# ── UI helpers ────────────────────────────────────────────────────────────────

#' Create a styled info/status badge for the sidebar data summary
#'
#' @param label  Text label
#' @param value  Value to display
#' @param colour  Bootstrap colour class (default "primary")
stat_badge <- function(label, value, colour = "primary") {
  tags$div(
    class = "d-flex justify-content-between align-items-center mb-1",
    tags$span(class = "text-muted small", label),
    tags$span(class = paste0("badge bg-", colour), value)
  )
}

#' Wrap a plot output with a spinner while it is computing
#'
#' @param output_id  The plotOutput id
#' @param ...        Additional args passed to plotOutput
spinner_plot <- function(output_id, ...) {
  shinycssloaders_available <- requireNamespace("shinycssloaders", quietly = TRUE)
  p <- plotOutput(output_id, ...)
  if (shinycssloaders_available)
    shinycssloaders::withSpinner(p, type = 6, color = ME_BLUE)
  else
    p
}

message("[global.R] MicrobialExplorer global environment loaded.")
