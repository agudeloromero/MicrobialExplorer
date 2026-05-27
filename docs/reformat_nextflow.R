#!/usr/bin/env Rscript
# =============================================================================
# reformat_nextflow.R
# Reformats Nextflow 16S pipeline output for upload into MicrobialExplorer
#
# Usage:
#   Rscript reformat_nextflow.R \
#     --otu   path/to/best_tax_merged_freq_tax.tsv \
#     --tax   path/to/taxonomy.tsv \
#     --outdir path/to/output_folder
#
# Output:
#   otu_table_clean.tsv   — OTU/ASV table ready for MicrobialExplorer
#   taxonomy_ranks.tsv    — Taxonomy table with separate rank columns
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# --- Parse arguments ----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, args) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx == length(args)) return(NULL)
  args[idx + 1]
}

otu_path  <- get_arg("--otu",    args)
tax_path  <- get_arg("--tax",    args)
out_dir   <- get_arg("--outdir", args)

if (is.null(otu_path) || is.null(tax_path) || is.null(out_dir)) {
  cat("Usage:\n")
  cat("  Rscript reformat_nextflow.R \\\n")
  cat("    --otu   path/to/best_tax_merged_freq_tax.tsv \\\n")
  cat("    --tax   path/to/taxonomy.tsv \\\n")
  cat("    --outdir path/to/output_folder\n")
  quit(status = 1)
}

if (!file.exists(otu_path)) stop("OTU file not found: ", otu_path)
if (!file.exists(tax_path)) stop("Taxonomy file not found: ", tax_path)
if (!dir.exists(out_dir))   dir.create(out_dir, recursive = TRUE)

# =============================================================================
# 1. OTU / ASV TABLE
# =============================================================================

cat("=== Processing OTU table ===\n")

otu_raw <- read.table(
  otu_path,
  sep              = "\t",
  header           = TRUE,
  comment.char     = "#",   # removes the #q2:types row
  stringsAsFactors = FALSE,
  check.names      = FALSE
)

# Remove non-count columns
cols_to_remove <- c("Sequence", "Taxon", "Confidence")
otu_counts <- otu_raw[, !colnames(otu_raw) %in% cols_to_remove]

# Clean sample names: remove .fastq suffix
colnames(otu_counts) <- gsub("\\.fastq$", "", colnames(otu_counts))

# Remove control samples
control_patterns <- c("Control", "control", "CONTROL", "blank", "Blank", "BLANK",
                       "negative", "Negative", "reagent", "Reagent")
is_control <- grepl(paste(control_patterns, collapse = "|"), colnames(otu_counts))
if (any(is_control)) {
  cat("  Removing control samples:", paste(colnames(otu_counts)[is_control], collapse = ", "), "\n")
  otu_counts <- otu_counts[, !is_control]
}

# Remove all-zero rows
row_sums <- rowSums(otu_counts[, -1])
n_zero   <- sum(row_sums == 0)
if (n_zero > 0) {
  cat("  Removing", n_zero, "all-zero taxa\n")
  otu_counts <- otu_counts[row_sums > 0, ]
}

out_otu <- file.path(out_dir, "otu_table_clean.tsv")
write.table(otu_counts, out_otu, sep = "\t", row.names = FALSE, quote = FALSE)

cat("  Samples:", ncol(otu_counts) - 1, "\n")
cat("  Taxa:   ", nrow(otu_counts), "\n")
cat("  Saved:  ", out_otu, "\n\n")

# =============================================================================
# 2. TAXONOMY TABLE
# =============================================================================

cat("=== Processing taxonomy table ===\n")

tax_raw <- read.table(
  tax_path,
  sep              = "\t",
  header           = TRUE,
  stringsAsFactors = FALSE,
  check.names      = FALSE
)

# Identify the ASV ID column and taxon string column
id_col  <- grep("^(Feature|feature|#OTU|id)\\s*(ID|id)?$",
                colnames(tax_raw), ignore.case = TRUE, value = TRUE)[1]
tax_col <- grep("^Taxon", colnames(tax_raw), ignore.case = TRUE, value = TRUE)[1]

if (is.na(id_col))  stop("Could not find Feature ID column in taxonomy file.")
if (is.na(tax_col)) stop("Could not find Taxon column in taxonomy file.")

cat("  ID column:  ", id_col, "\n")
cat("  Tax column: ", tax_col, "\n")

# Standardise domain prefix
tax_raw[[tax_col]] <- gsub("^d__", "k__", tax_raw[[tax_col]])

# Split into rank columns
tax_split <- tax_raw %>%
  select(all_of(c(id_col, tax_col))) %>%
  separate(
    col   = all_of(tax_col),
    into  = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep   = ";\\s*",
    fill  = "right",
    extra = "drop"
  )

# Remove rank prefixes (k__, p__, c__, etc.)
ranks    <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
prefixes <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
for (i in seq_along(ranks)) {
  if (ranks[i] %in% colnames(tax_split)) {
    tax_split[[ranks[i]]] <- trimws(gsub(paste0("^", prefixes[i]), "", tax_split[[ranks[i]]]))
    # Replace NA and empty with "Unclassified"
    tax_split[[ranks[i]]][is.na(tax_split[[ranks[i]]]) | tax_split[[ranks[i]]] == ""] <- "Unclassified"
  }
}

# Rename ID column to standard name
colnames(tax_split)[1] <- "Feature ID"

out_tax <- file.path(out_dir, "taxonomy_ranks.tsv")
write.table(tax_split, out_tax, sep = "\t", row.names = FALSE, quote = FALSE)

cat("  Taxa:  ", nrow(tax_split), "\n")
cat("  Ranks: ", paste(ranks, collapse = ", "), "\n")
cat("  Saved: ", out_tax, "\n\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("=== Done ===\n")
cat("Files ready for MicrobialExplorer upload:\n")
cat("  OTU table:   ", out_otu, "\n")
cat("  Taxonomy:    ", out_tax, "\n")
cat("  Tree:         phylotree_mafft_rooted.nwk (no reformatting needed)\n")
cat("  Metadata:     prepare your own — see template on the upload page\n\n")
