# =============================================================================
# validators.R — Input validation for all ViromeAnalyst data types
# =============================================================================
# Source this file before any user-facing data processing:
#   source("R/validators.R")
#
# All functions return list(valid, errors, warnings) so callers can
# decide whether to stop, warn, or proceed.
# =============================================================================

# source("R/utils.R")

# =============================================================================
# OTU / ASV TABLE VALIDATION
# =============================================================================

#' Validate an OTU or ASV count table
#'
#' @param otu_mat A matrix or data frame (taxa × samples)
#' @return list(valid, errors, warnings)
validate_otu_table <- function(otu_mat) {
  errors   <- c()
  warnings <- c()

  if (!is.matrix(otu_mat) && !is.data.frame(otu_mat))
    errors <- c(errors, "OTU table must be a matrix or data frame.")

  if (is.data.frame(otu_mat)) otu_mat <- as.matrix(otu_mat)

  if (any(otu_mat < 0, na.rm = TRUE))
    errors <- c(errors, "OTU table contains negative values. Counts must be >= 0.")

  if (all(otu_mat == 0, na.rm = TRUE))
    errors <- c(errors, "OTU table is entirely zeros.")

  if (nrow(otu_mat) < 5)
    errors <- c(errors, paste0("Too few taxa: ", nrow(otu_mat), " (minimum 5)."))

  if (ncol(otu_mat) < 5)
    errors <- c(errors, paste0("Too few samples: ", ncol(otu_mat), " (minimum 5)."))

  if (any(is.na(otu_mat)))
    errors <- c(errors, "OTU table contains NA values. Please replace with 0.")

  # Warnings
  sparsity <- mean(otu_mat == 0)
  if (sparsity > 0.95)
    warnings <- c(warnings, sprintf("Very high sparsity: %.1f%% zeros.", sparsity * 100))

  if (!all(otu_mat == round(otu_mat), na.rm = TRUE))
    warnings <- c(warnings, "OTU table contains non-integer values. DA methods expect counts.")

  n_zero_samples <- sum(colSums(otu_mat) == 0)
  if (n_zero_samples > 0)
    warnings <- c(warnings, paste0(n_zero_samples, " sample(s) have zero total reads."))

  list(valid = length(errors) == 0, errors = errors, warnings = warnings)
}


# =============================================================================
# TAXONOMY TABLE VALIDATION
# =============================================================================

#' Validate a taxonomy table
#'
#' @param tax_mat A matrix or data frame (taxa × ranks)
#' @return list(valid, errors, warnings)
validate_taxonomy_table <- function(tax_mat) {
  errors   <- c()
  warnings <- c()
  
  if (!is.matrix(tax_mat) && !is.data.frame(tax_mat)) {
    errors <- c(errors, "Taxonomy table must be a matrix or data frame.")
    return(list(valid = FALSE, errors = errors, warnings = warnings))
  }
  
  tax_mat <- as.data.frame(tax_mat, stringsAsFactors = FALSE)
  
  # Clean column names: remove spaces and hidden BOM characters
  colnames(tax_mat) <- trimws(colnames(tax_mat))
  colnames(tax_mat) <- gsub("^\ufeff", "", colnames(tax_mat))
  
  # If first column is an ID column from write.csv(row.names = TRUE), remove it
  first_col <- colnames(tax_mat)[1]
  possible_id_cols <- c("", "X", "Taxa", "OTU", "ASV", "FeatureID", "Feature.ID", "Feature ID", "featureid", "Unnamed: 0")
  


  if (first_col %in% possible_id_cols) {
    rownames(tax_mat) <- tax_mat[[1]]
    tax_mat <- tax_mat[, -1, drop = FALSE]
  }
  
  # Clean again after removing first column
  colnames(tax_mat) <- trimws(colnames(tax_mat))
  colnames(tax_mat) <- gsub("^\ufeff", "", colnames(tax_mat))
  
  available <- colnames(tax_mat)
  missing_ranks <- setdiff(c("Phylum", "Genus"), available)
  
  if (length(missing_ranks) > 0) {
    errors <- c(
      errors,
      paste(
        "Missing required rank columns:",
        paste(missing_ranks, collapse = ", "),
        "| Available columns:",
        paste(available, collapse = ", ")
      )
    )
  }
  
  pct_na_genus <- if ("Genus" %in% available) {
    round(100 * mean(is.na(tax_mat[, "Genus"]) | tax_mat[, "Genus"] == ""), 1)
  } else {
    NA
  }
  
  if (!is.na(pct_na_genus) && pct_na_genus > 50) {
    warnings <- c(
      warnings,
      paste0(pct_na_genus, "% of taxa have NA/empty values at Genus level.")
    )
  }
  
  list(valid = length(errors) == 0, errors = errors, warnings = warnings)
}


# =============================================================================
# METADATA VALIDATION
# =============================================================================

#' Validate a sample metadata table
#'
#' @param meta_df       A data frame (samples × variables)
#' @param required_vars Character vector of required column names
#' @return list(valid, errors, warnings)
validate_metadata <- function(meta_df, required_vars = NULL) {
  errors   <- c()
  warnings <- c()

  if (!is.data.frame(meta_df))
    errors <- c(errors, "Metadata must be a data frame.")

  if (nrow(meta_df) < 5)
    errors <- c(errors, paste0("Too few samples in metadata: ",
                                nrow(meta_df), " (minimum 5)."))

  if (!is.null(required_vars)) {
    missing <- setdiff(required_vars, colnames(meta_df))
    if (length(missing) > 0)
      errors <- c(errors, paste("Missing required metadata columns:",
                                 paste(missing, collapse = ", ")))
  }

  # Check for small groups in categorical variables
  for (col in colnames(meta_df)) {
    if (is.character(meta_df[[col]]) || is.factor(meta_df[[col]])) {
      group_counts <- table(meta_df[[col]])
      small_groups <- names(group_counts[group_counts < 3])
      if (length(small_groups) > 0)
        warnings <- c(warnings, sprintf(
          "Variable '%s' has groups with < 3 samples: %s",
          col, paste(small_groups, collapse = ", ")
        ))
    }
  }

  # Check for all-NA columns
  all_na_cols <- names(meta_df)[sapply(meta_df, function(x) all(is.na(x)))]
  if (length(all_na_cols) > 0)
    warnings <- c(warnings, paste("All-NA columns:", paste(all_na_cols, collapse = ", ")))

  list(valid = length(errors) == 0, errors = errors, warnings = warnings)
}


# =============================================================================
# PHYLOSEQ OBJECT VALIDATION
# =============================================================================

#' Validate a phyloseq object before analysis
#'
#' @param ps              A phyloseq object
#' @param required_vars   Required metadata column names
#' @param min_samples     Minimum sample count. Default = 5
#' @param min_taxa        Minimum taxon count. Default = 10
#' @return list(valid, errors, warnings)
validate_phyloseq <- function(ps,
                               required_vars = NULL,
                               min_samples   = 5,
                               min_taxa      = 10) {
  errors   <- c()
  warnings <- c()

  if (!inherits(ps, "phyloseq"))
    return(list(valid = FALSE,
                errors = "Input is not a phyloseq object.",
                warnings = character(0)))

  if (nsamples(ps) < min_samples)
    errors <- c(errors, paste0("Too few samples: ", nsamples(ps),
                                " (minimum ", min_samples, ")."))

  if (ntaxa(ps) < min_taxa)
    errors <- c(errors, paste0("Too few taxa: ", ntaxa(ps),
                                " (minimum ", min_taxa, ")."))

  if (min(sample_sums(ps)) == 0)
    errors <- c(errors, "One or more samples have zero total reads.")

  if (any(is.na(as.matrix(otu_table(ps)))))
    errors <- c(errors, "OTU table contains NA values.")

  if (any(as.matrix(otu_table(ps)) < 0))
    errors <- c(errors, "OTU table contains negative values.")

  # Required metadata
  if (!is.null(required_vars)) {
    available_vars <- sample_variables(ps)
    missing_vars   <- setdiff(required_vars, available_vars)
    if (length(missing_vars) > 0)
      errors <- c(errors, paste("Missing required metadata columns:",
                                 paste(missing_vars, collapse = ", ")))
  }

  # Warnings
  if (min(sample_sums(ps)) < 1000)
    warnings <- c(warnings, paste0("Low read count in some samples: ",
                                    min(sample_sums(ps)), " reads."))

  expected_ranks <- c("Phylum", "Genus")
  missing_ranks  <- setdiff(expected_ranks, rank_names(ps))
  if (length(missing_ranks) > 0)
    warnings <- c(warnings, paste("Missing recommended taxonomic ranks:",
                                   paste(missing_ranks, collapse = ", ")))

  # Check taxa_are_rows
  if (!taxa_are_rows(ps))
    warnings <- c(warnings, "taxa_are_rows is FALSE — modules assume TRUE. Consider transposing.")

  list(valid = length(errors) == 0, errors = errors, warnings = warnings)
}


# =============================================================================
# LONGITUDINAL DATA VALIDATION
# =============================================================================

#' Validate metadata for longitudinal analysis
#'
#' @param meta_df       A data frame of sample metadata
#' @param time_var      Column name for time variable
#' @param subject_var   Column name for subject identifier
#' @param min_subjects  Minimum number of subjects. Default = 5
#' @param min_timepoints Minimum time points per subject. Default = 2
#' @return list(valid, errors, warnings)
validate_longitudinal <- function(meta_df,
                                   time_var       = "timepoint",
                                   subject_var    = "subject_id",
                                   min_subjects   = 5,
                                   min_timepoints = 2) {
  errors   <- c()
  warnings <- c()

  # Check required columns
  for (col in c(time_var, subject_var)) {
    if (!col %in% colnames(meta_df))
      errors <- c(errors, paste0("Required column '", col, "' not found in metadata."))
  }

  if (length(errors) > 0)
    return(list(valid = FALSE, errors = errors, warnings = warnings))

  # Time must be numeric
  time_vals <- suppressWarnings(as.numeric(as.character(meta_df[[time_var]])))
  if (any(is.na(time_vals)))
    errors <- c(errors, paste0("Column '", time_var, "' contains non-numeric values."))

  # Subject count
  n_subjects <- n_distinct(meta_df[[subject_var]])
  if (n_subjects < min_subjects)
    errors <- c(errors, paste0("Too few subjects: ", n_subjects,
                                " (minimum ", min_subjects, ")."))

  # Time points per subject
  tp_per_subj <- meta_df %>%
    group_by(.data[[subject_var]]) %>%
    summarise(n_tp = n(), .groups = "drop")

  low_tp <- sum(tp_per_subj$n_tp < min_timepoints)
  if (low_tp > 0)
    warnings <- c(warnings, paste0(low_tp, " subject(s) have fewer than ",
                                    min_timepoints, " time points."))

  list(valid = length(errors) == 0, errors = errors, warnings = warnings)
}


# =============================================================================
# PICRUST2 OUTPUT VALIDATION
# =============================================================================

#' Validate a PICRUSt2 output directory
#'
#' @param picrust_dir Path to PICRUSt2 output directory
#' @return list(valid, errors, warnings, found_files)
validate_picrust2_dir <- function(picrust_dir) {
  errors   <- c()
  warnings <- c()

  if (!dir.exists(picrust_dir))
    return(list(valid = FALSE,
                errors = paste0("PICRUSt2 directory not found: ", picrust_dir),
                warnings = character(0), found_files = character(0)))

  # Expected output files
  expected_files <- list(
    pathway = c("pathways_out/path_abun_unstrat.tsv.gz",
                "pathways_out/path_abun_unstrat.tsv"),
    ko      = c("KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",
                "KO_metagenome_out/pred_metagenome_unstrat.tsv"),
    ec      = c("EC_metagenome_out/pred_metagenome_unstrat.tsv.gz",
                "EC_metagenome_out/pred_metagenome_unstrat.tsv")
  )

  found_files <- c()
  for (type in names(expected_files)) {
    found <- FALSE
    for (f in expected_files[[type]]) {
      if (file.exists(file.path(picrust_dir, f))) {
        found       <- TRUE
        found_files <- c(found_files, setNames(f, type))
        break
      }
    }
    if (!found)
      warnings <- c(warnings, paste0("PICRUSt2 ", type, " output file not found."))
  }

  if (length(found_files) == 0)
    errors <- c(errors, "No PICRUSt2 output files found in directory.")

  list(valid = length(errors) == 0, errors = errors,
       warnings = warnings, found_files = found_files)
}


# =============================================================================
# CONVENIENCE WRAPPER — VALIDATE EVERYTHING
# =============================================================================

#' Run all validations and print a summary report
#'
#' @param otu_mat       OTU table matrix
#' @param tax_mat       Taxonomy table matrix
#' @param meta_df       Metadata data frame
#' @param required_vars Required metadata columns
#' @return TRUE if all valid, stops with error if not
validate_all_inputs <- function(otu_mat, tax_mat, meta_df,
                                 required_vars = NULL) {

  cat("=== Input validation ===\n")
  all_valid <- TRUE

  # OTU table
  cat("  OTU table... ")
  v_otu <- validate_otu_table(otu_mat)
  if (v_otu$valid) {
    cat("✓\n")
  } else {
    cat("✗\n")
    for (e in v_otu$errors) cat("    ERROR:", e, "\n")
    all_valid <- FALSE
  }
  for (w in v_otu$warnings) cat("    WARN:", w, "\n")

  # Taxonomy
  cat("  Taxonomy table... ")
  v_tax <- validate_taxonomy_table(tax_mat)
  if (v_tax$valid) {
    cat("✓\n")
  } else {
    cat("✗\n")
    for (e in v_tax$errors) cat("    ERROR:", e, "\n")
    all_valid <- FALSE
  }
  for (w in v_tax$warnings) cat("    WARN:", w, "\n")

  # Metadata
  cat("  Metadata... ")
  v_meta <- validate_metadata(meta_df, required_vars = required_vars)
  if (v_meta$valid) {
    cat("✓\n")
  } else {
    cat("✗\n")
    for (e in v_meta$errors) cat("    ERROR:", e, "\n")
    all_valid <- FALSE
  }
  for (w in v_meta$warnings) cat("    WARN:", w, "\n")

  # Sample alignment
  cat("  Sample alignment... ")
  otu_samples  <- colnames(otu_mat)
  meta_samples <- rownames(meta_df)
  n_shared     <- length(intersect(otu_samples, meta_samples))

  if (n_shared == 0) {
    cat("✗\n")
    cat("    ERROR: No shared samples between OTU table and metadata.\n")
    cat("    OTU table sample names:    ", paste(head(otu_samples, 3), collapse=", "), "...\n")
    cat("    Metadata sample names:     ", paste(head(meta_samples, 3), collapse=", "), "...\n")
    all_valid <- FALSE
  } else {
    cat("✓ (", n_shared, "shared samples)\n")
    if (n_shared < length(otu_samples))
      cat("    WARN:", length(otu_samples) - n_shared,
          "OTU samples not in metadata.\n")
    if (n_shared < length(meta_samples))
      cat("    WARN:", length(meta_samples) - n_shared,
          "metadata samples not in OTU table.\n")
  }

  cat("\n")

  if (!all_valid) {
    stop("Input validation failed. Please fix the errors listed above before proceeding.")
  }

  cat("  All inputs valid. Proceeding with import.\n\n")
  invisible(TRUE)
}

cat("[validators.R] Input validation functions loaded.\n")
