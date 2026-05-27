# =============================================================================
# precompute.R — Cache management and pre-computation for ViromeAnalyst
# =============================================================================
# Run precompute_all() ONCE when a user uploads data.
# All subsequent analyses load from cache (< 1 second).
#
# Usage:
#   source("R/precompute.R")
#   precompute_all(ps, group_var = "disease_status",
#                  output_dir = "cache/")
#   cache <- load_cache("cache/")
# =============================================================================

# source("R/utils.R")

# =============================================================================
# PRE-COMPUTATION WRAPPER
# =============================================================================

#' Run all expensive computations once and save to cache.
#'
#' Call this function once when a user uploads their data.
#' Shiny then loads from cache on all subsequent user interactions.
#'
#' @param ps          A filtered phyloseq object (from Module 1 QC)
#' @param group_var   Primary grouping variable in metadata
#' @param rank        Primary taxonomic rank. Default = "Genus"
#' @param output_dir  Base cache directory. Default = "cache/"
#' @param seed        Random seed for reproducibility. Default = 42
#' @param methods     Distance methods to compute. Default = c("bray","jaccard")
#' @param rare_depth  Rarefaction depth. Default = min sample depth
#' @return Named list of all cached objects (also saved as RDS files)

precompute_all <- function(ps,
                            group_var   = NULL,
                            rank        = "Genus",
                            output_dir  = "cache/",
                            seed        = 42,
                            methods     = c("bray", "jaccard"),
                            rare_depth  = NULL) {

  t_start <- proc.time()["elapsed"]
  cat("\n=== Pre-computation started ===\n")
  cat("  Samples:", nsamples(ps), "\n")
  cat("  Taxa:   ", ntaxa(ps), "\n")
  cat("  Rank:   ", rank, "\n\n")

  # Create cache subdirectories
  subdirs <- c("qc", "composition", "alpha", "beta",
               "differential", "network", "ml", "correlation")
  for (d in subdirs) {
    dir.create(file.path(output_dir, d),
               recursive = TRUE, showWarnings = FALSE)
  }

  cache   <- list()
  results <- list()

  # --- 1. Agglomerated phyloseq objects ------------------------------------
  cache_step("Agglomerating to Phylum", {
    ps_phylum <- tax_glom(ps, taxrank = "Phylum", NArm = FALSE)
    taxa_names(ps_phylum) <- as.character(tax_table(ps_phylum)[, "Phylum"])
    saveRDS(ps_phylum, file.path(output_dir, "ps_phylum.rds"))
    results$ps_phylum <- ps_phylum
  })

  cache_step("Agglomerating to Genus", {
    ps_genus <- tax_glom(ps, taxrank = rank, NArm = FALSE)
    new_names <- as.character(tax_table(ps_genus)[, rank])
    new_names[is.na(new_names)] <- paste0("Unknown_", seq_len(sum(is.na(new_names))))
    dup_idx <- duplicated(new_names)
    new_names[dup_idx] <- paste0(new_names[dup_idx], "_dup", seq_len(sum(dup_idx)))
    taxa_names(ps_genus) <- new_names
    saveRDS(ps_genus, file.path(output_dir, paste0("ps_", tolower(rank), ".rds")))
    results$ps_genus <- ps_genus
  })

  # --- 2. Relative abundance transformed objects ---------------------------
  cache_step("Relative abundance transformation", {
    ps_rel <- transform_sample_counts(results$ps_genus, function(x) x / sum(x))
    saveRDS(ps_rel, file.path(output_dir, "ps_relative.rds"))
    results$ps_relative <- ps_rel
  })

  cache_step("CLR transformation", {
    otu_mat <- as.matrix(otu_table(results$ps_genus))
    if (!taxa_are_rows(results$ps_genus)) otu_mat <- t(otu_mat)
    clr_mat <- clr_transform(otu_mat)
    saveRDS(clr_mat, file.path(output_dir, "clr_matrix.rds"))
    results$clr_matrix <- clr_mat
  })

  # --- 3. Alpha diversity --------------------------------------------------
  cache_step("Alpha diversity (rarefaction)", {
    if (!pkg_available("vegan")) stop("vegan required for alpha diversity")
    library(vegan)

    if (is.null(rare_depth)) rare_depth <- min(sample_sums(ps))

    set.seed(seed)
    otu_t      <- t(as.matrix(otu_table(ps)))
    if (taxa_are_rows(ps)) otu_t <- t(otu_t)

    rare_mat   <- rrarefy(otu_t[rowSums(otu_t) >= rare_depth, ],
                           sample = rare_depth)

    diversity_df <- data.frame(
      sample     = rownames(rare_mat),
      observed   = rowSums(rare_mat > 0),
      shannon    = round(diversity(rare_mat, index = "shannon"), 4),
      simpson    = round(diversity(rare_mat, index = "simpson"), 4),
      pielou     = NA,
      stringsAsFactors = FALSE
    )
    diversity_df$pielou <- ifelse(
      diversity_df$observed > 1,
      round(diversity_df$shannon / log(diversity_df$observed), 4), 0
    )

    meta_df <- data.frame(sample_data(ps)) %>%
      rownames_to_column("sample")
    diversity_df <- left_join(diversity_df, meta_df, by = "sample")

    saveRDS(diversity_df, file.path(output_dir, "alpha", "diversity_df.rds"))
    results$diversity_df <- diversity_df
  })

  # --- 4. Beta diversity ---------------------------------------------------
  cache_step("Distance matrices", {
    if (!pkg_available("vegan")) stop("vegan required")
    library(vegan)

    otu_mat <- as.matrix(otu_table(ps))
    if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)

    set.seed(seed)
    if (!is.null(rare_depth)) {
      otu_rare <- rrarefy(otu_mat[rowSums(otu_mat) >= rare_depth, ],
                           sample = rare_depth)
    } else {
      otu_rare <- otu_mat
    }

    distances <- list(
      bray    = vegdist(otu_rare, method = "bray"),
      jaccard = vegdist(otu_rare, method = "jaccard", binary = TRUE)
    )

    # Aitchison distance (CLR-based Euclidean)
    clr_rare <- results$clr_matrix[rownames(otu_rare), , drop = FALSE]
    distances$aitchison <- dist(clr_rare)

    saveRDS(distances, file.path(output_dir, "beta", "distances.rds"))
    results$distances <- distances
  })

  cache_step("PCoA ordinations", {
    ordinations <- list()
    for (method in names(results$distances)) {
      d   <- results$distances[[method]]
      pco <- cmdscale(d, k = 3, eig = TRUE)
      eig_pos <- pco$eig[pco$eig > 0]
      var_exp <- round(eig_pos / sum(eig_pos) * 100, 1)
      ordinations[[method]] <- list(
        coords  = as.data.frame(pco$points) %>%
          setNames(paste0("PC", 1:3)) %>%
          rownames_to_column("sample"),
        var_exp = var_exp,
        method  = "PCoA",
        dist    = method
      )
    }
    saveRDS(ordinations, file.path(output_dir, "beta", "ordinations.rds"))
    results$ordinations <- ordinations
  })

  # --- 5. Taxonomic composition --------------------------------------------
  cache_step("Composition summary tables", {
    comp_summary <- list()
    for (r in c("Phylum", "Family", rank)) {
      ps_agg  <- tryCatch(tax_glom(ps, taxrank = r, NArm = FALSE),
                           error = function(e) NULL)
      if (is.null(ps_agg)) next
      ps_rel  <- transform_sample_counts(ps_agg, function(x) x / sum(x))
      taxa_names(ps_rel) <- as.character(tax_table(ps_rel)[, r])

      otu_mat <- as.matrix(otu_table(ps_rel))
      if (!taxa_are_rows(ps_rel)) otu_mat <- t(otu_mat)

      comp_summary[[r]] <- data.frame(
        taxon      = rownames(otu_mat),
        mean_abund = rowMeans(otu_mat),
        prevalence = rowSums(otu_mat > 0) / ncol(otu_mat),
        stringsAsFactors = FALSE
      ) %>% arrange(desc(mean_abund))
    }
    saveRDS(comp_summary, file.path(output_dir, "composition", "summaries.rds"))
    results$comp_summary <- comp_summary
  })

  # --- 6. Metadata summary -------------------------------------------------
  cache_step("Metadata summary", {
    meta_df   <- data.frame(sample_data(ps))
    num_vars  <- names(meta_df)[sapply(meta_df, is.numeric)]
    cat_vars  <- names(meta_df)[sapply(meta_df, function(x) is.character(x)|is.factor(x))]

    meta_summary <- list(
      n_samples    = nsamples(ps),
      num_vars     = num_vars,
      cat_vars     = cat_vars,
      group_var    = group_var,
      group_levels = if (!is.null(group_var) && group_var %in% cat_vars)
        unique(meta_df[[group_var]]) else NULL,
      sample_sums  = sample_sums(ps)
    )
    saveRDS(meta_summary, file.path(output_dir, "meta_summary.rds"))
    results$meta_summary <- meta_summary
  })

  # --- Done ----------------------------------------------------------------
  t_secs <- round(proc.time()["elapsed"] - t_start, 1)
  cat(sprintf("\n=== Pre-computation complete (%.1fs) ===\n", t_secs))
  cat("  Cache saved to:", output_dir, "\n")
  cat("  Files created:", length(list.files(output_dir, recursive = TRUE)), "\n\n")

  invisible(results)
}


# =============================================================================
# CACHE LOADER
# =============================================================================

#' Load all cached objects for Shiny use
#'
#' @param output_dir Cache directory. Default = "cache/"
#' @return Named list of all cached objects

load_cache <- function(output_dir = "cache/") {

  cat("=== Loading cache from:", output_dir, "===\n")

  safe_load <- function(path) {
    full_path <- file.path(output_dir, path)
    if (file.exists(full_path)) {
      tryCatch(readRDS(full_path),
               error = function(e) {
                 cat("  [WARN] Failed to load:", path, "\n")
                 NULL
               })
    } else {
      cat("  [MISS]", path, "\n")
      NULL
    }
  }

  cache <- list(
    ps_phylum      = safe_load("ps_phylum.rds"),
    ps_genus       = safe_load("ps_genus.rds"),
    ps_relative    = safe_load("ps_relative.rds"),
    clr_matrix     = safe_load("clr_matrix.rds"),
    diversity_df   = safe_load("alpha/diversity_df.rds"),
    distances      = safe_load("beta/distances.rds"),
    ordinations    = safe_load("beta/ordinations.rds"),
    comp_summary   = safe_load("composition/summaries.rds"),
    meta_summary   = safe_load("meta_summary.rds")
  )

  n_loaded <- sum(!sapply(cache, is.null))
  cat(sprintf("  Loaded %d/%d objects from cache.\n\n", n_loaded, length(cache)))

  invisible(cache)
}


# =============================================================================
# CACHE UTILITIES
# =============================================================================

#' Check if cache is fresh (newer than the source phyloseq)
#'
#' @param output_dir  Cache directory
#' @param ps_file     Path to the source phyloseq RDS file
#' @return TRUE if cache is fresh, FALSE if stale or missing
is_cache_fresh <- function(output_dir = "cache/", ps_file = NULL) {
  cache_files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
  if (length(cache_files) == 0) return(FALSE)
  if (is.null(ps_file))        return(TRUE)
  if (!file.exists(ps_file))   return(TRUE)

  ps_mtime     <- file.mtime(ps_file)
  cache_mtimes <- file.mtime(cache_files)
  all(cache_mtimes > ps_mtime)
}

#' Clear the cache directory
#'
#' @param output_dir Cache directory to clear
#' @param confirm    Whether to ask for confirmation. Default = TRUE
clear_cache <- function(output_dir = "cache/", confirm = TRUE) {
  if (confirm) {
    response <- readline(paste0("Clear cache in '", output_dir, "'? (y/n): "))
    if (tolower(response) != "y") {
      cat("Cache not cleared.\n")
      return(invisible(FALSE))
    }
  }
  files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
  file.remove(files)
  cat("Cache cleared:", length(files), "files removed.\n")
  invisible(TRUE)
}

#' Get cache status summary
#'
#' @param output_dir Cache directory
cache_status <- function(output_dir = "cache/") {
  files    <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
  n_files  <- length(files)
  size_mb  <- round(sum(file.size(files)) / 1e6, 1)
  cat(sprintf("Cache: %d files, %.1f MB in '%s'\n", n_files, size_mb, output_dir))
  invisible(list(n_files = n_files, size_mb = size_mb))
}


# =============================================================================
# INTERNAL HELPER
# =============================================================================

cache_step <- function(label, expr) {
  t0 <- proc.time()["elapsed"]
  cat(sprintf("  %-45s", paste0(label, "...")))
  tryCatch({
    eval(expr, envir = parent.frame())
    t1 <- proc.time()["elapsed"]
    cat(sprintf(" OK  (%.1fs)\n", t1 - t0))
  }, error = function(e) {
    cat(sprintf(" FAIL\n    %s\n", conditionMessage(e)))
  })
}

cat("[precompute.R] Cache management functions loaded.\n")
