# =============================================================================
# MICROBIOME PIPELINE — COMPREHENSIVE TEST SUITE
# =============================================================================
# Description : Sequential test runner for all 10 analysis modules.
#               Creates a simulated dataset, runs every function in every
#               module, and reports PASS/FAIL with timing and diagnostics.
# Usage       : source("test_all_modules.R")
#               Or from terminal: Rscript test_all_modules.R
# Output      : Console report + test_results/test_report.html
# Author      : Patricia
# =============================================================================

# --- 0. SETUP AND CONFIGURATION ---------------------------------------------

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║          MICROBIOME PIPELINE — TEST SUITE                   ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Set options
options(warn = 1)
options(error = NULL)
# Add these two:
suppressWarnings(suppressMessages(library(tidytree)))
options(lifecycle_verbosity = "quiet")

# Create output directory
OUTPUT_DIR   <- "test_results"
CACHE_DIR    <- file.path(OUTPUT_DIR, "cache")
PLOTS_DIR    <- file.path(OUTPUT_DIR, "plots")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_DIR,  recursive = TRUE, showWarnings = FALSE)

# --- Load required packages -------------------------------------------------
required_pkgs <- c(
  "phyloseq", "ggplot2", "dplyr", "tidyr", "vegan", "patchwork",
  "scales", "tibble", "RColorBrewer", "stringr", "forcats",
  "randomForest", "caret", "pROC", "e1071", "glmnet",
  "Hmisc", "ggrepel", "zoo", "igraph", "ggraph", "tidygraph"
)

optional_pkgs <- c(
  "lme4", "lmerTest", "DESeq2", "ANCOMBC", "ALDEx2",
  "SpiecEasi", "picante", "broom.mixed"
)

load_packages <- function(pkgs, optional = FALSE) {
  results <- sapply(pkgs, function(pkg) {
    tryCatch({
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
      TRUE
    }, error = function(e) FALSE)
  })
  missing <- names(results[!results])
  if (length(missing) > 0) {
    if (optional) {
      cat("  Optional packages not installed (skipping those tests):\n")
      cat("  ", paste(missing, collapse = ", "), "\n\n")
    } else {
      stop("Required packages missing: ", paste(missing, collapse = ", "),
           "\nInstall with: install.packages(c('",
           paste(missing, collapse = "','"), "'))")
    }
  }
  results
}

cat("Loading required packages...\n")
pkg_status   <- load_packages(required_pkgs, optional = FALSE)
opt_status   <- load_packages(optional_pkgs, optional = TRUE)

AVAILABLE_PKGS <- c(names(pkg_status[pkg_status]),
                    names(opt_status[opt_status]))

# Redirect tidytree/phyloseq class conflict messages to null
local({
  old <- getOption("warn")
  options(warn = -1)
  try(methods::removeClass("phylo", where = ".GlobalEnv"), silent = TRUE)
  try(methods::removeClass("phylo", where = "package:tidytree"), silent = TRUE)
  options(warn = old)
})

# Suppress phylo class conflict noise from tidytree
assignInNamespace("phylo", getClass("phylo", where = "phyloseq"), "phyloseq")

# --- Source all modules -----------------------------------------------------
MODULE_FILES <- list(
  qc            = "R/microbiome_qc.R",
  composition   = "R/microbiome_composition.R",
  alpha         = "R/microbiome_alpha_diversity.R",
  beta          = "R/microbiome_beta_diversity.R",
  differential  = "R/microbiome_differential_abundance.R",
  functional    = "R/microbiome_functional.R",
  network       = "R/microbiome_network.R",
  ml            = "R/microbiome_ml_classification.R",
  correlation   = "R/microbiome_correlation.R",
  longitudinal  = "R/microbiome_longitudinal.R"
)

cat("Sourcing module files...\n")
source_status <- sapply(names(MODULE_FILES), function(mod) {
  path <- MODULE_FILES[[mod]]
  if (!file.exists(path)) {
    cat("  [MISSING]", path, "\n")
    return(FALSE)
  }
  tryCatch({
    source(path)
    cat("  [OK]     ", path, "\n")
    TRUE
  }, error = function(e) {
    cat("  [ERROR]  ", path, ":", conditionMessage(e), "\n")
    FALSE
  })
})

cat("\n")

# =============================================================================
# TEST INFRASTRUCTURE
# =============================================================================

# Global results tracker
TEST_RESULTS <- list()
TEST_START   <- Sys.time()

#' Record a test result
record_test <- function(module, test_name, passed, time_secs,
                         message = "", warning_msg = "") {
  key <- paste0(module, "::", test_name)
  TEST_RESULTS[[key]] <<- list(
    module      = module,
    test        = test_name,
    passed      = passed,
    time        = round(time_secs, 2),
    message     = message,
    warning     = warning_msg,
    timestamp   = Sys.time()
  )
}

#' Run a test with timing and error handling
run_test <- function(module, test_name, expr,
                      expected_class = NULL,
                      expected_names = NULL,
                      min_rows       = NULL,
                      check_fn       = NULL) {

  t_start  <- proc.time()["elapsed"]
  warn_msg <- ""
  result   <- NULL
  passed   <- FALSE
  message  <- ""

  tryCatch({
    result <- withCallingHandlers(
      eval(expr),
      warning = function(w) {
        warn_msg <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    )

    # Class check
    if (!is.null(expected_class)) {
      if (!inherits(result, expected_class) &&
          !any(class(result) %in% expected_class)) {
        stop("Expected class '", paste(expected_class, collapse="/"),
             "' got '", paste(class(result), collapse="/"), "'")
      }
    }

    # Names check
    if (!is.null(expected_names)) {
      missing_names <- setdiff(expected_names, names(result))
      if (length(missing_names) > 0) {
        stop("Missing expected names: ", paste(missing_names, collapse = ", "))
      }
    }

    # Min rows check
    if (!is.null(min_rows)) {
      actual_rows <- if (is.data.frame(result)) nrow(result)
                     else if (is.list(result) && "data" %in% names(result))
                       nrow(result$data)
                     else NA
      if (!is.na(actual_rows) && actual_rows < min_rows) {
        stop("Result has ", actual_rows, " rows; expected at least ", min_rows)
      }
    }

    # Custom check function
    if (!is.null(check_fn)) {
      check_result <- check_fn(result)
      if (!isTRUE(check_result)) {
        stop("Custom check failed: ", check_result)
      }
    }

    passed  <- TRUE
    message <- if (nchar(warn_msg) > 0) paste("WARN:", warn_msg) else "OK"

  }, error = function(e) {
    passed  <<- FALSE
    message <<- conditionMessage(e)
  })

  t_end    <- proc.time()["elapsed"]
  t_secs   <- t_end - t_start

  status   <- if (passed) "\033[32m[PASS]\033[0m" else "\033[31m[FAIL]\033[0m"
  time_str <- if (t_secs > 60) paste0(round(t_secs/60, 1), "m")
               else paste0(round(t_secs, 1), "s")

  cat(sprintf("  %-6s %-45s %6s\n", status, test_name, time_str))
  if (!passed) cat("         └─", message, "\n")
  if (nchar(warn_msg) > 0 && passed) cat("         └─ Warning:", warn_msg, "\n")

  record_test(module, test_name, passed, t_secs, message, warn_msg)
  invisible(result)
}

# Section header printer
print_section <- function(module_number, module_name) {
  cat("\n")
  cat(paste0(strrep("─", 64), "\n"))
  cat(sprintf(" Module %02d: %s\n", module_number, module_name))
  cat(paste0(strrep("─", 64), "\n"))
}


# =============================================================================
# SIMULATE TEST DATASET
# =============================================================================

print_section(0, "Simulating Test Dataset")

run_test("setup", "Create simulated phyloseq object", {
  
  source("R/simulate_microbiome_data.R")
  
  sim <- simulate_microbiome_data(
    output_dir   = file.path(CACHE_DIR, "simulated"),
    dataset_type = "network_enriched",
    seed         = 42,
    n_samples    = 60,
    n_taxa       = 200,
    n_subjects   = 20,
    n_timepoints = 3,
    add_tree     = TRUE
  )
  
  SIM_FILES <<- sim$files
  PS_TEST   <<- sim$ps
  
  PS_TEST
  
}, expected_class = "phyloseq",
check_fn = function(ps) {
  if (nsamples(ps) < 10) return("Too few samples")
  if (ntaxa(ps) < 10)    return("Too few taxa")
  if (is.null(phy_tree(ps, errorIfNULL = FALSE))) return("Tree missing")
  TRUE
})

cat(sprintf("  Dataset: %d samples x %d taxa\n", nsamples(PS_TEST), ntaxa(PS_TEST)))
cat(sprintf("  Groups: %s\n", paste(unique(sample_data(PS_TEST)$disease_status),
                                    collapse = ", ")))


# =============================================================================
# MODULE 1 — QC
# =============================================================================

print_section(1, "Quality Control (microbiome_qc.R)")

if (source_status["qc"]) {

  M1_PS <- run_test("qc", "qc_sequencing_depth()", {
    qc_sequencing_depth(PS_TEST, min_reads = 100, group_var = "disease_status")
  }, expected_class = "list",
     expected_names = c("plot", "summary", "flagged", "data"))

  run_test("qc", "qc_rarefaction_curves()", {
    qc_rarefaction_curves(PS_TEST, step = 200, group_var = "disease_status",
                           n_samples = 20)
  }, expected_class = "list",
     expected_names = c("plot", "data"))

  M1_FILTERED <- run_test("qc", "qc_filter_taxa()", {
    qc_filter_taxa(PS_TEST, min_prevalence = 0.10, min_abundance = 5)
  }, expected_class = "list",
     expected_names = c("ps_filtered", "plot", "summary"),
     check_fn = function(r) {
       if (ntaxa(r$ps_filtered) == 0) return("No taxa after filtering")
       if (ntaxa(r$ps_filtered) >= ntaxa(PS_TEST))
         return("Filter had no effect")
       TRUE
     })

  M1_SAMPLE <- run_test("qc", "qc_filter_samples()", {
    qc_filter_samples(PS_TEST, min_reads = 50, min_taxa = 5)
  }, expected_class = "list",
     expected_names = c("ps_filtered", "removed_samples"))

  # Save filtered object for downstream modules
  PS_CLEAN <- M1_FILTERED$ps_filtered
  saveRDS(PS_CLEAN, file.path(CACHE_DIR, "ps_clean.rds"))
  cat(sprintf("  Cached: ps_clean.rds (%d samples, %d taxa)\n",
              nsamples(PS_CLEAN), ntaxa(PS_CLEAN)))

} else {
  cat("  [SKIP] Module file not found\n")
  PS_CLEAN <- PS_TEST
}


# =============================================================================
# MODULE 2 — COMPOSITION
# =============================================================================

print_section(2, "Taxonomic Composition (microbiome_composition.R)")

if (source_status["composition"]) {

  M2_AGG <- run_test("composition", "agglomerate_taxa() — Phylum", {
    agglomerate_taxa(PS_CLEAN, rank = "Phylum", transform = "relative", top_n = 8)
  }, expected_class = "phyloseq")

  run_test("composition", "agglomerate_taxa() — Genus CLR", {
    agglomerate_taxa(PS_CLEAN, rank = "Genus", transform = "clr", top_n = 20)
  }, expected_class = "phyloseq")

  run_test("composition", "plot_composition_bars()", {
    ps_phy <- agglomerate_taxa(PS_CLEAN, rank = "Phylum",
                                transform = "relative", top_n = 8)
    plot_composition_bars(ps_phy, rank = "Phylum",
                           group_var = "disease_status")
  })

  run_test("composition", "plot_mean_composition()", {
    ps_phy <- agglomerate_taxa(PS_CLEAN, rank = "Phylum",
                                transform = "relative", top_n = 8)
    plot_mean_composition(ps_phy, rank = "Phylum",
                           group_var = "disease_status", top_n = 6)
  })

  run_test("composition", "plot_abundance_heatmap()", {
    plot_abundance_heatmap(PS_CLEAN, rank = "Genus", top_n = 20,
                            group_var = "disease_status", transform = "clr")
  })

  run_test("composition", "make_abundance_table()", {
    make_abundance_table(PS_CLEAN, rank = "Genus",
                          group_var = "disease_status")
  }, expected_class = "list",
     expected_names = c("overall"))

  run_test("composition", "identify_core_microbiome()", {
    identify_core_microbiome(PS_CLEAN, prevalence_cuts = c(0.5, 0.75),
                              min_abundance = 0.001)
  }, expected_class = "list",
     expected_names = c("plot", "core_taxa", "prevalence"))

  run_test("composition", "calculate_fb_ratio()", {
    calculate_fb_ratio(PS_CLEAN, group_var = "disease_status")
  }, expected_class = "list",
     expected_names = c("plot", "data"))

} else {
  cat("  [SKIP] Module file not found\n")
}


# =============================================================================
# MODULE 3 — ALPHA DIVERSITY
# =============================================================================

print_section(3, "Alpha Diversity (microbiome_alpha_diversity.R)")

if (source_status["alpha"]) {

  M3_DIV <- run_test("alpha", "calculate_alpha_diversity()", {
    calculate_alpha_diversity(PS_CLEAN, rarefaction = TRUE,
                               rare_depth = 50, n_rare = 3)
  }, expected_class = "data.frame",
     min_rows = 5,
     check_fn = function(df) {
       required_cols <- c("sample", "observed", "shannon", "simpson", "pielou")
       missing <- setdiff(required_cols, colnames(df))
       if (length(missing) > 0) return(paste("Missing cols:", paste(missing, collapse=",")))
       if (any(df$shannon < 0, na.rm = TRUE)) return("Negative Shannon values")
       if (any(df$simpson < 0 | df$simpson > 1, na.rm = TRUE))
         return("Simpson out of [0,1] range")
       TRUE
     })

  run_test("alpha", "plot_alpha_diversity()", {
    plot_alpha_diversity(M3_DIV,
                         metrics   = c("observed", "shannon", "pielou"),
                         group_var = "disease_status",
                         test      = "wilcox")
  }, expected_class = "list",
     expected_names = c("plot", "stats"))

  run_test("alpha", "plot_diversity_gradient()", {
    plot_diversity_gradient(M3_DIV, metric = "shannon",
                             x_var = "age", group_var = "disease_status")
  })

  run_test("alpha", "plot_diversity_correlations()", {
    plot_diversity_correlations(M3_DIV, method = "spearman")
  })

  run_test("alpha", "summarise_alpha_stats()", {
    summarise_alpha_stats(M3_DIV, group_var = "disease_status", test = "wilcox")
  }, expected_class = "list",
     expected_names = c("descriptive", "tests"))

  saveRDS(M3_DIV, file.path(CACHE_DIR, "alpha_diversity.rds"))

} else {
  cat("  [SKIP] Module file not found\n")
}


# =============================================================================
# MODULE 4 — BETA DIVERSITY
# =============================================================================

print_section(4, "Beta Diversity (microbiome_beta_diversity.R)")

if (source_status["beta"]) {

  M4_DISTS <- run_test("beta", "compute_distances() — bray + jaccard", {
    compute_distances(PS_CLEAN,
                       methods     = c("bray", "jaccard"),
                       rarefaction = TRUE,
                       rare_depth  = 50)
  }, expected_class = "list",
     check_fn = function(r) {
       if (!"bray" %in% names(r)) return("bray distance missing")
       if (!inherits(r$bray, "dist")) return("bray is not a dist object")
       TRUE
     })

  M4_PCOA <- run_test("beta", "run_ordination() — PCoA", {
    run_ordination(M4_DISTS$bray, method = "PCoA", k = 3)
  }, expected_class = "list",
     expected_names = c("coords", "var_exp", "method"))

  run_test("beta", "run_ordination() — NMDS", {
    run_ordination(M4_DISTS$bray, method = "NMDS")
  }, expected_class = "list",
     expected_names = c("coords", "stress", "method"))

  run_test("beta", "plot_ordination()", {
    plot_ordination(M4_PCOA, PS_CLEAN,
                     group_var = "disease_status", ellipse = TRUE)
  })

  run_test("beta", "run_permanova()", {
    run_permanova(M4_DISTS$bray, PS_CLEAN,
                   formula_rhs  = "disease_status",
                   permutations = 99)
  }, expected_class = "list",
     expected_names = c("permanova", "betadisper"))

  run_test("beta", "plot_betadisper()", {
    plot_betadisper(M4_DISTS$bray, PS_CLEAN,
                     group_var = "disease_status")
  })

  run_test("beta", "plot_distance_heatmap()", {
    plot_distance_heatmap(M4_DISTS$bray, PS_CLEAN,
                           group_var = "disease_status")
  })

  run_test("beta", "run_mantel_test()", {
    run_mantel_test(M4_DISTS$bray, M4_DISTS$jaccard,
                     permutations = 99)
  }, expected_class = "list",
     expected_names = c("result", "plot"))

  saveRDS(list(distances = M4_DISTS, pcoa = M4_PCOA),
          file.path(CACHE_DIR, "beta_diversity.rds"))

} else {
  cat("  [SKIP] Module file not found\n")
}


# =============================================================================
# MODULE 5 — DIFFERENTIAL ABUNDANCE
# =============================================================================

print_section(5, "Differential Abundance (microbiome_differential_abundance.R)")

if (source_status["differential"]) {

  M5_PREP <- run_test("differential", "prepare_da_data()", {
    prepare_da_data(PS_CLEAN, rank = "Genus",
                    min_prevalence = 0.10, min_count = 5)
  }, expected_class = "phyloseq",
     check_fn = function(ps) {
       if (ntaxa(ps) == 0) return("No taxa after DA preparation")
       TRUE
     })

  M5_RF_DA <- run_test("differential", "run_ancombc() [skip if not installed]", {
    if (!"ANCOMBC" %in% AVAILABLE_PKGS) {
      data.frame(taxon = "placeholder", comparison = "test",
                 lfc = 0, q_value = 1, diff_abund = FALSE, method = "skipped")
    } else {
      run_ancombc(M5_PREP, formula = "disease_status",
                   group_var = "disease_status", alpha = 0.05)
    }
  }, expected_class = "data.frame")

  run_test("differential", "run_deseq2() [skip if not installed]", {
    if (!"DESeq2" %in% AVAILABLE_PKGS) {
      cat("     [INFO] DESeq2 not installed — skipping\n")
      NULL
    } else {
      run_deseq2(M5_PREP, group_var = "disease_status",
                  reference = "Healthy", alpha = 0.05, lfc_threshold = 0.5)
    }
  })

  # Test visualisation functions with synthetic results
  M5_SYNTH <- data.frame(
    taxon      = paste0("Genus_", seq_len(50)),
    comparison = "IBD_vs_Healthy",
    lfc        = rnorm(50, sd = 1.5),
    se         = abs(rnorm(50, mean = 0.3, sd = 0.1)),
    q_value    = c(runif(15, 0, 0.04), runif(35, 0.05, 1)),
    diff_abund = c(rep(TRUE, 15), rep(FALSE, 35)),
    method     = "synthetic",
    stringsAsFactors = FALSE
  )

  run_test("differential", "plot_volcano()", {
    plot_volcano(M5_SYNTH, alpha = 0.05, lfc_threshold = 0.5,
                  top_n_label = 10)
  })

  run_test("differential", "plot_effect_sizes()", {
    plot_effect_sizes(M5_SYNTH, alpha = 0.05, lfc_threshold = 0.5,
                       top_n = 15, se_col = "se")
  })

  run_test("differential", "plot_da_heatmap()", {
    sig_taxa <- M5_SYNTH$taxon[M5_SYNTH$diff_abund]
    plot_da_heatmap(PS_CLEAN, sig_taxa, group_var = "disease_status",
                     rank = "Genus", transform = "zscore")
  })

  run_test("differential", "build_consensus()", {
    build_consensus(
      result_list = list(method_a = M5_SYNTH,
                          method_b = M5_SYNTH %>%
                            mutate(diff_abund = sample(c(TRUE,FALSE), 50,
                                                        replace=TRUE,
                                                        prob=c(0.3,0.7)),
                                   method = "method_b")),
      alpha = 0.05, min_methods = 2
    )
  }, expected_class = "list",
     expected_names = c("consensus", "all_results"))

} else {
  cat("  [SKIP] Module file not found\n")
}


# =============================================================================
# MODULE 6 — FUNCTIONAL PREDICTION
# =============================================================================

print_section(6, "Functional Prediction (microbiome_functional.R)")

if (source_status["functional"]) {

  # Simulate PICRUSt2-like output
  N_PATHWAYS <- 200
  N_KO       <- 500

  M6_PATH_MAT <- matrix(
    rnbinom(N_PATHWAYS * nsamples(PS_CLEAN), mu = 100, size = 0.5),
    nrow = N_PATHWAYS,
    dimnames = list(
      paste0("PWY_", seq_len(N_PATHWAYS)),
      sample_names(PS_CLEAN)
    )
  )

  M6_KO_MAT <- matrix(
    rnbinom(N_KO * nsamples(PS_CLEAN), mu = 50, size = 0.5),
    nrow = N_KO,
    dimnames = list(
      paste0("K", str_pad(seq_len(N_KO), 5, pad = "0")),
      sample_names(PS_CLEAN)
    )
  )

  M6_META <- as.data.frame(sample_data(PS_CLEAN))

  run_test("functional", "plot_pathway_abundance()", {
    plot_pathway_abundance(M6_PATH_MAT, M6_META,
                            group_var = "disease_status", top_n = 20)
  }, expected_class = "list",
     expected_names = c("plot", "summary", "top_pathways"))

  run_test("functional", "analyse_functional_diversity()", {
    analyse_functional_diversity(M6_PATH_MAT, M6_META,
                                  group_var = "disease_status",
                                  feature_type = "Pathway")
  }, expected_class = "list",
     expected_names = c("diversity", "p_alpha", "p_beta"))

  run_test("functional", "test_functional_da()", {
    test_functional_da(M6_PATH_MAT, M6_META,
                        group_var = "disease_status",
                        feature_type = "Pathway", alpha = 0.05)
  }, expected_class = "list",
     expected_names = c("results", "p_bubble", "n_sig"))

  run_test("functional", "summarise_kegg_categories()", {
    summarise_kegg_categories(M6_KO_MAT, M6_META,
                               group_var = "disease_status")
  }, expected_class = "list",
     expected_names = c("plot", "summary", "matrix"))

  run_test("functional", "plot_functional_heatmap()", {
    plot_functional_heatmap(M6_PATH_MAT, M6_META,
                             group_var = "disease_status",
                             top_n = 30, feature_type = "Pathway")
  })

  run_test("functional", "plot_nsti_quality()", {
    simulated_p2 <- list(metadata = M6_META)
    plot_nsti_quality(simulated_p2, group_var = "disease_status")
  }, expected_class = "list",
     expected_names = c("plot", "data", "median_nsti"))

} else {
  cat("  [SKIP] Module file not found\n")
}


# =============================================================================
# MODULE 7 — NETWORK ANALYSIS
# =============================================================================

print_section(7, "Network Analysis (microbiome_network.R)")

if (source_status["network"]) {

  M7_PREP <- run_test("network", "prepare_network_data()", {
    prepare_network_data(PS_CLEAN, rank = "Genus",
                          min_prevalence = 0.30, max_taxa = 40)
  }, expected_class = "phyloseq")

  M7_NET <- run_test("network", "build_network() - spearman", {
    build_network(M7_PREP, method = "spearman",
                  cor_threshold = 0.3, p_threshold = 0.05)
  }, expected_class = "list",
  expected_names = c("graph", "cor_mat", "adj_mat"),
  check_fn = function(r) {
    if (!inherits(r$graph, "igraph")) return("graph is not igraph")
    if (vcount(r$graph) == 0) return("Empty network (0 nodes)")
    TRUE
  })

  M7_TOPO <- run_test("network", "calculate_network_topology()", {
    calculate_network_topology(M7_NET$graph, ps = M7_PREP, rank = "Genus")
  }, expected_class = "list",
  expected_names = c("global", "nodes", "community"),
  check_fn = function(r) {
    if (is.null(r$nodes)) return("nodes table is NULL")
    if (!is.data.frame(r$nodes)) return("nodes is not a data.frame")
    if (NROW(r$nodes) == 0) return("No node metrics computed")
    if (!"hub_score" %in% names(r$nodes)) return("hub_score missing")
    TRUE
  })

  run_test("network", "plot_network()", {
    plot_network(M7_NET$graph, M7_TOPO$nodes,
                  layout = "fr", colour_by = "module",
                  size_by = "degree", label_hubs = TRUE)
  })

  run_test("network", "analyse_hub_taxa()", {
    analyse_hub_taxa(M7_NET$graph, M7_TOPO$nodes, top_n = 10)
  }, expected_class = "list",
     expected_names = c("hub_df", "plot"))

  run_test("network", "analyse_modules()", {
    analyse_modules(M7_NET$graph, M7_TOPO$nodes,
                     ps = M7_PREP, rank = "Phylum")
  }, expected_class = "list",
     expected_names = c("module_stats", "plot_comp"))

  run_test("network", "compare_random_network()", {
    compare_random_network(M7_NET$graph, n_random = 20)
  }, expected_class = "list",
     expected_names = c("sw_index", "plot"),
     check_fn = function(r) {
       if (is.na(r$sw_index)) return("sw_index is NA")
       TRUE
     })

} else {
  cat("  [SKIP] Module file not found\n")
}


# =============================================================================
# MODULE 8 — MACHINE LEARNING
# =============================================================================

print_section(8, "Machine Learning (microbiome_ml_classification.R)")

if (source_status["ml"]) {

  M8_DATA <- run_test("ml", "prepare_ml_data()", {
    prepare_ml_data(PS_CLEAN, group_var = "disease_status",
                    rank = "Genus", transform = "clr",
                    test_fraction = 0.25, max_features = 50)
  }, expected_class = "list",
     expected_names = c("train", "test", "features", "labels"),
     check_fn = function(r) {
       if (nrow(r$train) == 0) return("Empty training set")
       if (!"label" %in% colnames(r$train)) return("label column missing")
       TRUE
     })

  M8_RF <- run_test("ml", "train_random_forest()", {
    train_random_forest(M8_DATA, n_trees = 50,
                         cv_folds = 3, cv_repeats = 1, tune_mtry = FALSE)
  }, expected_class = "list",
     expected_names = c("model", "predictions", "probabilities", "confusion"),
     check_fn = function(r) {
       acc <- r$confusion$overall["Accuracy"]
       if (is.na(acc)) return("Accuracy is NA")
       TRUE
     })

  M8_LASSO <- run_test("ml", "train_lasso()", {
    train_lasso(M8_DATA, cv_folds = 5, alpha = 1)
  }, expected_class = "list",
     expected_names = c("model", "selected_features", "n_selected"))

  run_test("ml", "plot_feature_importance()", {
    plot_feature_importance(M8_RF, M8_LASSO, top_n = 20, ps = PS_CLEAN)
  }, expected_class = "list",
     expected_names = c("plot", "plots", "importance"))

  run_test("ml", "plot_roc_curves()", {
    plot_roc_curves(
      model_results = list(RF = M8_RF, LASSO = M8_LASSO),
      ml_data       = M8_DATA
    )
  }, expected_class = "list",
     expected_names = c("roc_objects", "auc_table", "plot"),
     check_fn = function(r) {
       if (nrow(r$auc_table) == 0) return("Empty AUC table")
       if (any(r$auc_table$auc < 0 | r$auc_table$auc > 1))
         return("AUC out of [0,1] range")
       TRUE
     })

  run_test("ml", "plot_confusion_matrix()", {
    plot_confusion_matrix(M8_RF$confusion, model_name = "RandomForest_Test")
  }, expected_class = "list",
     expected_names = c("plot", "p_cm"))

  run_test("ml", "plot_learning_curve()", {
    plot_learning_curve(M8_DATA, fractions = c(0.3, 0.6, 1.0),
                         cv_folds = 3, n_trees = 50)
  }, expected_class = "list",
     expected_names = c("plot", "data"))

  run_test("ml", "assess_feature_stability()", {
    assess_feature_stability(M8_DATA, n_boot = 10, top_n = 15, n_trees = 50)
  }, expected_class = "list",
     expected_names = c("stability", "plot"))

} else {
  cat("  [SKIP] Module file not found\n")
}


# =============================================================================
# MODULE 9 — CORRELATION AND ASSOCIATION
# =============================================================================

print_section(9, "Correlation and Association (microbiome_correlation.R)")

if (source_status["correlation"]) {

  M9_CORR <- run_test("correlation", "correlate_taxa_metadata()", {
    correlate_taxa_metadata(
      PS_CLEAN,
      rank       = "Genus",
      meta_vars  = c("age", "bmi", "crp", "calprotectin"),
      transform  = "clr",
      method     = "spearman",
      alpha      = 0.05
    )
  }, expected_class = "data.frame",
     min_rows = 1,
     check_fn = function(df) {
       required_cols <- c("taxon", "metadata", "rho", "q_value", "significant")
       missing <- setdiff(required_cols, colnames(df))
       if (length(missing) > 0)
         return(paste("Missing columns:", paste(missing, collapse = ",")))
       if (any(abs(df$rho) > 1, na.rm = TRUE))
         return("Correlation |rho| > 1")
       TRUE
     })

  run_test("correlation", "plot_taxa_metadata_heatmap()", {
    plot_taxa_metadata_heatmap(M9_CORR, top_n = 20, alpha = 0.05)
  })

  run_test("correlation", "plot_association_bubbles()", {
    plot_association_bubbles(M9_CORR, top_n = 40, alpha = 0.05)
  })

  run_test("correlation", "plot_top_associations()", {
    plot_top_associations(PS_CLEAN, M9_CORR, rank = "Genus",
                           top_n = 6, group_var = "disease_status")
  })

  M9_TT <- run_test("correlation", "compute_taxa_taxa_correlation()", {
    compute_taxa_taxa_correlation(PS_CLEAN, rank = "Genus",
                                   top_n = 25, transform = "clr",
                                   alpha = 0.05)
  }, expected_class = "list",
     expected_names = c("cor_matrix", "p_matrix", "plot", "n_sig"),
     check_fn = function(r) {
       if (!is.matrix(r$cor_matrix)) return("cor_matrix is not a matrix")
       if (any(abs(r$cor_matrix) > 1 + 1e-6, na.rm = TRUE))
         return("Correlation > 1 in matrix")
       TRUE
     })

  run_test("correlation", "plot_metadata_correlations()", {
    plot_metadata_correlations(PS_CLEAN,
                                meta_vars = c("age", "bmi", "crp",
                                              "calprotectin"),
                                method = "spearman")
  }, expected_class = "list",
     expected_names = c("cor_matrix", "p_matrix", "plot"))

  # Cross-domain test with synthetic external data
  N_EXT <- 30
  EXT_MAT <- matrix(
    rnorm(N_EXT * nsamples(PS_CLEAN)),
    nrow     = N_EXT,
    dimnames = list(paste0("Metabolite_", seq_len(N_EXT)),
                    sample_names(PS_CLEAN))
  )

  run_test("correlation", "correlate_cross_domain()", {
    correlate_cross_domain(PS_CLEAN, external_mat = EXT_MAT,
                            rank = "Genus", alpha = 0.05,
                            top_n_taxa = 15, top_n_ext = 15)
  }, expected_class = "list",
     expected_names = c("results", "plot", "n_sig"))

} else {
  cat("  [SKIP] Module file not found\n")
}


# =============================================================================
# MODULE 10 — LONGITUDINAL ANALYSIS
# =============================================================================

print_section(10, "Longitudinal Analysis (microbiome_longitudinal.R)")

if (source_status["longitudinal"]) {

  M10_PREP <- run_test("longitudinal", "prepare_longitudinal_data()", {
    prepare_longitudinal_data(
      PS_CLEAN,
      time_var    = "timepoint",
      subject_var = "subject_id",
      group_var   = "disease_status",
      rank        = "Genus",
      min_timepoints = 2
    )
  }, expected_class = "list",
     expected_names = c("ps", "ps_raw", "metadata", "diversity",
                         "features", "time_points", "n_subjects"),
     check_fn = function(r) {
       if (r$n_subjects == 0) return("No subjects found")
       if (length(r$time_points) < 2) return("Fewer than 2 time points")
       if (nrow(r$diversity) == 0) return("Empty diversity data frame")
       TRUE
     })

  run_test("longitudinal", "plot_diversity_trajectories()", {
    plot_diversity_trajectories(M10_PREP,
                                 metrics = c("shannon", "observed"),
                                 test    = TRUE)
  }, expected_class = "list",
     expected_names = c("plot", "stats"))

  run_test("longitudinal", "plot_beta_trajectories()", {
    plot_beta_trajectories(M10_PREP, distance = "bray",
                            rarefaction = TRUE, rare_depth = 50)
  }, expected_class = "list",
     expected_names = c("plot", "pcoa_coords", "dist_baseline"))

  run_test("longitudinal", "analyse_stability()", {
    analyse_stability(M10_PREP, distance = "bray")
  }, expected_class = "list",
     expected_names = c("plot", "stability_df", "subject_stability"),
     check_fn = function(r) {
       if (nrow(r$stability_df) == 0) return("Empty stability data frame")
       if (any(r$stability_df$distance < 0, na.rm = TRUE))
         return("Negative distances")
       TRUE
     })

  run_test("longitudinal", "plot_composition_over_time()", {
    plot_composition_over_time(M10_PREP, rank = "Phylum",
                                top_n = 6, facet_by = "group")
  }, expected_class = "list",
     expected_names = c("plot", "p_area", "data"))

  run_test("longitudinal", "run_lme_over_time() [skip if lme4 missing]", {
    if (!all(c("lme4", "lmerTest") %in% AVAILABLE_PKGS)) {
      cat("     [INFO] lme4/lmerTest not installed — skipping\n")
      list(results = data.frame(), plot = ggplot(), n_sig = 0)
    } else {
      run_lme_over_time(M10_PREP, interaction = TRUE, alpha = 0.05, top_n = 10)
    }
  }, expected_class = "list")

  run_test("longitudinal", "analyse_intervention_response()", {
    analyse_intervention_response(
      M10_PREP,
      pre_timepoints     = c(0, 4),
      post_timepoints    = c(8, 12, 16),
      intervention_label = "Test Intervention"
    )
  }, expected_class = "list",
     expected_names = c("plot", "p_paired", "p_response", "responders"))

  run_test("longitudinal", "detect_changepoints()", {
    detect_changepoints(M10_PREP, metric = "shannon", window = 2)
  }, expected_class = "list",
     expected_names = c("plot", "data", "n_changepoints"))

} else {
  cat("  [SKIP] Module file not found\n")
}


# =============================================================================
# GENERATE TEST REPORT
# =============================================================================

print_section(99, "Test Report")

# Compile results
results_df <- bind_rows(lapply(TEST_RESULTS, function(r) {
  data.frame(
    Module   = r$module,
    Test     = r$test,
    Status   = ifelse(r$passed, "PASS", "FAIL"),
    Time_s   = r$time,
    Message  = r$message,
    stringsAsFactors = FALSE
  )
}))

# Summary statistics
n_total   <- nrow(results_df)
n_pass    <- sum(results_df$Status == "PASS")
n_fail    <- sum(results_df$Status == "FAIL")
n_skip    <- sum(results_df$Status == "SKIP")
pct_pass  <- round(100 * n_pass / max(n_total, 1), 1)
total_time <- round(as.numeric(difftime(Sys.time(), TEST_START, units = "secs")), 1)

# --- Console summary --------------------------------------------------------
cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat(sprintf("║  RESULTS: %d/%d passed (%.1f%%)  |  Total time: %.1fs      \n",
            n_pass, n_total, pct_pass, total_time))
cat("╠══════════════════════════════════════════════════════════════╣\n")

# Per-module summary
module_summary <- results_df %>%
  group_by(Module) %>%
  summarise(
    Total    = n(),
    Passed   = sum(Status == "PASS"),
    Failed   = sum(Status == "FAIL"),
    Time_s   = round(sum(Time_s), 1),
    .groups  = "drop"
  ) %>%
  mutate(Status = ifelse(Failed == 0, "✓ PASS", "✗ FAIL"))

for (i in seq_len(nrow(module_summary))) {
  row    <- module_summary[i, ]
  icon   <- if (row$Status == "✓ PASS") "\033[32m✓\033[0m" else "\033[31m✗\033[0m"
  cat(sprintf("║  %s %-18s  %d/%d passed  (%5.1fs)              \n",
              icon, row$Module, row$Passed, row$Total, row$Time_s))
}

cat("╠══════════════════════════════════════════════════════════════╣\n")

# Failed tests detail
if (n_fail > 0) {
  cat("║  FAILED TESTS:\n")
  failed_tests <- results_df %>% filter(Status == "FAIL")
  for (i in seq_len(nrow(failed_tests))) {
    cat(sprintf("║  ✗ [%s] %s\n",
                failed_tests$Module[i], failed_tests$Test[i]))
    cat(sprintf("║      └─ %s\n",
                str_trunc(failed_tests$Message[i], 55)))
  }
  cat("╠══════════════════════════════════════════════════════════════╣\n")
}

cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# --- Save results to CSV ---------------------------------------------------
write.csv(results_df,
          file.path(OUTPUT_DIR, "test_results.csv"),
          row.names = FALSE)

# --- Generate HTML report ---------------------------------------------------
html_rows <- apply(results_df, 1, function(row) {
  status_colour <- if (row["Status"] == "PASS") "#27ae60" else
                   if (row["Status"] == "FAIL") "#e74c3c" else "#95a5a6"
  sprintf(
    "<tr>
       <td>%s</td>
       <td>%s</td>
       <td style='color:%s;font-weight:bold'>%s</td>
       <td>%.2fs</td>
       <td style='font-size:0.85em;color:#7f8c8d'>%s</td>
     </tr>",
    row["Module"], row["Test"], status_colour, row["Status"],
    as.numeric(row["Time_s"]),
    str_trunc(row["Message"], 80)
  )
})

html_content <- sprintf(
'<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>Microbiome Pipeline — Test Report</title>
  <style>
    body  { font-family: -apple-system, sans-serif; max-width: 1100px;
            margin: 2rem auto; padding: 0 1rem; background: #f8f9fa; }
    h1    { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 0.5rem; }
    .summary { display: flex; gap: 1rem; margin: 1.5rem 0; flex-wrap: wrap; }
    .card { background: white; border-radius: 8px; padding: 1rem 1.5rem;
            box-shadow: 0 2px 8px rgba(0,0,0,0.08); flex: 1; min-width: 120px; }
    .card-val { font-size: 2rem; font-weight: 700; }
    .card-lab { font-size: 0.8rem; color: #95a5a6; text-transform: uppercase; }
    .pass { color: #27ae60; }
    .fail { color: #e74c3c; }
    .pct  { color: #3498db; }
    table { width: 100%%; border-collapse: collapse; background: white;
            border-radius: 8px; overflow: hidden;
            box-shadow: 0 2px 8px rgba(0,0,0,0.08); }
    th    { background: #2c3e50; color: white; padding: 0.75rem 1rem;
            text-align: left; font-size: 0.85rem; }
    td    { padding: 0.6rem 1rem; border-bottom: 1px solid #ecf0f1;
            font-size: 0.9rem; }
    tr:hover { background: #f8f9fa; }
    .footer { margin-top: 2rem; color: #95a5a6; font-size: 0.8rem; }
  </style>
</head>
<body>
  <h1>Microbiome Pipeline — Test Report</h1>
  <p>Generated: %s | Total runtime: %.1fs</p>

  <div class="summary">
    <div class="card">
      <div class="card-val">%d</div>
      <div class="card-lab">Total tests</div>
    </div>
    <div class="card">
      <div class="card-val pass">%d</div>
      <div class="card-lab">Passed</div>
    </div>
    <div class="card">
      <div class="card-val fail">%d</div>
      <div class="card-lab">Failed</div>
    </div>
    <div class="card">
      <div class="card-val pct">%.1f%%</div>
      <div class="card-lab">Pass rate</div>
    </div>
  </div>

  <table>
    <thead>
      <tr>
        <th>Module</th><th>Test</th><th>Status</th>
        <th>Time</th><th>Message</th>
      </tr>
    </thead>
    <tbody>
      %s
    </tbody>
  </table>

  <div class="footer">
    R version: %s | Platform: %s
  </div>
</body>
</html>',
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  total_time,
  n_total, n_pass, n_fail, pct_pass,
  paste(html_rows, collapse = "\n      "),
  R.version$version.string,
  R.version$platform
)

report_path <- file.path(OUTPUT_DIR, "test_report.html")
writeLines(html_content, report_path)

cat("  Test results saved to:\n")
cat("   ", file.path(OUTPUT_DIR, "test_results.csv"), "\n")
cat("   ", report_path, "(open in browser)\n\n")

# --- Session info -----------------------------------------------------------
capture.output(sessionInfo(),
               file = file.path(OUTPUT_DIR, "session_info.txt"))
cat("  Session info saved to:", file.path(OUTPUT_DIR, "session_info.txt"), "\n\n")

# --- Final verdict ----------------------------------------------------------
if (n_fail == 0) {
  cat("\033[32m✓ ALL TESTS PASSED — Ready to integrate with Shiny\033[0m\n\n")
} else {
  cat("\033[31m✗", n_fail, "test(s) failed — Fix these before Shiny integration\033[0m\n")
  cat("  See test_results/test_report.html for details\n\n")
}

# Return results invisibly for programmatic use
invisible(list(
  results    = results_df,
  n_pass     = n_pass,
  n_fail     = n_fail,
  pct_pass   = pct_pass,
  total_time = total_time,
  all_passed = n_fail == 0
))
