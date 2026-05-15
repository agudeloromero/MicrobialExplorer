# =============================================================================
# MICROBIOME ANALYSIS PLATFORM
# SETUP AND TESTING PHASE — MASTER PLAN
# =============================================================================
# Version    : 1.0
# Author     : Patricia
# Phase      : Pre-Shiny Setup and Testing
# Goal       : Verified, tested, documented pipeline ready for Shiny
#              integration
# =============================================================================


# =============================================================================
# SECTION 1 — PROJECT FOLDER STRUCTURE
# =============================================================================

# Run this block ONCE to create the entire project structure

create_project_structure <- function(base_dir = "ViromeAnalyst") {

  dirs <- c(

    # --- Core application directories ---------------------------------------
    "R",                          # All module scripts live here
    "app",                        # Shiny app (to be built later)
    "app/ui",                     # Shiny UI components
    "app/server",                 # Shiny server components
    "app/modules",                # Shiny module files

    # --- Data directories ---------------------------------------------------
    "data/raw",                   # User-uploaded raw data (gitignored)
    "data/example",               # Example datasets for demo
    "data/reference",             # Reference databases (gitignored)

    # --- Cache directories --------------------------------------------------
    "cache",                      # Pre-computed objects (gitignored)
    "cache/qc",
    "cache/composition",
    "cache/alpha",
    "cache/beta",
    "cache/differential",
    "cache/functional",
    "cache/network",
    "cache/ml",
    "cache/correlation",
    "cache/longitudinal",

    # --- Output directories -------------------------------------------------
    "output",                     # Analysis outputs (gitignored)
    "output/plots",
    "output/tables",
    "output/reports",

    # --- Testing directories ------------------------------------------------
    "tests",                      # Test scripts
    "tests/results",              # Test output
    "tests/results/cache",
    "tests/results/plots",

    # --- Documentation directories ------------------------------------------
    "docs",                       # Documentation
    "docs/module_contracts",      # Module interface specifications

    # --- Configuration ------------------------------------------------------
    "config",                     # Configuration files

    # --- Logs ---------------------------------------------------------------
    "logs"                        # Runtime logs
  )

  cat("Creating project structure in:", base_dir, "\n\n")

  for (d in dirs) {
    full_path <- file.path(base_dir, d)
    dir.create(full_path, recursive = TRUE, showWarnings = FALSE)
    cat("  Created:", full_path, "\n")
  }

  # --- Create placeholder files -------------------------------------------

  # .gitignore
  gitignore_content <- c(
    "# R",
    ".Rhistory",
    ".RData",
    ".Rproj.user",
    "",
    "# Data (never commit raw data)",
    "data/raw/",
    "data/reference/",
    "",
    "# Cache (generated files)",
    "cache/",
    "",
    "# Output (generated files)",
    "output/",
    "",
    "# Logs",
    "logs/",
    "",
    "# renv (keep lock file, ignore library)",
    "renv/library/",
    ".Rprofile"
  )
  writeLines(gitignore_content, file.path(base_dir, ".gitignore"))

  # README skeleton
  readme_content <- c(
    "# ViromeAnalyst",
    "",
    "Interactive microbiome and virome analysis platform.",
    "",
    "## Structure",
    "- `R/`         — Analysis modules",
    "- `app/`       — Shiny application",
    "- `data/`      — Input data",
    "- `tests/`     — Test suite",
    "- `docs/`      — Documentation",
    "",
    "## Quick Start",
    "```r",
    "renv::restore()            # Restore package environment",
    "source('tests/test_all_modules.R')  # Run test suite",
    "shiny::runApp('app/')      # Launch app",
    "```"
  )
  writeLines(readme_content, file.path(base_dir, "README.md"))

  cat("\nProject structure created.\n")
  cat("Next step: copy your R module files into", file.path(base_dir, "R/"), "\n")
}

# Run it:
# create_project_structure("ViromeAnalyst")


# =============================================================================
# SECTION 2 — COMPLETE FILE LIST WITH LOCATIONS
# =============================================================================

# Every file in the project, what it does, and where it lives.

FILE_MANIFEST <- data.frame(
  File = c(

    # Core R modules
    "R/utils.R",
    "R/microbiome_qc.R",
    "R/microbiome_composition.R",
    "R/microbiome_alpha_diversity.R",
    "R/microbiome_beta_diversity.R",
    "R/microbiome_differential_abundance.R",
    "R/microbiome_functional.R",
    "R/microbiome_network.R",
    "R/microbiome_ml_classification.R",
    "R/microbiome_correlation.R",
    "R/microbiome_longitudinal.R",

    # Validation and contracts
    "R/validators.R",
    "R/precompute.R",

    # Tests
    "tests/test_all_modules.R",
    "tests/test_data_generation.R",
    "tests/test_validators.R",

    # Configuration
    "config/pipeline_config.R",
    "config/ui_config.R",

    # Documentation
    "docs/module_contracts/contracts.R",
    "docs/TROUBLESHOOTING.md",

    # renv
    "renv.lock",
    ".Rprofile"
  ),
  Status = c(
    "WRITE NOW", # utils.R
    "DONE",      # qc
    "DONE",      # composition
    "DONE",      # alpha
    "DONE",      # beta
    "DONE",      # differential
    "DONE",      # functional
    "DONE",      # network
    "DONE",      # ml
    "DONE",      # correlation
    "DONE",      # longitudinal
    "WRITE NOW", # validators
    "WRITE NOW", # precompute
    "DONE",      # test_all_modules
    "WRITE NOW", # test_data_generation
    "WRITE NOW", # test_validators
    "WRITE NOW", # pipeline_config
    "WRITE LATER", # ui_config (Shiny phase)
    "WRITE NOW",   # contracts
    "WRITE NOW",   # troubleshooting
    "AUTO",        # renv.lock
    "AUTO"         # .Rprofile
  ),
  Description = c(
    "Shared theme, colours, helper functions — sourced by all modules",
    "Quality control and filtering pipeline",
    "Taxonomic composition analysis",
    "Alpha diversity metrics and visualisation",
    "Beta diversity ordination and statistics",
    "Differential abundance — ANCOM-BC, DESeq2, LASSO consensus",
    "Functional prediction — PICRUSt2 integration",
    "Co-occurrence network analysis",
    "Machine learning classification",
    "Correlation and association analysis",
    "Longitudinal trajectory analysis",
    "Input validation functions for all data types",
    "Pre-computation caching — run once on upload",
    "Master test suite — all modules",
    "Test dataset generation functions",
    "Validator unit tests",
    "Analysis parameters and defaults",
    "Shiny UI configuration — colours, labels, tooltips",
    "Module interface contracts and documentation",
    "Troubleshooting guide for common errors",
    "Package version lockfile — generated by renv",
    "renv auto-loader — generated by renv"
  ),
  stringsAsFactors = FALSE
)

# Print the manifest
cat("\n=== File Manifest ===\n")
for (i in seq_len(nrow(FILE_MANIFEST))) {
  status_symbol <- switch(FILE_MANIFEST$Status[i],
    "DONE"        = "[✓]",
    "WRITE NOW"   = "[→]",
    "WRITE LATER" = "[…]",
    "AUTO"        = "[⚙]"
  )
  cat(sprintf("  %s %-45s %s\n",
              status_symbol,
              FILE_MANIFEST$File[i],
              FILE_MANIFEST$Description[i]))
}


# =============================================================================
# SECTION 3 — SEQUENTIAL TESTING PLAN
# =============================================================================

# Work through these steps IN ORDER.
# Do not move to the next step until the current one passes completely.

TESTING_PHASES <- list(

  phase_1 = list(
    name        = "Environment Setup",
    description = "Install packages, set up renv, verify R version",
    steps       = c(
      "1.1  Verify R >= 4.2.0 is installed",
      "1.2  Install renv: install.packages('renv')",
      "1.3  Initialise renv: renv::init()",
      "1.4  Install required packages (see Section 4)",
      "1.5  Install optional packages (see Section 4)",
      "1.6  Snapshot environment: renv::snapshot()",
      "1.7  Verify: renv::status() shows no issues"
    ),
    pass_criteria = "renv::status() reports 'No issues found'",
    estimated_time = "30–60 minutes"
  ),

  phase_2 = list(
    name        = "File and Structure Verification",
    description = "Confirm all files are in the right place",
    steps       = c(
      "2.1  Run create_project_structure('ViromeAnalyst')",
      "2.2  Copy all R module files into R/",
      "2.3  Copy test script into tests/",
      "2.4  Verify: all 11 R files exist in R/",
      "2.5  Create utils.R (extract shared code)",
      "2.6  Add source('R/utils.R') to top of each module",
      "2.7  Verify: source each module file individually"
    ),
    pass_criteria = "All 11 source() calls succeed with no errors",
    estimated_time = "1–2 hours"
  ),

  phase_3 = list(
    name        = "Simulated Data Test",
    description = "Run full test suite on simulated dataset",
    steps       = c(
      "3.1  Run: source('tests/test_all_modules.R')",
      "3.2  Open: tests/results/test_report.html",
      "3.3  Fix any FAIL results (see Section 6 — Common Errors)",
      "3.4  Re-run test suite until 100% pass rate",
      "3.5  Save session info: capture.output(sessionInfo())",
      "3.6  Review timing — flag any function > 10 seconds"
    ),
    pass_criteria = "All tests PASS. test_report.html shows 100%",
    estimated_time = "2–4 hours (including bug fixing)"
  ),

  phase_4 = list(
    name        = "Real Data Test — Small Dataset",
    description = "Test pipeline on a small public microbiome dataset",
    steps       = c(
      "4.1  Download curatedMetagenomicData or use HMP dataset",
      "4.2  Subset to 30–50 samples",
      "4.3  Run Module 1 (QC) — inspect output manually",
      "4.4  Run Module 2 (Composition) — check plots make sense",
      "4.5  Run Module 3 (Alpha) — verify diversity values are realistic",
      "4.6  Run Module 4 (Beta) — check PCoA separates groups",
      "4.7  Run Module 5 (DA) — confirm known DA taxa are detected",
      "4.8  Run remaining modules in order",
      "4.9  Compare results to published paper using same dataset"
    ),
    pass_criteria = "Results are biologically plausible and consistent with literature",
    estimated_time = "4–8 hours"
  ),

  phase_5 = list(
    name        = "Your Own Data Test",
    description = "Test pipeline on your actual research dataset",
    steps       = c(
      "5.1  Prepare your OTU/ASV table, taxonomy, and metadata files",
      "5.2  Run validate_otu_table(), validate_metadata()",
      "5.3  Run Module 1 (QC) — check flagged samples",
      "5.4  Decide on min_reads and min_prevalence thresholds",
      "5.5  Run all modules sequentially",
      "5.6  Save all plots to output/plots/",
      "5.7  Review biological plausibility with domain knowledge",
      "5.8  Document any dataset-specific parameter choices"
    ),
    pass_criteria = "Full pipeline completes on your data with no errors",
    estimated_time = "Full day"
  ),

  phase_6 = list(
    name        = "Performance Profiling",
    description = "Identify slow functions before Shiny",
    steps       = c(
      "6.1  Install profvis: install.packages('profvis')",
      "6.2  Profile each module on your real dataset",
      "6.3  Flag any function taking > 5 seconds",
      "6.4  Implement caching for slow functions (see precompute.R)",
      "6.5  Test with cache: confirm < 2 seconds for cached results",
      "6.6  Profile memory use: library(pryr); mem_used()",
      "6.7  Ensure memory < 4 GB for full pipeline run"
    ),
    pass_criteria = "All interactive functions run in < 5 seconds with cache",
    estimated_time = "2–4 hours"
  ),

  phase_7 = list(
    name        = "Documentation and Cleanup",
    description = "Final preparation before Shiny integration",
    steps       = c(
      "7.1  Write docs/TROUBLESHOOTING.md",
      "7.2  Write docs/module_contracts/contracts.R",
      "7.3  Add @param documentation to all exported functions",
      "7.4  Remove debug cat() statements from production code",
      "7.5  Standardise function argument names across all modules",
      "7.6  Final renv::snapshot()",
      "7.7  Commit everything to Git",
      "7.8  Tag release: git tag v0.1.0-pretesting"
    ),
    pass_criteria = "Clean Git commit with no uncommitted changes",
    estimated_time = "2–3 hours"
  )
)

cat("\n=== Sequential Testing Phases ===\n\n")
for (ph_name in names(TESTING_PHASES)) {
  ph <- TESTING_PHASES[[ph_name]]
  cat(sprintf("PHASE: %s — %s\n", ph$name, ph$description))
  cat(sprintf("  Time estimate: %s\n", ph$estimated_time))
  cat(sprintf("  Pass criteria: %s\n", ph$pass_criteria))
  cat("  Steps:\n")
  for (step in ph$steps) cat("   ", step, "\n")
  cat("\n")
}


# =============================================================================
# SECTION 4 — COMPLETE PACKAGE LIST
# =============================================================================

PACKAGES <- list(

  # Bioconductor — install first
  bioconductor_required = c(
    "phyloseq",        # Core microbiome data structure
    "BiocManager",     # Bioconductor package manager
    "microbiome",      # Additional microbiome utilities
    "decontam",        # Contamination detection
    "Hmisc"            # Correlation with p-values
  ),

  bioconductor_optional = c(
    "DESeq2",          # Differential abundance
    "ANCOMBC",         # Differential abundance (recommended)
    "ALDEx2",          # Differential abundance
    "picante",         # Faith's phylogenetic diversity
    "Tax4Fun2"         # Functional prediction (alternative to PICRUSt2)
  ),

  # CRAN — standard install
  cran_required = c(
    "ggplot2",         # Visualisation
    "dplyr",           # Data manipulation
    "tidyr",           # Data tidying
    "patchwork",       # Combining plots
    "vegan",           # Ecological statistics
    "scales",          # Axis formatting
    "tibble",          # Modern data frames
    "RColorBrewer",    # Colour palettes
    "stringr",         # String manipulation
    "forcats",         # Factor handling
    "ggrepel",         # Non-overlapping text labels
    "zoo",             # Rolling statistics
    "igraph",          # Network analysis
    "ggraph",          # Network visualisation
    "tidygraph",       # Tidy network manipulation
    "randomForest",    # Random forest classification
    "caret",           # ML training framework
    "pROC",            # ROC curves
    "e1071",           # SVM (caret dependency)
    "glmnet",          # LASSO regression
    "rstatix",         # Tidy statistical tests
    "ggpubr",          # Publication-ready plots
    "corrplot",        # Correlation visualisation
    "profvis",         # Performance profiling
    "renv"             # Package environment management
  ),

  cran_optional = c(
    "lme4",            # Linear mixed-effects models
    "lmerTest",        # P-values for lme4 models
    "broom.mixed",     # Tidy lme4 output
    "xgboost",         # Gradient boosting (ML module)
    "DALEX",           # ML explainability
    "pryr",            # Memory profiling
    "testthat",        # Formal unit testing
    "covr"             # Code coverage
  ),

  # GitHub — install with remotes
  github_optional = c(
    "zdk123/SpiecEasi" # Sparse inverse covariance for network analysis
  )
)

# Installation script
cat("\n=== Package Installation Script ===\n\n")
cat("# Step 1: Install BiocManager\n")
cat("install.packages('BiocManager')\n\n")

cat("# Step 2: Install Bioconductor required\n")
cat("BiocManager::install(c(\n")
cat(" ", paste0("  '", PACKAGES$bioconductor_required, "'",
               collapse = ",\n  "), "\n))\n\n")

cat("# Step 3: Install CRAN required\n")
cat("install.packages(c(\n")
cat(" ", paste0("  '", PACKAGES$cran_required, "'",
               collapse = ",\n  "), "\n))\n\n")

cat("# Step 4: Install optional (run if you want full functionality)\n")
cat("BiocManager::install(c(\n")
cat(" ", paste0("  '", PACKAGES$bioconductor_optional, "'",
               collapse = ",\n  "), "\n))\n\n")
cat("install.packages(c(\n")
cat(" ", paste0("  '", PACKAGES$cran_optional, "'",
               collapse = ",\n  "), "\n))\n\n")

cat("# Step 5: Lock the environment\n")
cat("renv::snapshot()\n\n")


# =============================================================================
# SECTION 5 — utils.R CONTENT TO CREATE NOW
# =============================================================================

# This file must exist before any module can be sourced.
# Create R/utils.R with this content:

UTILS_CONTENT <- '
# =============================================================================
# utils.R — Shared utilities for all microbiome analysis modules
# Source this file at the top of every module: source("R/utils.R")
# =============================================================================

# --- Shared ggplot2 theme ---------------------------------------------------
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

# --- Shared colour palettes -------------------------------------------------
TAXA_PALETTE <- c(
  "#3498db","#e74c3c","#2ecc71","#f39c12","#9b59b6",
  "#1abc9c","#e67e22","#34495e","#e91e63","#00bcd4",
  "#8bc34a","#ff5722","#607d8b","#795548","#ffc107",
  "#673ab7","#009688","#ff9800","#4caf50","#f44336",
  "#b0bec5"
)

KEGG_COLOURS <- c(
  "Metabolism"                     = "#27ae60",
  "Genetic Information Processing" = "#3498db",
  "Environmental Information Proc" = "#9b59b6",
  "Cellular Processes"             = "#e67e22",
  "Organismal Systems"             = "#e74c3c",
  "Human Diseases"                 = "#c0392b",
  "Unknown"                        = "#bdc3c7"
)

# --- Helper functions -------------------------------------------------------
pkg_available <- function(pkg) requireNamespace(pkg, quietly = TRUE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_run <- function(name, expr) {
  cat("[START]", name, "\n")
  tryCatch({
    result <- eval(expr)
    cat("[OK]   ", name, "\n")
    result
  }, error = function(e) {
    cat("[FAIL] ", name, ":", conditionMessage(e), "\n")
    NULL
  })
}

check_phyloseq <- function(ps, label = "phyloseq") {
  cat("\\n===", label, "===\\n")
  cat("  Samples:", nsamples(ps), "\\n")
  cat("  Taxa:   ", ntaxa(ps), "\\n")
  cat("  Ranks:  ", paste(rank_names(ps), collapse=", "), "\\n")
  meta_df  <- as.data.frame(sample_data(ps))
  num_vars <- names(meta_df)[sapply(meta_df, is.numeric)]
  cat_vars <- names(meta_df)[sapply(meta_df, function(x) is.character(x)|is.factor(x))]
  cat("  Numeric vars:", paste(num_vars, collapse=", "), "\\n")
  cat("  Group vars:  ", paste(cat_vars, collapse=", "), "\\n\\n")
}
'

writeLines(UTILS_CONTENT, "R/utils.R")
cat("utils.R created at R/utils.R\n\n")


# =============================================================================
# SECTION 6 — PIPELINE CONFIGURATION FILE
# =============================================================================

# Create config/pipeline_config.R — centralises all analysis defaults
# Change defaults here, not inside the module files

CONFIG_CONTENT <- '
# =============================================================================
# config/pipeline_config.R — Global analysis parameters
# =============================================================================

PIPELINE_CONFIG <- list(

  # QC parameters
  qc = list(
    min_reads      = 1000,
    min_taxa       = 10,
    min_prevalence = 0.05,
    min_abundance  = 10,
    rarefaction    = TRUE
  ),

  # Taxonomic rank for primary analysis
  primary_rank = "Genus",

  # Abundance transformation
  transform = "clr",

  # Beta diversity
  beta = list(
    methods      = c("bray", "jaccard"),
    permutations = 999,
    rare_depth   = NULL     # NULL = auto (min sample depth)
  ),

  # Differential abundance
  differential = list(
    alpha         = 0.05,
    lfc_threshold = 1,
    min_methods   = 2,
    p_adjust      = "BH"
  ),

  # Machine learning
  ml = list(
    test_fraction = 0.2,
    n_trees       = 500,
    cv_folds      = 5,
    cv_repeats    = 3,
    n_boot        = 50
  ),

  # Network
  network = list(
    method        = "spearman",
    cor_threshold = 0.6,
    p_threshold   = 0.05,
    min_prevalence = 0.20,
    max_taxa      = 100
  ),

  # Visualisation
  viz = list(
    top_n_taxa    = 20,
    top_n_features = 30,
    figure_width  = 12,
    figure_height = 8,
    dpi           = 300
  ),

  # Reproducibility
  seed = 42
)
'

writeLines(CONFIG_CONTENT, "config/pipeline_config.R")
cat("pipeline_config.R created at config/pipeline_config.R\n\n")


# =============================================================================
# SECTION 7 — COMMON ERRORS AND FIXES REFERENCE
# =============================================================================

KNOWN_ERRORS <- list(

  e1 = list(
    error   = "Error in tax_glom: NAs in taxonomic rank",
    cause   = "Taxonomy table has NA at the requested rank",
    fix     = "Add NArm = FALSE to tax_glom(), or pre-fill NAs:
               tax_table(ps)[is.na(tax_table(ps)[,'Genus']), 'Genus'] <- 'Unknown'",
    modules = c("composition", "differential", "network", "ml", "correlation")
  ),

  e2 = list(
    error   = "Error: must have at least 2 groups",
    cause   = "A sample filter removed all samples in one group",
    fix     = "Check: table(sample_data(ps)$group_var) after filtering.
               Reduce min_reads or min_prevalence thresholds.",
    modules = c("alpha", "beta", "differential", "ml")
  ),

  e3 = list(
    error   = "Error in adonis2: object of class 'dist' expected",
    cause   = "Sample order in distance matrix differs from metadata",
    fix     = "Ensure metadata rows match distance matrix labels:
               meta_df <- meta_df[labels(dist_obj), , drop = FALSE]",
    modules = c("beta", "longitudinal")
  ),

  e4 = list(
    error   = "NaN / Inf values after CLR transformation",
    cause   = "Zeros in OTU table — log(0) is undefined",
    fix     = "Always add pseudocount: apply(otu + 0.5, 2, function(x) log(x) - mean(log(x)))",
    modules = c("all modules using CLR")
  ),

  e5 = list(
    error   = "Error in lmer: number of observations <= number of random effects",
    cause   = "Too few samples per subject for the random effects structure",
    fix     = "Require at least 3 time points per subject (min_timepoints = 3),
               or simplify the random effects structure",
    modules = c("longitudinal", "correlation")
  ),

  e6 = list(
    error   = "Error in randomForest: Can't have empty classes",
    cause   = "One class has 0 samples in the training set after split",
    fix     = "Use balance_classes = TRUE, or increase training fraction,
               or check class sizes: table(sample_data(ps)$group_var)",
    modules = c("ml")
  ),

  e7 = list(
    error   = "Error in phyloseq_to_deseq2: rownames not matching",
    cause   = "Sample names in OTU table and metadata don't match",
    fix     = "Align: ps <- prune_samples(sample_names(ps), ps)
               after confirming rownames(sample_data(ps)) == sample_names(ps)",
    modules = c("differential")
  ),

  e8 = list(
    error   = "network has 0 edges after filtering",
    cause   = "Correlation threshold too high, or too few samples",
    fix     = "Lower cor_threshold (try 0.4), reduce p_threshold (try 0.10),
               or increase min_prevalence to include fewer but better-observed taxa",
    modules = c("network")
  ),

  e9 = list(
    error   = "Error: no shared samples between phyloseq and external matrix",
    cause   = "Sample names use different format in the two objects",
    fix     = "Check: intersect(sample_names(ps), colnames(external_mat))
               Clean names: colnames(external_mat) <- trimws(colnames(external_mat))",
    modules = c("correlation", "functional")
  ),

  e10 = list(
    error   = "package 'ANCOMBC' not found",
    cause   = "Bioconductor package not installed",
    fix     = "BiocManager::install('ANCOMBC') — run_ancombc() will fall back
               to ANCOM-BC1 if ANCOMBC2 fails. DESeq2 can be used as alternative.",
    modules = c("differential")
  )
)

cat("=== Known Errors Reference ===\n\n")
for (e in KNOWN_ERRORS) {
  cat(sprintf("ERROR: %s\n", e$error))
  cat(sprintf("  Cause:   %s\n", e$cause))
  cat(sprintf("  Modules: %s\n", paste(e$modules, collapse = ", ")))
  cat(sprintf("  Fix:     %s\n\n", str_wrap(e$fix, width = 60,
                                              exdent = 11)))
}


# =============================================================================
# SECTION 8 — DECISION POINTS BEFORE SHINY
# =============================================================================

# These are questions you must answer before starting Shiny development.
# Document your decisions in config/pipeline_config.R.

DECISIONS <- list(

  d1 = list(
    question = "What is your primary taxonomic rank?",
    options  = c("Genus (recommended)", "Species", "Family"),
    impact   = "Affects computation time and biological resolution",
    decision = "Genus"
  ),

  d2 = list(
    question = "Which DA method is primary?",
    options  = c("ANCOM-BC2 (best)", "DESeq2 (widely used)", "Consensus of both"),
    impact   = "Determines which results table the app shows by default",
    decision = "Consensus — show both, highlight agreement"
  ),

  d3 = list(
    question = "How will you handle the virome module in V2?",
    options  = c("Separate upload (vOTU table)", "Integrated with bacterial",
                  "Separate tab in same app"),
    impact   = "Shapes the Shiny UI navigation structure",
    decision = "Separate tab — teased in V1, active in V2"
  ),

  d4 = list(
    question = "What is the maximum dataset size the app should accept?",
    options  = c("500 samples (conservative)", "1000 samples", "No limit"),
    impact   = "Determines server specs and caching strategy",
    decision = "500 samples for V1 — revisit after profiling"
  ),

  d5 = list(
    question = "Do users need to download figures?",
    options  = c("PDF", "PNG", "Both"),
    impact   = "Affects Shiny download button implementation",
    decision = "Both — PNG for presentations, PDF for papers"
  ),

  d6 = list(
    question = "Will the app store user data between sessions?",
    options  = c("No — fresh upload each session",
                  "Yes — user accounts with saved analyses"),
    impact   = "Authentication, database, privacy implications",
    decision = "No for V1 — each session is independent"
  ),

  d7 = list(
    question = "How will PICRUSt2 functional data be handled?",
    options  = c("Upload pre-computed PICRUSt2 output",
                  "Run PICRUSt2 on server",
                  "Skip in V1"),
    impact   = "Server requirements, user workflow",
    decision = "Upload pre-computed output in V1 — running on server is V2"
  )
)

cat("=== Decision Points ===\n\n")
for (d in DECISIONS) {
  cat(sprintf("Q: %s\n", d$question))
  cat("  Options:\n")
  for (opt in d$options) cat("    -", opt, "\n")
  cat(sprintf("  Suggested decision: %s\n\n", d$decision))
}


# =============================================================================
# SECTION 9 — PHASE COMPLETION CHECKLIST
# =============================================================================

# Print this checklist and tick items off as you complete them.

CHECKLIST <- c(
  # Environment
  "[ ] R >= 4.2.0 installed",
  "[ ] All required CRAN packages installed",
  "[ ] All required Bioconductor packages installed",
  "[ ] renv initialised and snapshot saved",
  "[ ] renv.lock committed to Git",

  # Structure
  "[ ] Project folder structure created",
  "[ ] All 10 module R files in R/",
  "[ ] utils.R created and sourced by all modules",
  "[ ] validators.R created and tested",
  "[ ] precompute.R created",
  "[ ] pipeline_config.R created with your defaults",

  # Testing — simulated data
  "[ ] test_all_modules.R runs to completion",
  "[ ] 100% pass rate on simulated data",
  "[ ] test_report.html reviewed and saved",
  "[ ] session_info.txt saved",

  # Testing — public data
  "[ ] Public microbiome dataset downloaded",
  "[ ] Module 1 (QC) passes on real data",
  "[ ] Module 2 (Composition) biologically plausible",
  "[ ] Module 3 (Alpha) values in expected ranges",
  "[ ] Module 4 (Beta) groups visually separate",
  "[ ] Module 5 (DA) known taxa detected",
  "[ ] Module 6 (Functional) NSTI values acceptable",
  "[ ] Module 7 (Network) non-empty network produced",
  "[ ] Module 8 (ML) AUC > 0.6 on structured data",
  "[ ] Module 9 (Correlation) significant associations found",
  "[ ] Module 10 (Longitudinal) trajectories plotted",

  # Testing — your own data
  "[ ] Your data imports cleanly",
  "[ ] Full pipeline runs on your data without errors",
  "[ ] Results are biologically meaningful",
  "[ ] Parameter choices documented",

  # Performance
  "[ ] profvis run on all slow functions",
  "[ ] precompute.R tested — cache loads in < 1 second",
  "[ ] Memory use < 4 GB for full pipeline",
  "[ ] No single interactive function > 5 seconds",

  # Documentation and cleanup
  "[ ] All functions have @param documentation",
  "[ ] Debug cat() statements removed from production code",
  "[ ] TROUBLESHOOTING.md written",
  "[ ] Module contracts documented",
  "[ ] Decision points in Section 8 answered",
  "[ ] Clean Git commit with tag v0.1.0-pretesting",

  # Final gate
  "[ ] Colleague has run the pipeline on their machine",
  "[ ] Pipeline reproduces on a fresh R session (renv::restore())",
  "[ ] READY TO START SHINY INTEGRATION"
)

cat("\n=== Phase Completion Checklist ===\n")
cat(paste(CHECKLIST, collapse = "\n"))
cat("\n\n")

cat(strrep("=", 64), "\n")
cat("PRINT THIS CHECKLIST. TICK ITEMS OFF IN ORDER.\n")
cat("Do not start Shiny until the last item is checked.\n")
cat(strrep("=", 64), "\n")
