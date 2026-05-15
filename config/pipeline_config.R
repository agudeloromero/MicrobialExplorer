# =============================================================================
# pipeline_config.R — Global analysis parameters for ViromeAnalyst
# =============================================================================
# This is the SINGLE SOURCE OF TRUTH for all analysis defaults.
# Change parameters here, not inside the module scripts.
#
# Usage:
#   source("config/pipeline_config.R")
#   PIPELINE_CONFIG$qc$min_reads
# =============================================================================

PIPELINE_CONFIG <- list(

  # ---------------------------------------------------------------------------
  # QC PARAMETERS (Module 1)
  # ---------------------------------------------------------------------------
  qc = list(
    min_reads      = 1000,      # Minimum reads per sample
    min_taxa       = 10,        # Minimum taxa per sample
    min_prevalence = 0.05,      # Minimum fraction of samples a taxon must appear in
    min_abundance  = 10,        # Minimum total reads across all samples
    rarefaction    = TRUE       # Whether to rarefy before diversity analysis
  ),

  # ---------------------------------------------------------------------------
  # TAXONOMIC RANK (used as default across all modules)
  # ---------------------------------------------------------------------------
  primary_rank = "Genus",       # Change to "Species" for higher resolution
                                 # or "Family" for lower resolution

  # ---------------------------------------------------------------------------
  # ABUNDANCE TRANSFORMATION
  # ---------------------------------------------------------------------------
  transform = "clr",            # "clr" (recommended), "relative", or "log10"

  # ---------------------------------------------------------------------------
  # BETA DIVERSITY (Module 4)
  # ---------------------------------------------------------------------------
  beta = list(
    methods      = c("bray", "jaccard"),  # Distance methods to compute
    permutations = 999,                    # PERMANOVA permutations (use 99 for testing)
    rare_depth   = NULL                    # NULL = auto (minimum sample depth)
  ),

  # ---------------------------------------------------------------------------
  # DIFFERENTIAL ABUNDANCE (Module 5)
  # ---------------------------------------------------------------------------
  differential = list(
    alpha         = 0.05,       # Significance threshold
    lfc_threshold = 1,          # Minimum absolute log fold change
    min_methods   = 2,          # Minimum DA methods for consensus
    p_adjust      = "BH",       # P-value adjustment method
    min_prevalence = 0.10,      # Min prevalence for DA testing
    min_count      = 10         # Min total count for DA testing
  ),

  # ---------------------------------------------------------------------------
  # MACHINE LEARNING (Module 8)
  # ---------------------------------------------------------------------------
  ml = list(
    test_fraction   = 0.2,      # Fraction of samples for test set
    n_trees         = 500,      # Random Forest trees
    n_trees_fast    = 200,      # Trees for interactive/Shiny use
    cv_folds        = 5,        # Cross-validation folds
    cv_repeats      = 3,        # CV repetitions
    n_boot          = 50,       # Bootstrap iterations for stability
    max_features    = 200       # Max features before variance-based selection
  ),

  # ---------------------------------------------------------------------------
  # NETWORK ANALYSIS (Module 7)
  # ---------------------------------------------------------------------------
  network = list(
    method         = "spearman",  # "spearman", "pearson", or "spiec-easi"
    cor_threshold  = 0.6,         # Minimum absolute correlation for an edge
    p_threshold    = 0.05,        # Maximum adjusted p-value for an edge
    min_prevalence = 0.20,        # Minimum taxon prevalence for network
    max_taxa       = 100          # Maximum taxa in network (by prevalence)
  ),

  # ---------------------------------------------------------------------------
  # CORRELATION ANALYSIS (Module 9)
  # ---------------------------------------------------------------------------
  correlation = list(
    method         = "spearman",  # Correlation method
    p_adjust       = "BH",        # Adjustment method
    alpha          = 0.05,        # Significance threshold
    min_prevalence = 0.10         # Minimum taxon prevalence
  ),

  # ---------------------------------------------------------------------------
  # LONGITUDINAL ANALYSIS (Module 10)
  # ---------------------------------------------------------------------------
  longitudinal = list(
    min_timepoints = 2,           # Minimum time points per subject
    min_subjects   = 5,           # Minimum subjects
    distance       = "bray"       # Distance for beta trajectory analysis
  ),

  # ---------------------------------------------------------------------------
  # VISUALISATION DEFAULTS
  # ---------------------------------------------------------------------------
  viz = list(
    top_n_taxa      = 20,         # Taxa in composition bars
    top_n_features  = 30,         # Features in heatmaps and importance plots
    top_n_da        = 25,         # Top DA taxa in effect size plots
    top_n_hub       = 20,         # Hub taxa to display
    figure_width    = 12,         # Default PDF width (inches)
    figure_height   = 8,          # Default PDF height (inches)
    dpi             = 300,        # Resolution for PNG export
    palette_default = "Set2",     # RColorBrewer palette for groups
    n_label_hubs    = TRUE        # Whether to label hub taxa in network
  ),

  # ---------------------------------------------------------------------------
  # SERVER AND PERFORMANCE
  # ---------------------------------------------------------------------------
  server = list(
    max_samples     = 500,        # Maximum samples accepted
    max_taxa        = 5000,       # Maximum taxa accepted
    cache_dir       = "cache/",   # Cache directory
    output_dir      = "output/",  # Output directory
    permutations_interactive = 99, # PERMANOVA perms for Shiny (fast)
    permutations_full        = 999 # PERMANOVA perms for download (accurate)
  ),

  # ---------------------------------------------------------------------------
  # REPRODUCIBILITY
  # ---------------------------------------------------------------------------
  seed = 42
)

# Print summary on load
cat("[pipeline_config.R] Configuration loaded.\n")
cat(sprintf("  Primary rank:  %s\n", PIPELINE_CONFIG$primary_rank))
cat(sprintf("  Transform:     %s\n", PIPELINE_CONFIG$transform))
cat(sprintf("  Alpha:         %.2f\n", PIPELINE_CONFIG$differential$alpha))
cat(sprintf("  Max samples:   %d\n", PIPELINE_CONFIG$server$max_samples))
cat(sprintf("  Seed:          %d\n", PIPELINE_CONFIG$seed))
