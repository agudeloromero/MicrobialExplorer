# ViromeAnalyst — Microbiome Analysis Platform
## Complete Pipeline Plan and Module Reference

---

> **Purpose of this document**
> This is the single reference for the complete ViromeAnalyst pipeline. It covers the project structure, every analysis module, the sequential testing plan, and what you need to do before starting Shiny integration. Keep it open while working.

---

## Table of Contents

1. [Project Overview](#1-project-overview)
2. [Folder Structure](#2-folder-structure)
3. [File Manifest](#3-file-manifest)
4. [Module Reference](#4-module-reference)
   - [Module 0 — Shared Utilities](#module-0--shared-utilities)
   - [Module 1 — Quality Control](#module-1--quality-control)
   - [Module 2 — Taxonomic Composition](#module-2--taxonomic-composition)
   - [Module 3 — Alpha Diversity](#module-3--alpha-diversity)
   - [Module 4 — Beta Diversity](#module-4--beta-diversity)
   - [Module 5 — Differential Abundance](#module-5--differential-abundance)
   - [Module 6 — Functional Prediction](#module-6--functional-prediction)
   - [Module 7 — Network Analysis](#module-7--network-analysis)
   - [Module 8 — Machine Learning](#module-8--machine-learning)
   - [Module 9 — Correlation and Association](#module-9--correlation-and-association)
   - [Module 10 — Longitudinal Analysis](#module-10--longitudinal-analysis)
5. [Supporting Scripts](#5-supporting-scripts)
6. [Package Installation](#6-package-installation)
7. [Sequential Testing Plan](#7-sequential-testing-plan)
8. [Performance Guidelines](#8-performance-guidelines)
9. [Common Errors and Fixes](#9-common-errors-and-fixes)
10. [Decision Points Before Shiny](#10-decision-points-before-shiny)
11. [Final Pre-Shiny Checklist](#11-final-pre-shiny-checklist)

---

## 1. Project Overview

ViromeAnalyst is a modular R-based microbiome and virome analysis platform. The pipeline is designed in two phases:

**Phase 1 (current) — Bacterial microbiome analysis (proof of concept)**
A complete, tested R pipeline covering quality control through machine learning, validated on real data before being integrated into a Shiny web application.

**Phase 2 — Virome-specific analysis**
A dedicated virome module (vOTU profiling, viral taxonomy, phage lifestyle prediction, AMG detection) added on top of the working bacterial foundation. This is what differentiates ViromeAnalyst from MicrobiomeAnalyst.

### Design Principles

- Every module is **self-contained**: each script can be sourced and run independently.
- Every module **returns objects**: no direct file writes inside functions. The caller (Shiny or the user) decides what to do with outputs.
- Every function **accepts a seed argument** for reproducibility.
- All modules **chain through a single filtered phyloseq RDS** object created by Module 1.
- **Cache-first**: expensive computations run once at data upload and are loaded from cache on all subsequent interactions.

---

## 2. Folder Structure

Run `create_project_structure("ViromeAnalyst")` from `setup_and_testing_plan.R` to create this automatically.

```
ViromeAnalyst/
│
├── R/                              # All analysis modules
│   ├── utils.R                     # Shared theme, colours, helpers
│   ├── validators.R                # Input validation functions
│   ├── precompute.R                # Cache management
│   ├── microbiome_qc.R             # Module 1
│   ├── microbiome_composition.R    # Module 2
│   ├── microbiome_alpha_diversity.R # Module 3
│   ├── microbiome_beta_diversity.R  # Module 4
│   ├── microbiome_differential_abundance.R # Module 5
│   ├── microbiome_functional.R     # Module 6
│   ├── microbiome_network.R        # Module 7
│   ├── microbiome_ml_classification.R # Module 8
│   ├── microbiome_correlation.R    # Module 9
│   └── microbiome_longitudinal.R   # Module 10
│
├── app/                            # Shiny application (Phase 2)
│   ├── ui/
│   ├── server/
│   └── modules/
│
├── config/
│   ├── pipeline_config.R           # All analysis defaults
│   └── ui_config.R                 # Shiny UI settings (Phase 2)
│
├── data/
│   ├── raw/                        # User data (gitignored)
│   ├── example/                    # Demo dataset
│   └── reference/                  # Reference databases (gitignored)
│
├── cache/                          # Pre-computed objects (gitignored)
│   ├── qc/
│   ├── composition/
│   ├── alpha/
│   ├── beta/
│   ├── differential/
│   ├── functional/
│   ├── network/
│   ├── ml/
│   ├── correlation/
│   └── longitudinal/
│
├── output/                         # Analysis outputs (gitignored)
│   ├── plots/
│   ├── tables/
│   └── reports/
│
├── tests/
│   ├── test_all_modules.R          # Master test suite
│   └── results/                    # Test output
│
├── docs/
│   └── VIROME_ANALYST_PLAN.md      # This document
│
├── logs/
├── renv.lock                       # Package versions (commit this)
├── .gitignore
└── README.md
```

---

## 3. File Manifest

| File | Status | Description |
|------|--------|-------------|
| `R/utils.R` | ✅ Write now | Shared theme, colours, helper functions |
| `R/validators.R` | ✅ Write now | Input validation for all data types |
| `R/precompute.R` | ✅ Write now | Cache management and pre-computation |
| `R/microbiome_qc.R` | ✅ Done | Quality control and filtering |
| `R/microbiome_composition.R` | ✅ Done | Taxonomic composition |
| `R/microbiome_alpha_diversity.R` | ✅ Done | Alpha diversity |
| `R/microbiome_beta_diversity.R` | ✅ Done | Beta diversity |
| `R/microbiome_differential_abundance.R` | ✅ Done | Differential abundance |
| `R/microbiome_functional.R` | ✅ Done | Functional prediction |
| `R/microbiome_network.R` | ✅ Done | Network analysis |
| `R/microbiome_ml_classification.R` | ✅ Done | Machine learning |
| `R/microbiome_correlation.R` | ✅ Done | Correlation and association |
| `R/microbiome_longitudinal.R` | ✅ Done | Longitudinal analysis |
| `config/pipeline_config.R` | ✅ Write now | Analysis defaults |
| `tests/test_all_modules.R` | ✅ Done | Master test suite |
| `setup_and_testing_plan.R` | ✅ Done | Setup utilities |

---

## 4. Module Reference

---

### Module 0 — Shared Utilities

**Script:** `R/utils.R`
**Status:** Write before sourcing any module
**Source line:** `source("R/utils.R")` — add to the top of every module

This file centralises all code shared across modules. It must exist before any module script is sourced. Defining shared elements in one place ensures visual consistency and avoids duplication.

**Contains:**

| Element | Type | Description |
|---------|------|-------------|
| `theme_microbiome()` | Function | Shared ggplot2 theme for all plots |
| `TAXA_PALETTE` | Vector | 21 colours for taxonomic bar plots |
| `KEGG_COLOURS` | Named vector | KEGG category colours |
| `pkg_available()` | Function | Silent package availability check |
| `` `%\|\|%` `` | Operator | Null-coalescing operator |
| `safe_run()` | Function | Error-catching wrapper with timing |
| `check_phyloseq()` | Function | Diagnostic summary of a phyloseq object |

**How to use:**
```r
source("R/utils.R")

# Check a phyloseq object before analysis
check_phyloseq(ps, label = "After QC")

# Safely run a module
result <- safe_run("Alpha diversity", quote(run_alpha_diversity(ps)))
```

---

### Module 1 — Quality Control

**Script:** `R/microbiome_qc.R`
**Input:** Raw OTU/ASV table, taxonomy table, metadata file (CSV/TSV)
**Output:** `cache/qc/phyloseq_qc_filtered.rds` — used by all downstream modules
**Run time:** 30 seconds – 3 minutes depending on dataset size

This is the entry point for the entire pipeline. Every downstream module reads the filtered phyloseq RDS produced here. The quality of this output determines the quality of all analyses that follow.

**Functions:**

| Function | Description |
|----------|-------------|
| `import_microbiome_data()` | Import OTU, taxonomy, and metadata into phyloseq |
| `qc_sequencing_depth()` | Assess and visualise read depth per sample |
| `qc_rarefaction_curves()` | Generate rarefaction curves |
| `qc_filter_taxa()` | Remove low-prevalence and low-abundance taxa |
| `qc_filter_samples()` | Remove samples below read/taxon thresholds |
| `qc_decontam()` | Detect and remove contaminant taxa |
| `qc_summary_report()` | Before/after filtering summary |
| `run_microbiome_qc()` | **Wrapper:** runs the complete QC pipeline |

**Typical usage:**
```r
source("R/utils.R")
source("R/microbiome_qc.R")

# Import data
ps_raw <- import_microbiome_data(
  otu_file  = "data/raw/otu_table.csv",
  tax_file  = "data/raw/taxonomy.csv",
  meta_file = "data/raw/metadata.csv"
)

# Run full QC
ps_clean <- run_microbiome_qc(
  otu_file   = "data/raw/otu_table.csv",
  tax_file   = "data/raw/taxonomy.csv",
  meta_file  = "data/raw/metadata.csv",
  min_reads  = 5000,
  group_var  = "disease_status",
  output_dir = "cache/qc"
)
```

**Key outputs saved:**
- `cache/qc/phyloseq_qc_filtered.rds` — the filtered phyloseq
- `cache/qc/01_sequencing_depth.pdf`
- `cache/qc/02_rarefaction_curves.pdf`
- `cache/qc/03_taxa_prevalence_abundance.pdf`

**What to check:** Median reads per sample, samples flagged below threshold, taxa removed percentage (aim for < 50% removed), rarefaction curves should plateau.

---

### Module 2 — Taxonomic Composition

**Script:** `R/microbiome_composition.R`
**Input:** `ps_clean` from Module 1
**Output:** Stacked bars, heatmaps, core microbiome, F:B ratio, abundance tables

This module answers the question: *what is in my samples, and how does composition differ between groups?* It is the most visually rich module and typically produces the figures that go directly into publications.

**Functions:**

| Function | Description |
|----------|-------------|
| `agglomerate_taxa()` | Collapse to a rank, apply transformation, merge rare taxa as "Other" |
| `plot_composition_bars()` | Stacked bar chart per sample with group annotation strip |
| `plot_mean_composition()` | Grouped bars with mean ± SE per group |
| `plot_abundance_heatmap()` | Hierarchically clustered heatmap with group annotation |
| `make_abundance_table()` | Tidy summary table with mean, SD, median, prevalence |
| `plot_phylum_drilldown()` | Genus-level breakdown within a selected phylum |
| `calculate_fb_ratio()` | Firmicutes:Bacteroidota ratio per sample |
| `identify_core_microbiome()` | Taxa present in ≥ X% of samples at multiple thresholds |
| `run_composition_analysis()` | **Wrapper:** runs the complete composition pipeline |

**Typical usage:**
```r
source("R/utils.R")
source("R/microbiome_composition.R")

ps <- readRDS("cache/qc/phyloseq_qc_filtered.rds")

results <- run_composition_analysis(
  ps        = ps,
  group_var = "disease_status",
  output_dir = "cache/composition"
)
```

**Key outputs:**
- Phylum and Genus level stacked bars
- CLR-transformed genus heatmap
- Core microbiome scatter plot (prevalence vs abundance)
- Firmicutes:Bacteroidota ratio violin plot
- `abundance_table_overall.csv` — supplementary data ready

---

### Module 3 — Alpha Diversity

**Script:** `R/microbiome_alpha_diversity.R`
**Input:** `ps_clean` from Module 1 (must be **raw counts**, not relative abundance)
**Output:** Diversity metrics per sample, statistical tests, rarefaction curves

Alpha diversity quantifies diversity *within* each sample. This module computes nine metrics simultaneously, averaging across multiple rarefaction iterations for stable estimates. Non-parametric tests are used throughout as microbiome diversity data is rarely normally distributed.

**Functions:**

| Function | Description |
|----------|-------------|
| `calculate_alpha_diversity()` | Compute Observed, Chao1, ACE, Shannon, Simpson, InvSimpson, Pielou, Fisher, Faith PD |
| `plot_alpha_diversity()` | Boxplots per metric with significance brackets |
| `plot_rarefaction_curves()` | Per-sample curves with SE ribbon and plateau detection |
| `plot_diversity_gradient()` | Diversity vs continuous metadata variable (Spearman) |
| `plot_longitudinal_diversity()` | Trajectories over time (uses longitudinal module data) |
| `summarise_alpha_stats()` | Publication-ready statistical summary table |
| `plot_diversity_correlations()` | Pairwise correlation matrix of all metrics |
| `run_alpha_diversity()` | **Wrapper:** runs the complete alpha diversity pipeline |

**Typical usage:**
```r
source("R/utils.R")
source("R/microbiome_alpha_diversity.R")

ps <- readRDS("cache/qc/phyloseq_qc_filtered.rds")

results <- run_alpha_diversity(
  ps         = ps,
  group_var  = "disease_status",
  rare_depth = 10000,
  output_dir = "cache/alpha"
)

# Access the diversity table directly
head(results$diversity_df)
```

**Important:** Always use raw counts (not relative abundance) as input. The function will error if relative abundance is detected.

---

### Module 4 — Beta Diversity

**Script:** `R/microbiome_beta_diversity.R`
**Input:** `ps_clean` from Module 1
**Output:** Distance matrices, ordination plots, PERMANOVA results, betadisper

Beta diversity quantifies differences in community composition *between* samples. This module is the statistical backbone for group comparison and is the source of the PCoA plots that appear in most microbiome papers.

**Functions:**

| Function | Description |
|----------|-------------|
| `compute_distances()` | Bray-Curtis, Jaccard, UniFrac, Weighted UniFrac, Aitchison |
| `run_ordination()` | PCoA, NMDS, or PCA with variance explained |
| `plot_ordination()` | Full ordination plot with ellipses, spiders, shape variables |
| `plot_pcoa_panel()` | Three-panel PCoA showing all axis pairs |
| `run_permanova()` | adonis2 with covariates, pairwise tests, betadisper |
| `plot_betadisper()` | Distance-to-centroid PCoA and boxplot |
| `plot_distance_heatmap()` | Hierarchically clustered pairwise distance heatmap |
| `plot_envfit()` | Metadata vectors overlaid on ordination |
| `run_mantel_test()` | Correlation between two distance matrices |
| `run_beta_diversity()` | **Wrapper:** runs the complete beta diversity pipeline |

**Typical usage:**
```r
source("R/utils.R")
source("R/microbiome_beta_diversity.R")

ps <- readRDS("cache/qc/phyloseq_qc_filtered.rds")

results <- run_beta_diversity(
  ps          = ps,
  group_var   = "disease_status",
  formula_rhs = "disease_status + age + sex",
  methods     = c("bray", "jaccard"),
  env_vars    = c("age", "bmi", "crp"),
  output_dir  = "cache/beta"
)

# PERMANOVA result
results$permanova$permanova
results$permanova$pairwise
```

**Always run betadisper alongside PERMANOVA.** If betadisper is significant, PERMANOVA differences may reflect dispersion rather than location — this must be acknowledged in publications.

---

### Module 5 — Differential Abundance

**Script:** `R/microbiome_differential_abundance.R`
**Input:** `ps_clean` from Module 1 (raw counts)
**Output:** DA results from multiple methods, consensus table, volcano plots, effect size plots

Identifying which specific taxa differ between groups is the most common analysis goal in microbiome research. This module runs three complementary methods and builds a consensus — taxa identified by multiple methods are more reliable than single-method results.

**Functions:**

| Function | Description |
|----------|-------------|
| `prepare_da_data()` | Agglomerate, filter, handle duplicates |
| `run_ancombc()` | ANCOM-BC2 with covariates (primary recommended method) |
| `run_deseq2()` | DESeq2 with poscounts size factors and LFC shrinkage |
| `run_aldex2()` | ALDEx2 Dirichlet-multinomial model |
| `build_consensus()` | Score taxa by agreement across methods |
| `plot_volcano()` | Volcano plot with top taxa labelled |
| `plot_effect_sizes()` | Ranked horizontal bar chart with error bars |
| `plot_da_heatmap()` | Z-score heatmap of significant taxa |
| `plot_da_boxplots()` | Individual boxplots for top DA taxa |
| `run_differential_abundance()` | **Wrapper:** runs all methods and consensus |

**Typical usage:**
```r
source("R/utils.R")
source("R/microbiome_differential_abundance.R")

ps <- readRDS("cache/qc/phyloseq_qc_filtered.rds")

results <- run_differential_abundance(
  ps            = ps,
  group_var     = "disease_status",
  formula       = "disease_status + age + sex",
  reference     = "Healthy",
  rank          = "Genus",
  alpha         = 0.05,
  lfc_threshold = 1,
  min_methods   = 2,
  output_dir    = "cache/differential"
)

# Consensus significant taxa
results$consensus$consensus %>% filter(consensus_da)
```

**Requires Bioconductor packages:** ANCOMBC, DESeq2 (BiocManager::install). ALDEx2 is optional. If packages are missing, functions fail gracefully with informative messages.

---

### Module 6 — Functional Prediction

**Script:** `R/microbiome_functional.R`
**Input:** PICRUSt2 output directory + sample metadata
**Output:** Pathway abundance plots, KEGG categories, functional DA, cross-domain heatmaps

This module processes the output of PICRUSt2 (which must be run externally on the command line) to translate 16S marker gene data into predicted functional profiles. It also supports cross-domain analysis linking microbiome to metabolomics or other external datasets.

**Functions:**

| Function | Description |
|----------|-------------|
| `import_picrust2()` | Import KO, EC, pathway, COG tables from PICRUSt2 output |
| `picrust2_to_phyloseq()` | Convert any PICRUSt2 matrix to phyloseq object |
| `plot_nsti_quality()` | Assess prediction reliability by NSTI score |
| `plot_pathway_abundance()` | Mean ± SE pathway abundances per group |
| `summarise_kegg_categories()` | Roll up KOs into KEGG Level 1 categories |
| `analyse_functional_diversity()` | Alpha and beta diversity of functional profiles |
| `test_functional_da()` | Wilcoxon/Kruskal-Wallis tests on pathway/KO abundances |
| `plot_functional_heatmap()` | Clustered heatmap of top functional features |
| `analyse_pathway_contributors()` | Which taxa drive differences in target pathways? |
| `run_functional_analysis()` | **Wrapper:** runs the complete functional pipeline |

**Typical usage:**
```r
source("R/utils.R")
source("R/microbiome_functional.R")

results <- run_functional_analysis(
  picrust_dir = "data/picrust2_output/",
  meta_file   = "data/raw/metadata.csv",
  group_var   = "disease_status",
  alpha       = 0.05,
  output_dir  = "cache/functional"
)
```

**Note:** NSTI > 0.15 indicates unreliable predictions. Always run `plot_nsti_quality()` first and exclude high-NSTI samples before interpreting results.

---

### Module 7 — Network Analysis

**Script:** `R/microbiome_network.R`
**Input:** `ps_clean` from Module 1
**Output:** Co-occurrence network, hub taxa, community modules, differential networks

Co-occurrence network analysis reveals which taxa tend to co-occur (positive edges) or exclude each other (negative edges). Hub taxa — those with many connections — are often keystone species that have outsized influence on the community.

**Functions:**

| Function | Description |
|----------|-------------|
| `prepare_network_data()` | Filter to prevalent taxa, cap at max_taxa |
| `build_network()` | Spearman/Pearson with BH correction, or SPIEC-EASI |
| `calculate_network_topology()` | Degree, betweenness, eigenvector, hub score, modularity |
| `plot_network()` | ggraph visualisation with multiple layout options |
| `plot_network_phylum()` | Cleaner phylum-coloured layout for publications |
| `analyse_hub_taxa()` | Composite hub score ranking and centrality scatter |
| `analyse_modules()` | Louvain community detection and phylum composition |
| `compare_networks()` | Differential network analysis between two groups |
| `compare_random_network()` | Small-world test against Erdős-Rényi random networks |
| `run_network_analysis()` | **Wrapper:** runs the complete network pipeline |

**Typical usage:**
```r
source("R/utils.R")
source("R/microbiome_network.R")

ps <- readRDS("cache/qc/phyloseq_qc_filtered.rds")

results <- run_network_analysis(
  ps            = ps,
  group_var     = "disease_status",
  rank          = "Genus",
  method        = "spearman",
  cor_threshold = 0.6,
  output_dir    = "cache/network"
)

# Hub taxa
results$topology$nodes %>% filter(is_hub)
# Small-world index
results$sw_index
```

**Tip:** If the network has 0 edges, lower `cor_threshold` (try 0.4) or reduce `p_threshold` (try 0.10). If the network is too dense (density > 0.3), raise `cor_threshold`.

---

### Module 8 — Machine Learning

**Script:** `R/microbiome_ml_classification.R`
**Input:** `ps_clean` from Module 1
**Output:** Trained models, ROC curves, feature importance, confusion matrices, biomarker panel

Machine learning classification identifies which combination of taxa best predicts group membership. This module uses two complementary approaches: Random Forest (powerful, handles non-linearity) and LASSO (sparse, automatic biomarker selection). Consensus features identified by both are the most robust biomarker candidates.

**Functions:**

| Function | Description |
|----------|-------------|
| `prepare_ml_data()` | CLR transform, prevalence filter, stratified train/test split |
| `train_random_forest()` | Repeated CV with mtry tuning, AUC/accuracy optimisation |
| `train_lasso()` | 10-fold CV for lambda, extracts non-zero coefficients |
| `plot_feature_importance()` | RF importance + LASSO coefficients + consensus panel |
| `plot_roc_curves()` | Overlaid ROC with 95% CI on AUC for multiple models |
| `plot_confusion_matrix()` | Heatmap with row-percentages and performance metrics |
| `plot_learning_curve()` | Training vs CV score to diagnose overfitting |
| `assess_feature_stability()` | Bootstrap importance distributions |
| `plot_multiclass_roc()` | One-vs-rest ROC for > 2 groups |
| `run_ml_classification()` | **Wrapper:** runs the complete ML pipeline |

**Typical usage:**
```r
source("R/utils.R")
source("R/microbiome_ml_classification.R")

ps <- readRDS("cache/qc/phyloseq_qc_filtered.rds")

results <- run_ml_classification(
  ps         = ps,
  group_var  = "disease_status",
  rank       = "Genus",
  n_trees    = 500,
  cv_folds   = 5,
  output_dir = "cache/ml"
)

# AUC
results$rf$confusion$overall["Accuracy"]
# Consensus biomarkers
results$importance$importance$rf %>% filter(is_hub)
```

**Requires:** `randomForest`, `caret`, `pROC`, `glmnet` (all CRAN).

---

### Module 9 — Correlation and Association

**Script:** `R/microbiome_correlation.R`
**Input:** `ps_clean` from Module 1 + optional external data matrix
**Output:** Taxa-metadata correlations, taxa-taxa correlations, cross-domain heatmaps

This module reveals which taxa are associated with clinical or biological metadata variables (age, BMI, disease scores), which taxa co-occur, and — critically for virome research — how microbiome profiles correlate with external datasets such as metabolomics or cytokines.

**Functions:**

| Function | Description |
|----------|-------------|
| `correlate_taxa_metadata()` | Spearman rho for every taxon-metadata pair with BH correction |
| `plot_taxa_metadata_heatmap()` | Clustered heatmap with significance asterisks |
| `plot_top_associations()` | Scatter plots for top significant associations |
| `plot_association_bubbles()` | Bubble chart of all associations (size = \|rho\|) |
| `compute_taxa_taxa_correlation()` | Pairwise taxa correlation matrix |
| `plot_metadata_correlations()` | Correlation matrix of all numeric metadata |
| `correlate_cross_domain()` | Microbiome × external data (metabolomics, cytokines) |
| `test_mixed_effects()` | lme4 models for repeated-measures designs |
| `run_correlation_analysis()` | **Wrapper:** runs the complete correlation pipeline |

**Typical usage:**
```r
source("R/utils.R")
source("R/microbiome_correlation.R")

ps <- readRDS("cache/qc/phyloseq_qc_filtered.rds")

results <- run_correlation_analysis(
  ps         = ps,
  group_var  = "disease_status",
  meta_vars  = c("age", "bmi", "crp", "calprotectin"),
  rank       = "Genus",
  alpha      = 0.05,
  output_dir = "cache/correlation"
)

# Significant associations
results$cor_taxa_meta %>% filter(significant) %>% head(20)
```

**Cross-domain example:**
```r
metabolites <- read.csv("data/metabolomics.csv", row.names = 1)

cross_res <- correlate_cross_domain(
  ps           = ps,
  external_mat = as.matrix(metabolites),
  rank         = "Genus",
  alpha        = 0.05
)
```

---

### Module 10 — Longitudinal Analysis

**Script:** `R/microbiome_longitudinal.R`
**Input:** `ps_clean` from Module 1 — **must have time and subject metadata columns**
**Output:** Diversity trajectories, community movement plots, stability analysis, LME models

Longitudinal analysis reveals how microbiome communities change over time within individuals. It requires at least a `timepoint` (numeric) and `subject_id` column in the metadata. This module is essential for intervention studies, clinical trials, and disease progression research.

**Functions:**

| Function | Description |
|----------|-------------|
| `prepare_longitudinal_data()` | Validate time/subject columns, compute per-sample diversity |
| `plot_diversity_trajectories()` | Individual lines + group mean ± SE + significance per timepoint |
| `plot_beta_trajectories()` | PCoA with arrows + distance-to-baseline over time |
| `analyse_stability()` | Consecutive-timepoint dissimilarity, violin and ranked bar plots |
| `plot_composition_over_time()` | Stacked area + relative change from baseline |
| `run_lme_over_time()` | lme4 models for every taxon, forest plot of time effects |
| `analyse_intervention_response()` | Pre/post comparison, responder classification |
| `detect_changepoints()` | Rolling z-score to flag abrupt community shifts |
| `run_longitudinal_analysis()` | **Wrapper:** runs the complete longitudinal pipeline |

**Typical usage:**
```r
source("R/utils.R")
source("R/microbiome_longitudinal.R")

ps <- readRDS("cache/qc/phyloseq_qc_filtered.rds")

results <- run_longitudinal_analysis(
  ps                 = ps,
  time_var           = "week",
  subject_var        = "patient_id",
  group_var          = "treatment",
  pre_timepoints     = c(0, 1),
  post_timepoints    = c(4, 8, 12),
  intervention_label = "Antibiotic course",
  output_dir         = "cache/longitudinal"
)
```

---

## 5. Supporting Scripts

### `R/utils.R`
Contains `theme_microbiome()`, `TAXA_PALETTE`, `KEGG_COLOURS`, `pkg_available()`, `` `%||%` ``, `safe_run()`, `check_phyloseq()`. **Must be sourced before any module.**

### `R/validators.R`
Contains `validate_otu_table()`, `validate_metadata()`, `validate_phyloseq()`, `check_contract()`. Used by Shiny to give informative error messages when users upload incorrect data.

### `R/precompute.R`
Contains `precompute_all()` and `load_cache()`. Runs expensive operations once at data upload and saves results to `cache/`. Shiny loads from cache instead of recomputing.

### `config/pipeline_config.R`
Single source of truth for all analysis defaults. Change parameters here, not in the module files. Contains `PIPELINE_CONFIG` list with settings for QC, beta diversity, differential abundance, ML, network, and visualisation.

### `tests/test_all_modules.R`
Comprehensive test suite that:
- Creates a simulated 60-sample × 150-taxa dataset
- Runs 60+ individual tests across all 10 modules
- Reports PASS/FAIL with timing for every test
- Generates `test_report.html` for human review
- Saves `session_info.txt` for reproducibility
- Returns a summary list for programmatic use

### `setup_and_testing_plan.R`
The master setup script. Contains `create_project_structure()` to build the folder tree, the full `FILE_MANIFEST`, `TESTING_PHASES`, `PACKAGES` lists, and the `KNOWN_ERRORS` reference.

---

## 6. Package Installation

### Step 1 — Install BiocManager
```r
install.packages("BiocManager")
```

### Step 2 — Required Bioconductor packages
```r
BiocManager::install(c(
  "phyloseq",
  "microbiome",
  "decontam",
  "Hmisc"
))
```

### Step 3 — Required CRAN packages
```r
install.packages(c(
  "ggplot2", "dplyr", "tidyr", "patchwork", "vegan",
  "scales", "tibble", "RColorBrewer", "stringr", "forcats",
  "ggrepel", "zoo", "igraph", "ggraph", "tidygraph",
  "randomForest", "caret", "pROC", "e1071", "glmnet",
  "rstatix", "ggpubr", "corrplot", "profvis", "renv"
))
```

### Step 4 — Optional packages (recommended for full functionality)
```r
# Differential abundance methods
BiocManager::install(c("DESeq2", "ANCOMBC", "ALDEx2"))

# Mixed-effects models
install.packages(c("lme4", "lmerTest", "broom.mixed"))

# Network (SPIEC-EASI)
remotes::install_github("zdk123/SpiecEasi")

# Phylogenetic diversity
BiocManager::install("picante")

# Profiling and testing
install.packages(c("pryr", "testthat"))
```

### Step 5 — Lock environment
```r
renv::init()
renv::snapshot()
```

---

## 7. Sequential Testing Plan

Work through these phases **in order**. Do not move to the next phase until all pass criteria for the current phase are met.

---

### Phase 1 — Environment Setup
**Time estimate:** 30–60 minutes
**Pass criteria:** `renv::status()` reports no issues

| Step | Action |
|------|--------|
| 1.1 | Verify R >= 4.2.0: `R.version$major` |
| 1.2 | `install.packages("renv")` |
| 1.3 | `renv::init()` |
| 1.4 | Install all required packages (see Section 6) |
| 1.5 | Install optional packages |
| 1.6 | `renv::snapshot()` |
| 1.7 | `renv::status()` — confirm no issues |

---

### Phase 2 — File and Structure Verification
**Time estimate:** 1–2 hours
**Pass criteria:** All 14 `source()` calls succeed with no errors

```r
# Source each module file individually and check for errors
files_to_test <- c(
  "R/utils.R",
  "R/validators.R",
  "R/precompute.R",
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
  "config/pipeline_config.R"
)

for (f in files_to_test) {
  cat("Sourcing", f, "... ")
  tryCatch({ source(f); cat("OK\n") },
           error = function(e) cat("FAIL:", e$message, "\n"))
}
```

---

### Phase 3 — Simulated Data Test
**Time estimate:** 2–4 hours
**Pass criteria:** 100% pass rate in `test_report.html`

```r
source("tests/test_all_modules.R")
# Open: tests/results/test_report.html
```

---

### Phase 4 — Real Data Test (Public Dataset)
**Time estimate:** 4–8 hours
**Pass criteria:** Results are biologically plausible

```r
# Download a well-annotated public dataset
BiocManager::install("curatedMetagenomicData")
library(curatedMetagenomicData)

# Use HMP IBD dataset
ps_public <- curatedMetagenomicData(
  "HMP_2019_ibdmdb.relative_abundance",
  dryrun = FALSE
) |> mergeData()

# Run each module sequentially
ps_qc    <- run_microbiome_qc(...)
run_composition_analysis(ps_qc, ...)
run_alpha_diversity(ps_qc, ...)
run_beta_diversity(ps_qc, ...)
run_differential_abundance(ps_qc, ...)
# etc.
```

**Check that results are consistent with the published paper that used this dataset.**

---

### Phase 5 — Your Own Data Test
**Time estimate:** Full day
**Pass criteria:** Full pipeline completes on your data with no errors

```r
# Validate your files first
source("R/validators.R")
otu_mat  <- read.csv("data/raw/otu.csv", row.names = 1)
meta_df  <- read.csv("data/raw/metadata.csv", row.names = 1)

v1 <- validate_otu_table(as.matrix(otu_mat))
v2 <- validate_metadata(meta_df,
                          required_vars = c("disease_status"))
if (!v1$valid) stop(paste(v1$errors, collapse="\n"))
if (!v2$valid) stop(paste(v2$errors, collapse="\n"))
if (length(v2$warnings) > 0) cat("Warnings:\n", paste(v2$warnings, collapse="\n"))
```

---

### Phase 6 — Performance Profiling
**Time estimate:** 2–4 hours
**Pass criteria:** All interactive functions run in < 5 seconds with cache

```r
library(profvis)

# Profile the slowest module
profvis({
  run_beta_diversity(ps_clean, group_var = "disease_status",
                      methods = c("bray", "jaccard"),
                      permutations = 999)
})

# Profile the pre-computation step
profvis({
  precompute_all(ps_clean, group_var = "disease_status",
                  output_dir = "cache/")
})
```

**Slow functions to watch:** `tax_glom()`, `metaMDS()`, `adonis2()`, `train_random_forest()`, `build_network()`, ANCOM-BC2.

---

### Phase 7 — Documentation and Cleanup
**Time estimate:** 2–3 hours
**Pass criteria:** Clean Git commit with tag `v0.1.0-pretesting`

```r
# Final snapshot
renv::snapshot()

# Final test run
source("tests/test_all_modules.R")

# Commit
# git add -A
# git commit -m "Pre-Shiny pipeline: all 10 modules tested, 100% pass rate"
# git tag v0.1.0-pretesting
```

---

## 8. Performance Guidelines

### Target response times for Shiny

| Operation | Target | Strategy |
|-----------|--------|----------|
| Data upload + pre-computation | < 60 sec | Run `precompute_all()` once, cache to RDS |
| Ordination plot | < 1 sec | Load pre-computed PCoA coords from cache |
| PERMANOVA | < 2 sec | Reduce to 499 perms in interactive mode |
| Composition bar plot | < 0.5 sec | Pre-agglomerated phyloseq in cache |
| Alpha diversity plot | < 0.5 sec | Load diversity_df from cache |
| Differential abundance | < 5 sec | Load DA results from cache |
| Random Forest | < 3 sec | Reduce to 200 trees for interactive |
| Network | < 3 sec | Pre-compute graph, only re-layout |

### Memory guidelines

```r
# Check memory use after loading data
library(pryr)
mem_used()  # Should be < 2 GB after pre-computation

# Check object sizes
object.size(ps_clean)   |> format(units = "MB")
object.size(distances)  |> format(units = "MB")
```

**Target:** Total memory use < 4 GB for a 500-sample dataset on a 16 GB RAM server.

---

## 9. Common Errors and Fixes

| Error | Cause | Fix |
|-------|-------|-----|
| `NAs in taxonomic rank` | NA at the requested rank in taxonomy table | Add `NArm = FALSE` to `tax_glom()` |
| `must have at least 2 groups` | A filter removed all samples in one group | Check `table(sample_data(ps)$group_var)` after filtering |
| `object of class 'dist' expected` | Sample order mismatch | `meta_df <- meta_df[labels(dist_obj), , drop = FALSE]` |
| `NaN after CLR` | `log(0)` is undefined | Always use `otu + 0.5` pseudocount before CLR |
| `number of obs <= random effects` | Too few samples per subject | Require `min_timepoints = 3` in longitudinal module |
| `Can't have empty classes` (RF) | One class has 0 samples in training | Set `balance_classes = TRUE` |
| `rownames not matching` (DESeq2) | Sample name format mismatch | `ps <- prune_samples(sample_names(ps), ps)` |
| `network has 0 edges` | Threshold too strict | Lower `cor_threshold` to 0.4 |
| `no shared samples` (cross-domain) | Name format difference | `trimws(colnames(external_mat))` |
| `package ANCOMBC not found` | Not installed | `BiocManager::install("ANCOMBC")` |

---

## 10. Decision Points Before Shiny

Answer these questions and document your answers in `config/pipeline_config.R`:

| Decision | Options | Recommendation |
|----------|---------|---------------|
| Primary taxonomic rank | Genus, Species, Family | **Genus** — best balance of resolution and stability |
| Primary DA method | ANCOM-BC2, DESeq2, Consensus | **Consensus** — more reliable, show agreement |
| Max dataset size V1 | 200, 500, 1000 samples | **500** — matches server specs |
| Virome module placement | Separate tab, integrated, teased | **Separate tab** — teased in V1, active in V2 |
| Download formats | PNG, PDF, both | **Both** — PNG for presentations, PDF for papers |
| Session persistence | Fresh each session, user accounts | **Fresh sessions** for V1 |
| PICRUSt2 handling | Upload output, run on server, skip | **Upload pre-computed** for V1 |

---

## 11. Final Pre-Shiny Checklist

Print this list and tick items off in order. Do not start Shiny until the last item is checked.

### Environment
- [ ] R >= 4.2.0 installed
- [ ] All required CRAN packages installed and loading
- [ ] All required Bioconductor packages installed and loading
- [ ] renv initialised and `renv::snapshot()` completed
- [ ] `renv.lock` committed to Git

### Structure and Files
- [ ] Project folder structure created
- [ ] All 10 module R files in `R/`
- [ ] `utils.R` created and sourced by all modules
- [ ] `validators.R` written and tested
- [ ] `precompute.R` written and tested
- [ ] `pipeline_config.R` created with your analysis defaults

### Testing — Simulated Data
- [ ] `tests/test_all_modules.R` runs to completion without crashing
- [ ] 100% pass rate on simulated data
- [ ] `test_report.html` reviewed — all green
- [ ] `session_info.txt` saved

### Testing — Public Data
- [ ] Public microbiome dataset downloaded
- [ ] Module 1 (QC) passes and output makes sense
- [ ] Module 2 (Composition) plots are biologically plausible
- [ ] Module 3 (Alpha) diversity values in expected ranges
- [ ] Module 4 (Beta) groups separate visually
- [ ] Module 5 (DA) detects known differentially abundant taxa
- [ ] Module 6 (Functional) NSTI values < 0.15 for most samples
- [ ] Module 7 (Network) produces a non-empty network
- [ ] Module 8 (ML) AUC > 0.6 on the structured public dataset
- [ ] Module 9 (Correlation) significant associations found
- [ ] Module 10 (Longitudinal) trajectories plotted (if time data available)

### Testing — Your Own Data
- [ ] Your data imports cleanly through `validate_*()` functions
- [ ] Full pipeline runs on your data without errors
- [ ] Results are biologically meaningful and consistent with your knowledge
- [ ] All parameter choices documented in `pipeline_config.R`

### Performance
- [ ] `profvis` run on all modules
- [ ] `precompute_all()` tested — cache loads in < 1 second
- [ ] Memory use < 4 GB for full pipeline on your data
- [ ] No single interactive operation > 5 seconds with cache

### Documentation and Cleanup
- [ ] All exported functions have `@param` documentation
- [ ] Debug `cat()` statements removed from production functions
- [ ] Decision points in Section 10 answered and documented
- [ ] Final `renv::snapshot()` completed
- [ ] All changes committed to Git
- [ ] Pipeline reproduces from scratch: `renv::restore()` on a new session
- [ ] **Git tag created:** `git tag v0.1.0-pretesting`

---

> **✓ READY TO START SHINY INTEGRATION**
> When all items above are checked, you have a verified, documented, reproducible analysis pipeline. The Shiny app is a user interface layer on top of this pipeline — not a rewrite of it.

---

*Last updated: May 2026 | ViromeAnalyst v0.1.0*
