# MicrobialExplorer — Shiny App Design Notes

## Compiled from development and testing sessions (May 2026)
## Updated after Phase 4 validation across 5 datasets

---

## 1. Architecture Overview

```
MicrobialExplorer/
├── R/                          ← 13 files (10 modules + 3 support)
│   ├── utils.R                 ← Shared theme, colours, helpers, tidytree suppression
│   ├── validators.R            ← Input validation functions
│   ├── precompute.R            ← Cache management
│   ├── microbiome_qc.R         ← Module 1: Quality control
│   ├── microbiome_composition.R ← Module 2: Taxonomic composition
│   ├── microbiome_alpha_diversity.R ← Module 3: Alpha diversity
│   ├── microbiome_beta_diversity.R  ← Module 4: Beta diversity
│   ├── microbiome_differential_abundance.R ← Module 5: DA
│   ├── microbiome_functional.R ← Module 6: Functional prediction
│   ├── microbiome_network.R    ← Module 7: Co-occurrence network
│   ├── microbiome_ml_classification.R ← Module 8: ML
│   ├── microbiome_correlation.R ← Module 9: Correlation
│   └── microbiome_longitudinal.R ← Module 10: Longitudinal
├── app/                        ← Shiny app (TO BUILD)
│   ├── global.R
│   ├── ui.R
│   ├── server.R
│   └── modules/
│       ├── mod_upload.R
│       ├── mod_qc.R
│       ├── mod_composition.R
│       ├── mod_alpha.R
│       ├── mod_beta.R
│       ├── mod_da.R
│       ├── mod_functional.R
│       ├── mod_network.R
│       ├── mod_ml.R
│       ├── mod_correlation.R
│       └── mod_longitudinal.R
├── tests/
│   ├── test_all_modules.R
│   └── simulate_microbiome_data.R
├── config/pipeline_config.R
├── cache/                      ← 10 subdirectories
├── data/                       ← raw, example, reference
├── output/                     ← plots, tables
├── docs/
└── logs/
```

**Design principle:** Each Shiny module wraps exactly one R pipeline module. The precompute cache runs once on upload, then every tab loads from cache for fast interaction.

**Project folder:** `/Users/pagudeloromero/Library/CloudStorage/OneDrive-TheKidsResearchInstituteAustralia/MicrobiomeWebSite/MicrobialExplorer`

---

## 2. Data Type Detection (Critical for Upload Module)

### Lesson from testing:
The enterotype dataset (relative abundance, 0–1 range) correctly triggered validation errors in count-dependent functions. The upload module MUST detect data type automatically.

### Detection logic:
```r
detect_data_type <- function(otu_mat) {
  vals <- as.numeric(otu_mat[otu_mat > 0])
  if (all(vals == floor(vals))) return("counts")
  if (max(vals) <= 1 && all(rowSums(otu_mat) <= 1.01)) return("relative")
  return("normalised")
}
```

### Additional detections needed:
```r
detect_features <- function(ps) {
  list(
    data_type     = detect_data_type(as.matrix(otu_table(ps))),
    has_taxonomy  = !is.null(tax_table(ps, errorIfNULL = FALSE)),
    has_tree      = !is.null(phy_tree(ps, errorIfNULL = FALSE)),
    avail_ranks   = tryCatch(rank_names(ps), error = function(e) character(0)),
    n_samples     = nsamples(ps),
    n_taxa        = ntaxa(ps),
    sample_vars   = sample_variables(ps),
    numeric_vars  = names(which(sapply(data.frame(sample_data(ps)), is.numeric)))
  )
}
```

### Module availability by data type:

| Module             | Counts | Relative | Normalised | Needs taxonomy | Needs tree |
|--------------------|--------|----------|------------|----------------|------------|
| QC depth           | ✓      | ✓ (warn) | ✓ (warn)   | No             | No         |
| QC rarefaction     | ✓      | ✗        | ✗          | No             | No         |
| QC filter taxa     | ✓      | ✓*       | ✓*         | Optional       | No         |
| Composition        | ✓      | ✓        | ✓          | Yes            | No         |
| Alpha diversity    | ✓      | ✗        | ✗          | No             | Optional (Faith's PD) |
| Beta diversity     | ✓      | ✓ (no rarefy) | ✓     | No             | No         |
| DESeq2             | ✓      | ✗        | ✗          | Optional       | No         |
| ALDEx2             | ✓      | ✗        | ✗          | Optional       | No         |
| Network            | ✓      | ✓        | ✓          | Optional       | No         |
| ML classification  | ✓      | ✓        | ✓          | Optional       | No         |
| Correlation        | ✓      | ✓        | ✓          | Optional       | No         |
| Longitudinal       | ✓      | ✓*       | ✓*         | Optional       | No         |

*With adjusted thresholds

---

## 3. Function Signatures Reference (CONFIRMED BY TESTING)

### Module 1 — QC:
```r
qc_sequencing_depth(ps)
qc_rarefaction_curves(ps, n_samples = 20)
qc_filter_taxa(ps, min_prevalence = 0.10, min_abundance = 5)
  # Returns: list with $ps_filtered, $plot, $summary
qc_filter_samples(ps, min_reads = 100)
  # Returns: list with $ps_filtered
```

### Module 2 — Composition:
```r
agglomerate_taxa(ps, rank = "Genus", transform = "relative", top_n = 20)
  # Returns: phyloseq (directly, not a list)
plot_composition_bars(ps, rank = "Phylum", group_var = NULL, facet_var = NULL,
                       sort_by = "group", show_legend = TRUE)
  # NO top_n parameter — handled inside agglomerate_taxa
plot_mean_composition(ps, rank = "Phylum", group_var = NULL)
calculate_fb_ratio(ps, group_var = NULL)
```

### Module 3 — Alpha Diversity:
```r
calculate_alpha_diversity(ps, rarefaction = TRUE, rare_depth = NULL, n_rare = 3)
  # Rejects relative abundance data — requires counts
  # Faith's PD computed when tree present (picante installed)
plot_alpha_diversity(diversity_df, metrics = c("observed","shannon","pielou"),
                      group_var = "disease_status")
plot_diversity_gradient(diversity_df, metric = "shannon", x_var = "age",
                         group_var = NULL, add_smooth = TRUE, smooth_method = "loess")
  # Parameter is x_var, NOT variable
summarise_alpha_stats(diversity_df, group_var = "disease_status")
```

### Module 4 — Beta Diversity:
```r
compute_distances(ps, methods = c("bray","jaccard"), rare_depth = NULL)
  # Returns: named list of dist objects ($bray, $jaccard)
run_ordination(dist_obj, method = "PCoA", k = 3)
  # Takes dist_obj only — NO ps parameter
plot_ordination(ordination_result, group_var = "disease_status")
run_permanova(dist_obj, ps, formula_rhs = "disease_status",
               permutations = 999, strata = NULL)
  # Parameter is formula_rhs, NOT group_var
run_mantel_test(dist1, dist2)
```

### Module 5 — Differential Abundance:
```r
prepare_da_data(ps, rank = "Genus", min_prevalence = 0.10, min_count = 10)
  # Returns: phyloseq (directly, not a list)
  # NO group_var parameter
run_deseq2(ps, group_var = "group", reference = NULL, alpha = 0.05,
            lfc_threshold = 1)
plot_volcano(deseq_result)
plot_effect_sizes(result_df, lfc_threshold = 1, ...)
  # lfc_threshold parameter was added during testing
build_consensus(...)
```

### Module 7 — Network:
```r
prepare_network_data(ps, rank = "Genus", min_prevalence = 0.20, max_taxa = 40)
build_network(ps_net, cor_threshold = 0.3, p_threshold = 0.05)
  # Returns: list with $graph, $cor_mat, $edge_df
calculate_network_topology(graph, cor_mat, ps_net)
  # Returns: list with $global, $nodes (data.frame), $community
  # NOTE: return key is $nodes NOT $node_df
analyse_hub_taxa(graph, node_df, top_n = 20)
analyse_modules(graph, node_df, ps = NULL, rank = "Phylum")
compare_random_network(graph, ...)
```

### Module 8 — ML:
```r
prepare_ml_data(ps, group_var = "disease_status")
  # NO top_n parameter
train_random_forest(ml_data)
train_lasso(ml_data)
plot_feature_importance(...)
plot_roc_curves(...)
```

### Module 9 — Correlation:
```r
correlate_taxa_metadata(ps, rank = "Genus", meta_vars = NULL, transform = "clr",
                         method = "spearman", p_adjust = "BH", alpha = 0.05)
  # Parameter is meta_vars, NOT metadata_vars
compute_taxa_taxa_correlation(ps, rank = "Genus", top_n = 25)
```

---

## 4. Parameter Name Mapping for Shiny

The UI will have a single "Group variable" selector. The server maps this to the correct parameter name:

```r
# In server.R
gv <- input$group_variable

# Module-specific mapping:
run_permanova(dist, ps, formula_rhs = gv)       # formula_rhs
run_deseq2(ps, group_var = gv)                   # group_var
plot_alpha_diversity(div_df, group_var = gv)      # group_var
plot_composition_bars(ps, group_var = gv)         # group_var
prepare_ml_data(ps, group_var = gv)              # group_var
correlate_taxa_metadata(ps, meta_vars = numeric_vars)  # meta_vars (different)
```

---

## 5. Known Edge Cases (ALL DISCOVERED DURING REAL DATA TESTING)

### A. sample_data S4 class issue:
`sample_data(ps)` returns S4 object. `as.data.frame()` still returns `sample_data` class.

**Fix applied in modules:**
```r
# Standard approach (works in most cases):
data.frame(sample_data(ps))

# For stubborn cases (e.g., KEGG metadata):
if (inherits(metadata, "sample_data")) {
  rn <- rownames(metadata)
  cn <- colnames(metadata)
  metadata <- data.frame(
    lapply(setNames(cn, cn), function(col) metadata[[col]]),
    row.names = rn, stringsAsFactors = FALSE
  )
}
```

### B. otu_table class issue:
`t(as.matrix(otu_table(ps)))` returns `otu_table`, not `matrix`. vegan functions fail.

**Fix applied:**
```r
otu_mat <- as.matrix(otu_table(ps))
otu_t <- matrix(t(otu_mat), nrow = ncol(otu_mat), ncol = nrow(otu_mat),
                dimnames = list(colnames(otu_mat), rownames(otu_mat)))
```

### C. tax_table validation error:
`as.data.frame(tax_table(ps))` can trigger `invalid class "taxonomyTable"` error.

**Fix applied (3 locations in network module):**
```r
tax_df <- tryCatch({
  raw_mat <- ps@tax_table@.Data
  df <- data.frame(raw_mat, stringsAsFactors = FALSE)
  rownames(df) <- taxa_names(ps)
  colnames(df) <- colnames(ps@tax_table)
  df %>% rownames_to_column("taxon") %>%
    select(taxon, any_of(c("Phylum", "Class", "Family", rank)))
}, error = function(e) NULL)
if (!is.null(tax_df)) node_df <- left_join(node_df, tax_df, by = "taxon")
```

### D. Genus name duplication:
Real data has duplicated genus names (Clostridium, Bacteroides from different families) and NAs (226/728 in GlobalPatterns).

**Fix applied in agglomerate_taxa:**
```r
new_names <- as.character(tax_table(ps_agg)[, rank])
new_names[is.na(new_names) | new_names == "NA"] <- paste0("Unknown_", rank, "_", seq_len(sum(...)))
new_names <- make.unique(new_names, sep = "_dup")
taxa_names(ps_agg) <- new_names
```

### E. phylo/tidytree class conflict:
`Found more than one class "phylo" in cache` warning — very noisy, harmless.

**Fix applied in utils.R:**
```r
suppressMessages(suppressWarnings({
  if (requireNamespace("tidytree", quietly = TRUE)) library(tidytree)
}))
```

### F. igraph list columns:
Network node_df can have list-type columns from igraph. Breaks `fct_reorder()`.

**Fix applied in analyse_hub_taxa and analyse_modules:**
```r
node_df <- node_df %>% mutate(across(where(is.list), ~as.vector(unlist(.))))
```

### G. Datasets without taxonomy:
soilrep has no tax_table. Composition and agglomeration modules fail.

**Shiny must:** Check `has_taxonomy` and disable taxonomy-dependent features.

### H. Relative abundance data:
enterotype is pre-normalized (0–1). Alpha diversity and DESeq2 correctly reject it.

**Shiny must:** Check `data_type` and disable count-dependent modules.

---

## 6. Package Availability

### Required (all installed and working):
ggplot2, dplyr, tidyr, patchwork, vegan, phyloseq, microbiome, decontam,
DESeq2, ALDEx2, scales, tibble, RColorBrewer, stringr, forcats, ggrepel,
zoo, igraph, ggraph, tidygraph, randomForest, caret, pROC, e1071, glmnet,
rstatix, ggpubr, corrplot, Hmisc, lme4, lmerTest, picante, ape

### Unavailable (IT restrictions — no Homebrew):
- **ANCOMBC** — blocked by CVXR → Rmpfr → mpfr.h (system C library)
- **SpiecEasi** — optional (network uses spearman correlation instead)

### Shiny startup should:
```r
# In global.R
ANCOMBC_AVAILABLE <- requireNamespace("ANCOMBC", quietly = TRUE)
PICANTE_AVAILABLE <- requireNamespace("picante", quietly = TRUE)

# Disable UI elements for unavailable methods
# Show message: "ANCOMBC not available — using DESeq2 + ALDEx2 consensus"
```

---

## 7. Caching Strategy

### Pre-compute on upload (runs once, ~30–60 seconds):
```r
precompute_cache <- function(ps, group_var, output_dir) {
  cache <- list()
  cache$data_type   <- detect_data_type(as.matrix(otu_table(ps)))
  cache$features    <- detect_features(ps)
  cache$ps_clean    <- qc_filter_taxa(ps)$ps_filtered
  if (cache$data_type == "counts") {
    cache$distances <- compute_distances(cache$ps_clean)
    cache$alpha     <- calculate_alpha_diversity(cache$ps_clean)
  }
  cache$ordination  <- run_ordination(cache$distances$bray)
  saveRDS(cache, file.path(output_dir, "precomputed_cache.rds"))
}
```

### Load from cache for each tab (< 2 seconds):
Each Shiny module reads from the cache rather than re-computing.

---

## 8. Server Specifications (BioCommons)

- **Contact:** Ziad Al Bkhetan
- **Specs:** 500GB storage, 32GB RAM, 8 CPUs (scalable for V2)
- **Architecture:** Pre-computed cache on upload, load from cache for interactions
- **No user accounts in V1** (fresh session each time)
- **Domain:** microbialexplorer.org (to register ~$15/year)

---

## 9. Shiny UI Layout

### Sidebar:
- Upload panel (OTU table, taxonomy, metadata, tree — all optional except OTU)
- Data summary card (samples, taxa, data type, tree present)
- Group variable selector (auto-populated from metadata columns)
- Taxonomic rank selector (auto-populated from taxonomy ranks, disabled if no taxonomy)
- Global parameters: prevalence threshold, abundance threshold

### Main panel — tabs:
1. **Dashboard** — Summary stats, data overview, QC status
2. **QC** — Sequencing depth, rarefaction curves, filtering controls
3. **Composition** — Stacked bars, mean composition, heatmap, F:B ratio
4. **Alpha Diversity** — Box plots, gradient, correlation matrix, stats table
5. **Beta Diversity** — PCoA/NMDS plot, PERMANOVA table, distance heatmap
6. **Differential Abundance** — Volcano, effect sizes, DA heatmap, consensus
7. **Functional** — Pathway abundance, KEGG categories, functional diversity
8. **Network** — Interactive network visualisation, topology, hub taxa
9. **ML** — ROC curves, feature importance, confusion matrix
10. **Correlation** — Taxa-metadata heatmap, bubble plot, scatter plots
11. **Longitudinal** — Trajectories, stability, composition over time
12. **Export** — Download all plots, tables, and reports

### Tab enable/disable logic:
```r
# Disable tabs based on data features
if (data_type != "counts") disable("Alpha Diversity", "Differential Abundance")
if (!has_taxonomy) disable("Composition")
if (!has_longitudinal) disable("Longitudinal")
if (n_groups < 2) disable("Differential Abundance", "ML")
if (n_groups != 2) disable("ML")  # Binary classification only in V1
```

---

## 10. V1 vs V2 Scope

### V1 (current — bacterial microbiome):
- All 10 analysis modules
- Shiny web interface on BioCommons
- Pre-computed caching
- Proof of concept for BioCommons grant

### V2 (future — virome expansion):
- vOTU profiling module
- Viral taxonomy (ICTV integration)
- Phage lifestyle prediction (temperate vs lytic)
- AMG detection module
- Host-virus interaction network
- Integration with viral databases (INPHARED, RefSeq Viral)

---

## 11. Testing Summary — 5 Datasets Validated

| Dataset        | Samples | Taxa   | Type     | Tree | Taxonomy | Modules tested | Key findings |
|----------------|---------|--------|----------|------|----------|----------------|-------------|
| Simulated      | 60      | 200    | Counts   | ✓    | K→S      | All 10         | 69/69 (100%) |
| GlobalPatterns | 9       | 5,281  | Counts   | ✓    | K→S      | QC, Comp, α, β, DA, Net, Corr | PERMANOVA p=0.002, F:B=2.21, 6 DA phyla, Faith's PD with tree |
| Enterotype     | 271     | 553    | Relative | ✗    | Genus only | Corr only  | Correctly rejects count-dependent modules, 243 sig pairs |
| soilrep        | 56      | 16,825 | Counts   | ✗    | None     | QC, α, β    | Handles no-taxonomy gracefully, 82.8% taxa filtered |
| RISK_CCFA      | 1,051   | 1,595  | Counts   | ✓    | K→S      | QC, Comp, α, β, DA, Net, Corr | Faith's PD p=0.005 (CD), PERMANOVA p=0.001, 70 DA taxa, F:B=0.81, 36-node network, hub: Blautia, Veillonella |

### Biological validation highlights:
- Body sites significantly different — PERMANOVA p=0.002 (GlobalPatterns)
- CD vs control separated — PERMANOVA p=0.001 (RISK_CCFA)
- Faith's PD reduced in Crohn's disease — p=0.005 (known biological finding)
- F:B ratio 0.81 in CD (reduced, expected) vs 2.21 in feces (healthy, expected)
- Hub genera in IBD network: Blautia, Veillonella, Parabacteroides, Ruminococcus
- 70 differentially abundant taxa in CD vs control via DESeq2
- Soil warming shows no significant effect — p=0.498 (consistent with literature)

### Edge cases validated:
- ✓ No taxonomy table (soilrep)
- ✓ Relative abundance input (enterotype)
- ✓ Single taxonomy rank (enterotype — Genus only)
- ✓ Duplicated genus names (GlobalPatterns — Clostridium, Bacteroides, etc.)
- ✓ High sparsity / 16K+ taxa (soilrep)
- ✓ Small sample size n=9 (GlobalPatterns)
- ✓ Large sample size n=1,051 (RISK_CCFA)
- ✓ With phylogenetic tree (GlobalPatterns, RISK_CCFA)
- ✓ Without phylogenetic tree (enterotype, soilrep)
- ✓ picante working for Faith's PD

---

*Last updated: May 2026*
*Project: MicrobialExplorer — microbialexplorer.org*
*Author: Patricia Agudelo-Romero*
*Status: Pipeline validated — Ready for Shiny development*
