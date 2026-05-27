# MicrobialExplorer — Shiny App

## Structure

```
app/
├── global.R          ← Package loading, source R modules, constants, detect_features()
├── ui.R              ← bslib page_navbar layout, 12 tabs, sidebar
├── server.R          ← Reactive data pipeline, tab enable/disable, module dispatch
├── modules/
│   ├── mod_upload.R  ← COMPLETE — file upload, demo loader, validation, precompute
│   ├── mod_qc.R      ← next to build
│   ├── mod_composition.R
│   ├── mod_alpha.R
│   ├── mod_beta.R
│   ├── mod_da.R
│   ├── mod_functional.R
│   ├── mod_network.R
│   ├── mod_ml.R
│   ├── mod_correlation.R
│   └── mod_longitudinal.R
└── www/
    └── custom.css    ← App styles
```

## Running the app

```r
# From the project root:
shiny::runApp("app")

# Or from within the app/ directory:
shiny::runApp()
```

## Dependencies (beyond the pipeline R packages)

```r
install.packages(c(
  "shiny",
  "shinyjs",
  "bslib",         # Bootstrap 5 layout
  "DT",            # Interactive tables
  "bsplus"         # bsCollapse for demo panel (optional, can replace)
))
```

## Module build order (recommended)

Each module follows the same pattern:
1. `mod_X_ui(id)` — sidebar controls + main panel layout
2. `mod_X_server(id, upload_data, input)` — reads from `upload_data$ps` and `upload_data$cache`
3. Wire into `server.R` by replacing the stub `renderUI` with `mod_X_server(...)`

Priority order:
1. `mod_qc.R` — most users will run QC first
2. `mod_composition.R`
3. `mod_alpha.R`
4. `mod_beta.R`
5. `mod_da.R`
6. `mod_network.R`
7. `mod_correlation.R`
8. `mod_ml.R`
9. `mod_functional.R`
10. `mod_longitudinal.R`

## Key design decisions

- **Precompute on upload**: `precompute_cache()` in `mod_upload.R` runs QC filtering,
  alpha diversity, and distance matrices once. Each module tab reads from `upload_data$cache`
  for near-instant response.

- **Feature detection**: `detect_features()` in `global.R` inspects the phyloseq object
  and returns flags (`data_type`, `has_taxonomy`, `has_tree`, etc.) used to enable/disable
  tabs and populate selectors.

- **Parameter name mapping**: The sidebar has one `group_variable` input. `server.R`
  passes this as `formula_rhs` to PERMANOVA and `group_var` to DESeq2/ML (see design notes §4).

- **ANCOMBC fallback**: Detected at startup. DA module uses DESeq2 + ALDEx2 consensus
  when ANCOMBC is unavailable.
