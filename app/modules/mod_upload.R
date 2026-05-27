# =============================================================================
# mod_upload.R — Upload & Data Import Module
# =============================================================================
# Handles:
#   1. File upload (OTU table, taxonomy, metadata, phylogenetic tree)
#   2. Demo dataset loading (GlobalPatterns, enterotype, soilrep)
#   3. phyloseq object construction with validation
#   4. Data type and feature detection
#   5. Precompute cache trigger
#
# Returns reactive values: ps, features, cache
# =============================================================================

# ── UI: sidebar upload panel ──────────────────────────────────────────────────
mod_upload_ui <- function(id) {
  ns <- NS(id)

  tagList(
    tags$h6(
      class = "fw-bold mb-2 text-muted text-uppercase small",
      icon("upload"), " Data input"
    ),

    # ── Demo dataset shortcut ───────────────────────────────────────────────
    actionButton(
      ns("load_demo"),
      label = tagList(icon("flask"), " Load demo data"),
      class = "btn btn-outline-secondary btn-sm w-100 mb-2"
    ),

    # Demo dataset options — simple dropdown, no bsplus needed
    selectInput(
      ns("demo_choice"),
      label = NULL,
      choices = c(
        "GlobalPatterns (body sites, n=9)"  = "GlobalPatterns",
        "Enterotype (gut microbiome, n=271)" = "enterotype",
        "Soilrep (soil warming, n=56)"       = "soilrep",
        "Network enriched (IBD/Healthy, n=60)" = "network_enriched"
      ),
      selected = "network_enriched"
    ),
    tags$p(class = "text-center text-muted small my-1", "— or upload your own —"),

    # ── OTU / ASV table ────────────────────────────────────────────────────
    fileInput(
      ns("otu_file"),
      label       = tooltip(
        tags$span("OTU / ASV table ", tags$span(class = "text-danger", "*"),
                  icon("circle-question", class = "small")),
        "Required. Rows = taxa, columns = samples. CSV or TSV with row names."
      ),
      accept      = ACCEPTED_OTU,
      placeholder = "taxa × samples (.csv/.tsv)"
    ),

    # ── Taxonomy table ─────────────────────────────────────────────────────
    fileInput(
      ns("tax_file"),
      label  = tooltip(
        tags$span("Taxonomy table ", icon("circle-question", class = "small")),
        "Optional. Rows = taxa (must match OTU table), columns = ranks (Kingdom … Species)."
      ),
      accept = ACCEPTED_TAX,
      placeholder = "taxa × ranks (.csv/.tsv)"
    ),

    # ── Sample metadata ────────────────────────────────────────────────────
    fileInput(
      ns("meta_file"),
      label  = tooltip(
        tags$span("Sample metadata ", icon("circle-question", class = "small")),
        "Optional. Rows = samples (must match OTU table columns), columns = variables."
      ),
      accept = ACCEPTED_META,
      placeholder = "samples × variables (.csv/.tsv)"
    ),

    # ── Phylogenetic tree ──────────────────────────────────────────────────
    fileInput(
      ns("tree_file"),
      label  = tooltip(
        tags$span("Phylogenetic tree ", icon("circle-question", class = "small")),
        "Optional. Newick or Nexus format. Enables Faith's PD in alpha diversity."
      ),
      accept = ACCEPTED_TREE,
      placeholder = "Newick / Nexus"
    ),

    # ── Build / import button ──────────────────────────────────────────────
    actionButton(
      ns("build_ps"),
      label = tagList(icon("play"), " Import data"),
      class = "btn btn-primary btn-sm w-100 mt-1",
      disabled = "disabled"
    ),

    # ── Validation feedback ────────────────────────────────────────────────
    uiOutput(ns("validation_output"))
  )
}


# ── UI: home / landing panel ──────────────────────────────────────────────────
mod_upload_landing_ui <- function(id) {
  ns <- NS(id)

  tagList(
    layout_columns(
      col_widths = c(8, 4),

      # Welcome card
      card(
        class = "border-0",
        card_body(
          tags$h2(
            class = "fw-bold mb-1",
            icon("bacteria", class = "text-success"), " MicrobialExplorer"
          ),
          tags$p(
            class = "lead text-muted",
            "An interactive platform for microbiome data analysis. ",
            "Upload your data or load a demo dataset to get started."
          ),
          hr(),
          tags$h5("What you can do"),
          tags$ul(
            tags$li("Quality control — depth, rarefaction, prevalence filtering"),
            tags$li("Taxonomic composition — stacked bars, heatmaps, F:B ratio"),
            tags$li("Alpha diversity — observed, Shannon, Faith's PD"),
            tags$li("Beta diversity — PCoA/NMDS, PERMANOVA, distance heatmaps"),
            tags$li("Differential abundance — DESeq2 + ALDEx2 consensus"),
            tags$li("Functional prediction — PICRUSt2 pathway visualisation"),
            tags$li("Co-occurrence networks — Spearman + hub taxon analysis"),
            tags$li("ML classification — Random Forest + LASSO, ROC curves"),
            tags$li("Correlation — taxa–metadata and taxa–taxa heatmaps"),
            tags$li("Longitudinal — trajectory and stability analysis")
          ),
          hr(),
          tags$p(
            class = "small text-muted",
            icon("circle-info"), " Upload your files using the panel on the left, ",
            "or click ", tags$strong("Load demo data"), " to explore with an example dataset."
          )
        )
      ),

      # Quick-start card
      card(
        card_header(icon("rocket"), " Quick start"),
        card_body(
          tags$ol(
            tags$li("Upload your OTU table (required)"),
            tags$li("Optionally add taxonomy, metadata, and tree"),
            tags$li('Click "Import data"'),
            tags$li("Select a group variable and rank"),
            tags$li("Navigate to any analysis tab")
          ),
          hr(),
          tags$h6("File formats"),
          tags$table(
            class = "table table-sm table-bordered small",
            tags$thead(tags$tr(tags$th("File"), tags$th("Format"))),
            tags$tbody(
              tags$tr(tags$td("OTU table"), tags$td("CSV / TSV, taxa × samples")),
              tags$tr(tags$td("Taxonomy"),  tags$td("CSV / TSV, taxa × ranks")),
              tags$tr(tags$td("Metadata"),  tags$td("CSV / TSV, samples × vars")),
              tags$tr(tags$td("Tree"),      tags$td("Newick / Nexus"))
            )
          ),
          hr(),
          tags$h6("Demo datasets"),
          tags$p(class = "small text-muted",
                 "GlobalPatterns · Enterotype · Soilrep"),
          tags$p(class = "small text-muted",
                 "(all from the phyloseq R package)")
        )
      )
    ),

    # Status card — shown after data is loaded
    uiOutput(ns("landing_status"))
  )
}


# ── Server ────────────────────────────────────────────────────────────────────
  mod_upload_server <- function(id) {
    moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Return values shared with the parent server
    rv <- reactiveValues(
      ps       = NULL,   # phyloseq object
      features = NULL,   # output of detect_features()
      cache    = NULL    # output of precompute_cache() [once implemented]
    )

    # ── Enable Import button once OTU file is chosen (or demo selected) ───
    observe({
      if (!is.null(input$otu_file) || input$load_demo > 0)
        shinyjs::enable("build_ps")
    })

    # ── Demo dataset loader ───────────────────────────────────────────────
    observeEvent(input$load_demo, {
      tryCatch({
        shinyjs::disable("build_ps")
        .show_status("Loading demo dataset…", type = "info")
        Sys.sleep(0.3)   # allow UI to update

        ps_demo <- load_demo_data(input$demo_choice)
        .finalise_import(ps_demo, source = input$demo_choice)
      }, error = function(e) {
        .show_status(paste("Error loading demo:", conditionMessage(e)), type = "danger")
      })
    })

    # ── Import button: build phyloseq from uploaded files ─────────────────
    observeEvent(input$build_ps, {

      req(input$otu_file)
      shinyjs::disable("build_ps")
      .show_status("Reading files…", type = "info")

      tryCatch({

        # ── 1. Read OTU table ────────────────────────────────────────────
        otu_raw <- read_table_file(input$otu_file$datapath)
        
        otu_mat <- as.matrix(otu_raw)
        storage.mode(otu_mat) <- "numeric"

        v_otu <- validate_otu_table(otu_mat)
        if (!v_otu$valid) {
          .show_validation_errors("OTU table", v_otu)
          shinyjs::enable("build_ps")
          return()
        }

        .show_status("Building phyloseq object…", type = "info")

        # ── 2. Read taxonomy (optional) ───────────────────────────────────
        TAX <- NULL
        if (!is.null(input$tax_file)) {
          tax_raw <- tryCatch({
            tmp <- read.table(input$tax_file$datapath, sep="\t", header=TRUE,
                              check.names=FALSE, comment.char="", quote="",
                              stringsAsFactors=FALSE)
            # Find ID column - first non-taxonomy column
            id_col <- which(grepl("^(Feature|feature|id|ID|#OTU)", colnames(tmp)))[1]
            if (is.na(id_col)) id_col <- 1
            rownames(tmp) <- as.character(tmp[[id_col]])
            tmp <- tmp[, -id_col, drop=FALSE]
            tmp
          }, error = function(e) read_table_file(input$tax_file$datapath))




          
          v_tax   <- validate_taxonomy_table(as.matrix(tax_raw))
          
          if (!v_tax$valid) {
            .show_validation_errors("Taxonomy table", v_tax)
            shinyjs::enable("build_ps")
            return()
          }
          # Align rows to OTU table
          shared_taxa <- intersect(rownames(otu_mat), rownames(tax_raw))
          if (length(shared_taxa) == 0) {
            .show_status(
              "No shared taxa between OTU table and taxonomy table. Check row names.",
              type = "danger"
            )
            shinyjs::enable("build_ps")
            return()
          }
          tax_mat <- as.matrix(tax_raw[shared_taxa, , drop = FALSE])
          otu_mat <- otu_mat[shared_taxa, , drop = FALSE]
          TAX     <- tax_table(tax_mat)
        }

        # ── 3. Read metadata (optional) ───────────────────────────────────
        SAM <- NULL
        if (!is.null(input$meta_file)) {
          meta_raw <- read_table_file(input$meta_file$datapath)
          
          v_meta   <- validate_metadata(meta_raw)
          
          if (!v_meta$valid) {
            .show_validation_errors("Metadata", v_meta)
            shinyjs::enable("build_ps")
            return()
          }
          # Align samples
          otu_samples  <- trimws(colnames(otu_mat))
          meta_samples <- trimws(rownames(meta_raw))
          
          colnames(otu_mat)    <- otu_samples
          rownames(meta_raw)   <- meta_samples
          
          shared_samp <- intersect(otu_samples, meta_samples)
          
          if (length(shared_samp) == 0) {
            .show_status(
              paste0(
                "No shared samples between OTU table and metadata. ",
                "OTU columns: ", paste(head(otu_samples, 3), collapse = ", "), " … ",
                "Metadata rows: ", paste(head(meta_samples, 3), collapse = ", "), " …"
              ),
              type = "danger"
            )
            shinyjs::enable("build_ps")
            return()
          }
          
          otu_mat  <- otu_mat[, shared_samp, drop = FALSE]
          meta_raw <- meta_raw[shared_samp, , drop = FALSE]
          SAM      <- sample_data(meta_raw)
        }

        # ── 4. Read tree (optional) ───────────────────────────────────────
        PHY <- NULL
        if (!is.null(input$tree_file)) {
          PHY <- tryCatch(
            ape::read.tree(input$tree_file$datapath),
            error = function(e) {
              tryCatch(
                ape::read.nexus(input$tree_file$datapath),
                error = function(e2) NULL
              )
            }
          )
          if (is.null(PHY)) {
            .show_status(
              "Could not parse the tree file. Check it is valid Newick or Nexus format.",
              type = "warning"
            )
          }
        }

        # ── 5. Align optional tree to OTU taxa ────────────────────────────
        if (!is.null(PHY)) {
          shared_tips <- intersect(PHY$tip.label, rownames(otu_mat))
          
          if (length(shared_tips) == 0) {
            .show_status(
              "Tree tips do not match OTU taxa. Tree will be ignored.",
              type = "warning"
            )
            PHY <- NULL
          } else {
            PHY <- ape::keep.tip(PHY, shared_tips)
            
            # Reorder OTU/taxonomy to match tree tip order
            otu_mat <- otu_mat[PHY$tip.label, , drop = FALSE]
            
            if (!is.null(TAX)) {
              tax_mat <- tax_mat[PHY$tip.label, , drop = FALSE]
              TAX <- tax_table(tax_mat)
            }
          }
        }
        
        # ── 6. Build OTU table component ─────────────────────────────────
        OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
        
        # ── 7. Assemble phyloseq object ───────────────────────────────────
        components <- list(OTU)
        if (!is.null(TAX)) components <- c(components, list(TAX))
        if (!is.null(SAM)) components <- c(components, list(SAM))
        if (!is.null(PHY)) components <- c(components, list(phy_tree(PHY)))
        
        ps <- do.call(phyloseq, components)

        .finalise_import(ps, source = "user upload")

      }, error = function(e) {
        msg <- paste("Unexpected error during import:", conditionMessage(e))
        message("[mod_upload] ", msg)
        traceback()
        
        .show_status(msg, type = "danger")
        shinyjs::enable("build_ps")
      })
      
    })   # closes observeEvent(input$build_ps, ...)
    
    # ── Internal: finalise import ...
    .finalise_import <- function(ps, source = "unknown") {

      # Validate the assembled phyloseq object
      v_ps <- validate_phyloseq(ps)
      if (!v_ps$valid) {
        .show_validation_errors(paste0("phyloseq object (", source, ")"), v_ps)
        shinyjs::enable("build_ps")
        return()
      }

      # Detect features
      feat  <- detect_features(ps)
      rv$ps       <- ps
      rv$features <- feat

      # Trigger precompute cache in background (non-blocking via future if available)
      # For V1 this runs synchronously; replace with future::future() for V2
      cache <- tryCatch(
        precompute_cache(ps, feat),
        error = function(e) {
          message("[mod_upload] precompute_cache failed: ", conditionMessage(e))
          NULL
        }
      )
      rv$cache <- cache

      # Build success message
      warn_msgs <- v_ps$warnings
      if (length(warn_msgs) > 0) {
        .show_validation_warnings(
          paste0("Data loaded from ", source),
          warn_msgs
        )
      } else {
        .show_status(
          paste0(
            "Loaded ", feat$n_samples, " samples × ",
            format(feat$n_taxa, big.mark = ","), " taxa from ", source, "."
          ),
          type = "success"
        )
      }

      shinyjs::enable("build_ps")
    }

    # ── Validation feedback helpers ───────────────────────────────────────
    .show_status <- function(msg, type = "info") {
      output$validation_output <- renderUI({
        tags$div(
          class = paste0("alert alert-", type, " p-2 small mt-2"),
          msg
        )
      })
    }

    .show_validation_errors <- function(source_label, v_result) {
      output$validation_output <- renderUI({
        tags$div(
          class = "alert alert-danger p-2 small mt-2",
          tags$strong(icon("xmark"), " ", source_label, " errors:"),
          tags$ul(
            class = "mb-0 ps-3 mt-1",
            lapply(v_result$errors, tags$li)
          ),
          if (length(v_result$warnings) > 0)
            tagList(
              tags$hr(class = "my-1"),
              tags$strong("Warnings:"),
              tags$ul(
                class = "mb-0 ps-3",
                lapply(v_result$warnings, tags$li)
              )
            )
        )
      })
    }

    .show_validation_warnings <- function(header, warnings) {
      output$validation_output <- renderUI({
        tags$div(
          class = "alert alert-warning p-2 small mt-2",
          tags$strong(icon("triangle-exclamation"), " ", header),
          tags$ul(
            class = "mb-0 ps-3 mt-1",
            lapply(warnings, tags$li)
          )
        )
      })
    }

    # ── Landing page status card ──────────────────────────────────────────
    output$landing_status <- renderUI({
      ps   <- rv$ps
      feat <- rv$features
      req(ps)

      avail <- get_module_availability(feat)
      n_avail <- sum(avail)

      card(
        class = "border-success mt-3",
        card_header(
          class = "bg-success text-white",
          icon("circle-check"), " Data loaded successfully"
        ),
        card_body(
          layout_columns(
            col_widths = c(4, 4, 4),
            tags$div(
              tags$strong("Samples"), tags$br(),
              tags$span(class = "fs-4 fw-bold text-success", feat$n_samples)
            ),
            tags$div(
              tags$strong("Taxa"), tags$br(),
              tags$span(class = "fs-4 fw-bold text-success",
                        format(feat$n_taxa, big.mark = ","))
            ),
            tags$div(
              tags$strong("Modules available"), tags$br(),
              tags$span(class = "fs-4 fw-bold text-success",
                        paste0(n_avail, " / 10"))
            )
          ),
          tags$p(
            class = "text-muted small mt-2 mb-0",
            "Navigate using the tabs above to begin your analysis."
          )
        )
      )
    })

    # ── Return reactive values to parent server ───────────────────────────
    return(rv)
  })
}


# =============================================================================
# precompute_cache() — runs once on upload, stores results for fast tab loading
# =============================================================================
# This is a thin wrapper here; the heavy lifting is in precompute.R.
# It gracefully skips steps that are not applicable to the data type.
# =============================================================================

precompute_cache <- function(ps, features) {

  cache <- list()
  cache$features  <- features
  cache$data_type <- features$data_type

  # ── QC filtering (always run) ─────────────────────────────────────────────
  cache$ps_clean <- tryCatch({
    qc_filter_taxa(ps)$ps_filtered
  }, error = function(e) ps)

  ps_c <- cache$ps_clean

  # ── Count-only precomputes ─────────────────────────────────────────────────
  if (features$data_type == "counts") {

    cache$alpha <- tryCatch(
      calculate_alpha_diversity(ps_c, rarefaction = TRUE),
      error = function(e) NULL
    )

    cache$distances <- tryCatch(
      compute_distances(ps_c, methods = c("bray", "jaccard")),
      error = function(e) NULL
    )

  } else {
    # Beta diversity works on relative/normalised data too (no rarefaction)
    cache$distances <- tryCatch(
      compute_distances(ps_c, methods = c("bray", "jaccard"), rare_depth = NULL),
      error = function(e) NULL
    )
  }

  # ── Ordination (if distances available) ───────────────────────────────────
  if (!is.null(cache$distances$bray)) {
    cache$ordination <- tryCatch(
      run_ordination(cache$distances$bray, method = "PCoA"),
      error = function(e) NULL
    )
  }

  cache
}

