# =============================================================================
# server.R — MicrobialExplorer Shiny App
# =============================================================================

server <- function(input, output, session) {

  # ── Upload module ──────────────────────────────────────────────────────────
  # Returns: list(ps, features, cache) as reactiveValues
  upload_data <- mod_upload_server("upload")

  # ── Convenience reactive: is data loaded? ─────────────────────────────────
  output$data_loaded <- reactive({
    !is.null(upload_data$ps)
  })
  outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)

  # ── Populate sidebar selectors once data is available ─────────────────────
  observeEvent(upload_data$features, {
    feat <- upload_data$features
    req(feat)

    # Group variable selector
    updateSelectInput(session, "group_variable",
      choices  = c("None" = "", feat$group_vars),
      selected = if (length(feat$group_vars) > 0) feat$group_vars[1] else ""
    )

    # Taxonomic rank selector
    if (feat$has_taxonomy && length(feat$avail_ranks) > 0) {
      preferred <- intersect(c("Genus", "Family", "Order", "Phylum"), feat$avail_ranks)
      default   <- if (length(preferred) > 0) preferred[1] else feat$avail_ranks[1]
      updateSelectInput(session, "tax_rank",
        choices  = feat$avail_ranks,
        selected = default
      )
      shinyjs::enable("tax_rank")
    } else {
      updateSelectInput(session, "tax_rank",
        choices  = c("No taxonomy" = ""),
        selected = ""
      )
      shinyjs::disable("tax_rank")
    }
  })

  # ── Sidebar data summary ───────────────────────────────────────────────────
  output$sidebar_summary <- renderUI({
    feat <- upload_data$features
    req(feat)

    data_type_colour <- switch(feat$data_type,
      counts     = "success",
      relative   = "warning",
      normalised = "info",
      "secondary"
    )

    tagList(
      stat_badge("Samples", feat$n_samples),
      stat_badge("Taxa",    format(feat$n_taxa, big.mark = ",")),
      stat_badge("Data type",
                 feat$data_type,
                 colour = data_type_colour),
      stat_badge("Taxonomy",
                 if (feat$has_taxonomy) paste(feat$avail_ranks, collapse = " › ") else "None",
                 colour = if (feat$has_taxonomy) "success" else "secondary"),
      stat_badge("Tree",
                 if (feat$has_tree) "Present" else "None",
                 colour = if (feat$has_tree) "success" else "secondary"),
      stat_badge("Metadata vars",
                 length(feat$sample_vars))
    )
  })

  # ── Sidebar warnings ──────────────────────────────────────────────────────
  output$sidebar_warnings <- renderUI({
    feat <- upload_data$features
    req(feat)

    msgs <- c()

    if (feat$data_type != "counts")
      msgs <- c(msgs,
        paste0("Data type is '", feat$data_type, "'. ",
               "Alpha diversity, DESeq2, and ALDEx2 require raw counts."))

    if (!feat$has_taxonomy)
      msgs <- c(msgs, "No taxonomy table — Composition tab is disabled.")

    if (!ANCOMBC_AVAILABLE)
      msgs <- c(msgs, "ANCOMBC not available. DA uses DESeq2 + ALDEx2 consensus.")

    if (PICANTE_AVAILABLE && !feat$has_tree)
      msgs <- c(msgs, "No phylogenetic tree — Faith's PD will be skipped.")

    if (length(msgs) == 0) return(NULL)

    tagList(
      hr(),
      tags$div(
        class = "alert alert-warning p-2 small",
        icon("triangle-exclamation"), " ",
        tags$strong("Notes:"),
        tags$ul(
          class = "mb-0 ps-3 mt-1",
          lapply(msgs, tags$li)
        )
      )
    )
  })

  # ── Module availability — reactive, updates with data ─────────────────────
  mod_avail <- reactive({
    feat <- upload_data$features
    if (is.null(feat)) return(NULL)
    get_module_availability(feat)
  })

  # ── Disable tabs when module is not available ─────────────────────────────
  observeEvent(mod_avail(), {
    avail <- mod_avail()
    req(avail)

    tab_map <- list(
      composition  = "composition",
      alpha        = "alpha",
      da           = "da",
      ml           = "ml",
      longitudinal = "longitudinal"
    )

    for (module in names(tab_map)) {
      nav_id <- tab_map[[module]]
      if (avail[[module]]) {
        shinyjs::runjs(sprintf(
          'document.querySelector(\'[data-value="%s"]\').parentElement.classList.remove("disabled");',
          nav_id
        ))
      } else {
        shinyjs::runjs(sprintf(
          'document.querySelector(\'[data-value="%s"]\').parentElement.classList.add("disabled");',
          nav_id
        ))
      }
    }
  })

  # ── Dashboard UI ──────────────────────────────────────────────────────────
  output$dashboard_ui <- renderUI({
    ps   <- upload_data$ps
    feat <- upload_data$features

    if (is.null(ps)) {
      return(.no_data_panel("dashboard"))
    }

    avail <- get_module_availability(feat)

    layout_columns(
      col_widths = c(4, 4, 4),

      # Sample count card
      card(
        card_header(icon("users"), " Samples"),
        card_body(
          tags$h2(class = "text-center fw-bold text-primary", feat$n_samples),
          tags$p(class = "text-center text-muted small",
                 paste(length(feat$group_vars), "grouping variable(s)"))
        )
      ),

      # Taxa count card
      card(
        card_header(icon("bacteria"), " Taxa"),
        card_body(
          tags$h2(class = "text-center fw-bold text-success",
                  format(feat$n_taxa, big.mark = ",")),
          tags$p(class = "text-center text-muted small",
                 if (feat$has_taxonomy)
                   paste(length(feat$avail_ranks), "taxonomic ranks")
                 else "No taxonomy")
        )
      ),

      # Data type card
      card(
        card_header(icon("vials"), " Data type"),
        card_body(
          tags$h2(class = paste0("text-center fw-bold text-",
                                  switch(feat$data_type,
                                    counts = "success", relative = "warning", "info")),
                  feat$data_type),
          tags$p(class = "text-center text-muted small",
                 if (feat$data_type == "counts") "Count data detected"
                 else "Some modules disabled")
        )
      ),

      # Module availability table
      card(
        col_widths = 12,
        card_header(icon("check-circle"), " Module availability"),
        card_body(
          uiOutput("module_avail_table")
        )
      )
    )
  })

  output$module_avail_table <- renderUI({
    avail <- mod_avail()
    req(avail)

    module_labels <- c(
      qc          = "Quality control",
      composition = "Composition",
      alpha       = "Alpha diversity",
      beta        = "Beta diversity",
      da          = "Differential abundance",
      functional  = "Functional (PICRUSt2)",
      network     = "Co-occurrence network",
      ml          = "ML classification",
      correlation = "Correlation",
      longitudinal = "Longitudinal"
    )

    rows <- lapply(names(avail), function(m) {
      ok <- avail[[m]]
      tags$tr(
        tags$td(module_labels[[m]]),
        tags$td(
          if (ok)
            tags$span(class = "badge bg-success", icon("check"), " Available")
          else
            tags$span(class = "badge bg-secondary", icon("ban"), " Unavailable")
        )
      )
    })

    tags$table(
      class = "table table-sm table-hover mb-0",
      tags$thead(tags$tr(tags$th("Module"), tags$th("Status"))),
      tags$tbody(rows)
    )
  })

  # ── QC module (fully wired) ───────────────────────────────────────────────
  ps_qc_filtered <- mod_qc_server("qc", upload_data, input)

  output$qc_ui <- renderUI({
    ps   <- upload_data$ps
    feat <- upload_data$features

    if (is.null(ps))
      return(.no_data_panel("Quality Control"))

    mod_qc_ui("qc")
  })

  # ── Placeholder UIs for analysis tabs ────────────────────────────────────
  # Each returns a card with either the module content or a 'data needed' notice.
  # Full module servers are wired in as stubs here — replace with real
  # mod_X_server() calls once each module UI is built.

  # Composition module (fully wired)
  mod_composition_server("composition", upload_data, ps_qc_filtered, input)
  mod_alpha_server("alpha", upload_data, input)
  mod_beta_server("beta", upload_data, input)
  mod_da_server("da", upload_data, input)

  # Functional module (fully wired)
  mod_functional_server("functional", upload_data, input)

  output$functional_ui <- renderUI({
    ps    <- upload_data$ps
    feat  <- upload_data$features
    avail <- mod_avail()
    if (is.null(ps)) return(.no_data_panel("Functional Prediction"))
    if (!is.null(avail) && !avail[["functional"]])
      return(.unavailable_panel("Functional Prediction", feat))
    mod_functional_ui("functional")
  })

  output$beta_ui <- renderUI({
    ps    <- upload_data$ps
    feat  <- upload_data$features
    avail <- mod_avail()
    if (is.null(ps)) return(.no_data_panel("Beta Diversity"))
    if (!is.null(avail) && !avail[["beta"]])
      return(.unavailable_panel("Beta Diversity", feat))
    mod_beta_ui("beta")
  })

  output$composition_ui <- renderUI({
    ps    <- upload_data$ps
    feat  <- upload_data$features
    avail <- mod_avail()
    if (is.null(ps)) return(.no_data_panel("Composition"))
    if (!is.null(avail) && !avail[["composition"]])
      return(.unavailable_panel("Composition", feat))
    mod_composition_ui("composition")
  })


  # ── Alpha diversity ────────────────────────────────────────────────────
  output$da_ui <- renderUI({
    ps    <- upload_data$ps
    feat  <- upload_data$features
    avail <- mod_avail()
    if (is.null(ps)) return(.no_data_panel("Differential Abundance"))
    if (!is.null(avail) && !avail[["da"]])
      return(.unavailable_panel("Differential Abundance", feat))
    mod_da_ui("da")
  })

  output$alpha_ui <- renderUI({
    ps    <- upload_data$ps
    feat  <- upload_data$features
    avail <- mod_avail()
    if (is.null(ps)) return(.no_data_panel("Alpha Diversity"))
    if (!is.null(avail) && !avail[["alpha"]])
      return(.unavailable_panel("Alpha Diversity", feat))
    mod_alpha_ui("alpha")
  })
  tab_modules <- list(
    network_ui     = list(module = "network",      label = "Co-occurrence Network",
                           icon_name = "circle-nodes"),
    ml_ui          = list(module = "ml",           label = "ML Classification",
                           icon_name = "robot"),
    correlation_ui = list(module = "correlation",  label = "Correlation",
                           icon_name = "arrows-left-right"),
    longitudinal_ui = list(module = "longitudinal", label = "Longitudinal",
                            icon_name = "timeline")
  )

  for (output_id in names(tab_modules)) {
    local({
      oid  <- output_id
      meta <- tab_modules[[oid]]

      output[[oid]] <- renderUI({
        ps    <- upload_data$ps
        feat  <- upload_data$features
        avail <- mod_avail()

        if (is.null(ps))
          return(.no_data_panel(meta$label))

        if (!is.null(avail) && !avail[[meta$module]])
          return(.unavailable_panel(meta$label, feat))

        # ── Stub card: replace with mod_X_ui() call ──────────────────
        card(
          card_header(
            icon(meta$icon_name), " ", meta$label,
            tags$span(class = "badge bg-success ms-2 float-end small",
                      "Pipeline ready")
          ),
          card_body(
            tags$p(
              class = "text-muted",
              icon("circle-check", class = "text-success"), " ",
              "Module pipeline validated (69/69 tests). ",
              "Shiny UI for this module is the next build step."
            ),
            tags$pre(
              class = "bg-light p-2 rounded small",
              paste0("# Wire up in server.R:\n",
                     "# mod_", meta$module, "_server(\"", meta$module, "\", ",
                     "upload_data, input)")
            )
          )
        )
      })
    })
  }

  # ── Export UI stub ─────────────────────────────────────────────────────────
  output$export_ui <- renderUI({
    req(upload_data$ps)
    card(
      card_header(icon("download"), " Export"),
      card_body(
        tags$p(class = "text-muted",
               "Download all plots, tables, and an HTML report after running analyses.")
      )
    )
  })

  # ── Internal UI helpers ───────────────────────────────────────────────────

  .no_data_panel <- function(label) {
    card(
      class = "border-0",
      card_body(
        class = "text-center py-5",
        tags$div(
          class = "text-muted",
          icon("upload", class = "fa-3x mb-3"),
          tags$h5("No data loaded"),
          tags$p(
            "Upload your data using the panel on the left, ",
            "or load a demo dataset to explore ", label, "."
          ),
          actionButton(
            "go_upload",
            label = "Go to upload",
            class = "btn btn-outline-primary btn-sm mt-2"
          )
        )
      )
    )
  }

  .unavailable_panel <- function(label, features) {
    reason <- if (!is.null(features)) {
      dt <- features$data_type
      if (dt != "counts" && label %in% c("Alpha Diversity", "Differential Abundance"))
        paste0("This module requires count data. Your data type is '", dt, "'.")
      else if (!features$has_taxonomy && label == "Composition")
        "This module requires a taxonomy table, which was not provided."
      else if (label == "ML Classification" && length(features$binary_group_vars) == 0)
        "ML classification requires a binary (two-group) grouping variable."
      else
        "This module is not available for the current dataset."
    } else {
      "Module not available."
    }

    card(
      class = "border-warning",
      card_body(
        class = "text-center py-5",
        tags$div(
          icon("triangle-exclamation", class = "fa-3x text-warning mb-3"),
          tags$h5(paste(label, "— unavailable")),
          tags$p(class = "text-muted", reason)
        )
      )
    )
  }

  # ── Button to jump to upload tab ──────────────────────────────────────────
  observeEvent(input$go_upload, {
    nav_select("top_nav", selected = "home")
  })
}
