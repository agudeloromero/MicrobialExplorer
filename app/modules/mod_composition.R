# =============================================================================
# mod_composition.R — Taxonomic Composition Module
# =============================================================================
# Wraps microbiome_composition.R functions:
#   agglomerate_taxa()          → used by all plots
#   plot_composition_bars()     → Stacked bars sub-tab
#   plot_mean_composition()     → Mean composition sub-tab
#   plot_abundance_heatmap()    → Heatmap sub-tab
#   calculate_fb_ratio()        → F:B ratio sub-tab
#   identify_core_microbiome()  → Core microbiome sub-tab
#   make_abundance_table()      → Abundance table (all sub-tabs)
# =============================================================================


# ── UI ────────────────────────────────────────────────────────────────────────
mod_composition_ui <- function(id) {
  ns <- NS(id)

  navset_card_tab(
    id = ns("comp_tabs"),

    # ── 1. Stacked bars ──────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("chart-bar"), " Stacked bars"),
      value = "bars",

      layout_columns(
        col_widths = c(3, 9),

        card(
          class = "border-0 bg-light",
          card_body(
            tags$h6(class = "fw-bold", "Parameters"),
            numericInput(ns("bars_top_n"), "Top N taxa",
                         value = 20, min = 5, max = 50, step = 5),
            selectInput(ns("bars_sort"), "Sort samples by",
                        choices = c("Group" = "group",
                                    "Dominant taxon" = "dominant_taxon",
                                    "None" = "none")),
            uiOutput(ns("bars_facet_ui")),
            hr(),
            downloadButton(ns("dl_bars"), "Download plot",
                           class = "btn btn-outline-secondary btn-sm w-100")
          )
        ),

        card(
          card_body(
            plotOutput(ns("bars_plot"), height = "550px")
          )
        )
      )
    ),

    # ── 2. Mean composition ──────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("chart-column"), " Mean composition"),
      value = "mean_comp",

      layout_columns(
        col_widths = c(3, 9),

        card(
          class = "border-0 bg-light",
          card_body(
            tags$h6(class = "fw-bold", "Parameters"),
            numericInput(ns("mean_top_n"), "Top N taxa",
                         value = 10, min = 3, max = 20, step = 1),
            hr(),
            uiOutput(ns("mean_comp_note")),
            hr(),
            downloadButton(ns("dl_mean"), "Download plot",
                           class = "btn btn-outline-secondary btn-sm w-100")
          )
        ),

        card(
          card_body(
            plotOutput(ns("mean_plot"), height = "500px")
          )
        )
      )
    ),

    # ── 3. Heatmap ───────────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("grid"), " Heatmap"),
      value = "heatmap",

      layout_columns(
        col_widths = c(3, 9),

        card(
          class = "border-0 bg-light",
          card_body(
            tags$h6(class = "fw-bold", "Parameters"),
            numericInput(ns("heat_top_n"), "Top N taxa",
                         value = 30, min = 5, max = 60, step = 5),
            selectInput(ns("heat_transform"), "Transformation",
                        choices = c("CLR" = "clr",
                                    "Log10" = "log10",
                                    "Relative abundance" = "relative"),
                        selected = "clr"),
            uiOutput(ns("heat_group_ui")),
            checkboxInput(ns("heat_cluster_rows"), "Cluster rows (taxa)",
                          value = TRUE),
            hr(),
            downloadButton(ns("dl_heat"), "Download plot",
                           class = "btn btn-outline-secondary btn-sm w-100")
          )
        ),

        card(
          card_body(
            plotOutput(ns("heat_plot"), height = "600px")
          )
        )
      )
    ),

    # ── 4. F:B ratio ─────────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("scale-balanced"), " F:B ratio"),
      value = "fb",

      layout_columns(
        col_widths = c(3, 9),

        card(
          class = "border-0 bg-light",
          card_body(
            tags$h6(class = "fw-bold", "About"),
            tags$p(class = "text-muted small",
                   "Firmicutes:Bacteroidota ratio per sample. ",
                   "Coloured by the selected group variable."),
            tags$p(class = "text-muted small fst-italic",
                   "Note: clinical significance of this ratio remains debated."),
            hr(),
            tags$h6(class = "fw-bold", "Summary"),
            uiOutput(ns("fb_summary")),
            hr(),
            downloadButton(ns("dl_fb"), "Download plot",
                           class = "btn btn-outline-secondary btn-sm w-100")
          )
        ),

        card(
          card_body(
            plotOutput(ns("fb_plot"), height = "500px")
          )
        )
      )
    ),

    # ── 5. Core microbiome ───────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("bullseye"), " Core microbiome"),
      value = "core",

      layout_columns(
        col_widths = c(3, 9),

        card(
          class = "border-0 bg-light",
          card_body(
            tags$h6(class = "fw-bold", "Parameters"),
            tags$p(class = "text-muted small",
                   "Taxa present at ≥ min abundance in the given fraction of samples."),
            numericInput(ns("core_min_abund"),
                         "Min relative abundance",
                         value = 0.001, min = 0.0001, max = 0.01,
                         step = 0.0005),
            hr(),
            tags$h6(class = "fw-bold", "Core taxa (≥ 75% samples)"),
            uiOutput(ns("core_taxa_list")),
            hr(),
            downloadButton(ns("dl_core"), "Download plot",
                           class = "btn btn-outline-secondary btn-sm w-100")
          )
        ),

        card(
          card_body(
            plotOutput(ns("core_plot"), height = "500px"),
            hr(),
            tags$h6(class = "fw-bold mt-2", "Core summary table"),
            DT::dataTableOutput(ns("core_table"))
          )
        )
      )
    ),

    # ── 6. Abundance table ───────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("table"), " Abundance table"),
      value = "table",

      layout_columns(
        col_widths = c(3, 9),

        card(
          class = "border-0 bg-light",
          card_body(
            tags$h6(class = "fw-bold", "Parameters"),
            uiOutput(ns("table_rank_ui")),
            hr(),
            downloadButton(ns("dl_table_overall"), "Download overall (.csv)",
                           class = "btn btn-outline-secondary btn-sm w-100 mb-2"),
            downloadButton(ns("dl_table_group"), "Download by group (.csv)",
                           class = "btn btn-outline-secondary btn-sm w-100")
          )
        ),

        card(
          card_body(
            tags$h6(class = "fw-bold mb-2", "Overall abundance summary"),
            DT::dataTableOutput(ns("abund_table"))
          )
        )
      )
    )
  )
}


# ── Server ────────────────────────────────────────────────────────────────────
mod_composition_server <- function(id, upload_data, ps_qc, parent_input) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Working phyloseq: prefer QC-filtered, fall back to raw ────────────
    ps_work <- reactive({
      ps_f <- if (is.function(ps_qc)) ps_qc() else ps_qc
      if (!is.null(ps_f)) ps_f else {
        req(upload_data$ps)
        upload_data$ps
      }
    })

    group_var <- reactive({
      gv <- parent_input$group_variable
      if (is.null(gv) || gv == "") NULL else gv
    })

    # Facet by — rendered as uiOutput so it exists before observe fires
    output$bars_facet_ui <- renderUI({
      feat <- upload_data$features
      req(feat)
      # Use group_vars (same source as sidebar group variable selector)
      # fall back to all sample_vars if needed
      facet_choices <- if (length(feat$group_vars) > 0) feat$group_vars
                       else if (length(feat$cat_vars) > 0) feat$cat_vars
                       else feat$sample_vars
      all_vars <- c("None" = "", facet_choices)
      selectInput(ns("bars_facet"), "Facet by",
                  choices = all_vars, selected = "")
    })

    # ── Rank selectors rendered as uiOutput so they populate correctly ────
    .rank_selector <- function(input_id, label = "Rank") {
      feat <- upload_data$features
      req(feat, feat$has_taxonomy)
      ranks     <- feat$avail_ranks
      preferred <- intersect(c("Genus", "Family", "Order", "Phylum"), ranks)
      default   <- if (length(preferred) > 0) preferred[1] else ranks[1]
      selectInput(ns(input_id), label, choices = ranks, selected = default)
    }

    output$table_rank_ui <- renderUI({ .rank_selector("table_rank") })

    # Heatmap group/sort selector
    output$heat_group_ui <- renderUI({
      feat <- upload_data$features
      req(feat)
      choices <- c("None" = "", feat$group_vars)
      selectInput(ns("heat_group"), "Sort columns by",
                  choices  = choices,
                  selected = "")
    })

    # Helper: resolve rank
    # bars/mean: always use global sidebar (no local selector)
    # heatmap/table: use local selector, fall back to global
    resolve_rank <- function(local_val = NULL) {
      if (!is.null(local_val) && nzchar(local_val)) return(local_val)
      grank <- parent_input$tax_rank
      if (!is.null(grank) && nzchar(grank)) return(grank)
      feat <- upload_data$features
      if (!is.null(feat) && length(feat$avail_ranks) > 0) {
        preferred <- intersect(c("Genus", "Family", "Order", "Phylum"), feat$avail_ranks)
        return(if (length(preferred) > 0) preferred[1] else feat$avail_ranks[1])
      }
      return("Genus")
    }



    # ── 1. STACKED BARS ───────────────────────────────────────────────────
    ps_bars <- reactive({
      req(ps_work())
      rank_use <- resolve_rank()   # uses global sidebar
      tryCatch(
        agglomerate_taxa(ps_work(),
                         rank      = rank_use,
                         transform = "relative",
                         top_n     = input$bars_top_n),
        error = function(e) { message("[mod_comp] bars agg: ", e$message); NULL }
      )
    })

    bars_plot_obj <- reactive({
      req(ps_bars())
      rank_use <- resolve_rank()
      facet <- if (is.null(input$bars_facet) || input$bars_facet == "") NULL
               else input$bars_facet
      tryCatch(
        plot_composition_bars(ps_bars(),
                              rank      = rank_use,
                              group_var = group_var(),
                              facet_var = facet,
                              sort_by   = input$bars_sort),
        error = function(e) { message("[mod_comp] bars plot: ", e$message); NULL }
      )
    })

    output$bars_plot <- renderPlot({
      req(bars_plot_obj())
      bars_plot_obj()
    })

    output$dl_bars <- downloadHandler(
      filename = function() paste0("composition_bars_", Sys.Date(), ".pdf"),
      content  = function(file) {
        req(bars_plot_obj())
        ggplot2::ggsave(file, bars_plot_obj(), width = 14, height = 7)
      }
    )

    # ── 2. MEAN COMPOSITION ───────────────────────────────────────────────
    ps_mean <- reactive({
      req(ps_work())
      rank_use <- resolve_rank()   # uses global sidebar
      tryCatch(
        agglomerate_taxa(ps_work(),
                         rank      = rank_use,
                         transform = "relative",
                         top_n     = input$mean_top_n),
        error = function(e) { message("[mod_comp] mean agg: ", e$message); NULL }
      )
    })

    mean_plot_obj <- reactive({
      req(ps_mean(), group_var())
      rank_use <- resolve_rank()
      tryCatch(
        plot_mean_composition(ps_mean(),
                              rank      = rank_use,
                              group_var = group_var(),
                              top_n     = input$mean_top_n),
        error = function(e) { message("[mod_comp] mean plot: ", e$message); NULL }
      )
    })

    output$mean_plot <- renderPlot({
      req(mean_plot_obj())
      mean_plot_obj()
    })

    output$mean_comp_note <- renderUI({
      if (is.null(group_var())) {
        tags$div(
          class = "alert alert-warning p-2 small",
          icon("triangle-exclamation"),
          " Select a group variable in the sidebar to enable this plot."
        )
      }
    })

    output$dl_mean <- downloadHandler(
      filename = function() paste0("mean_composition_", Sys.Date(), ".pdf"),
      content  = function(file) {
        req(mean_plot_obj())
        ggplot2::ggsave(file, mean_plot_obj(), width = 10, height = 7)
      }
    )

    # ── 3. HEATMAP ────────────────────────────────────────────────────────
    heat_plot_obj <- reactive({
      req(ps_work())
      rank_use <- resolve_rank()   # uses global sidebar
      # heat_group drives both annotation colour and within-group clustering
      hg <- if (is.null(input$heat_group) || input$heat_group == "") NULL
            else input$heat_group
      tryCatch(
        plot_abundance_heatmap(ps_work(),
                               rank            = rank_use,
                               top_n           = input$heat_top_n,
                               group_var       = hg,
                               transform       = input$heat_transform,
                               cluster_samples = !is.null(hg),
                               cluster_rows    = isTRUE(input$heat_cluster_rows)),
        error = function(e) { message("[mod_comp] heatmap: ", e$message); NULL }
      )
    })

    output$heat_plot <- renderPlot({
      req(heat_plot_obj())
      heat_plot_obj()
    })

    output$dl_heat <- downloadHandler(
      filename = function() paste0("abundance_heatmap_", Sys.Date(), ".pdf"),
      content  = function(file) {
        req(heat_plot_obj())
        ggplot2::ggsave(file, heat_plot_obj(), width = 14, height = 10)
      }
    )

    # ── 4. F:B RATIO ──────────────────────────────────────────────────────
    fb_result <- reactive({
      req(ps_work())
      tryCatch(
        calculate_fb_ratio(ps_work(), group_var = group_var()),
        error = function(e) { message("[mod_comp] fb ratio: ", e$message); NULL }
      )
    })

    output$fb_plot <- renderPlot({
      req(fb_result())
      fb_result()$plot
    })

    output$fb_summary <- renderUI({
      req(fb_result())
      df <- fb_result()$data
      med_fb <- round(median(df$FB_ratio, na.rm = TRUE), 2)
      tagList(
        stat_badge("Median F:B ratio", med_fb,
                   colour = if (!is.na(med_fb) && med_fb > 1) "warning" else "success"),
        stat_badge("Samples with ratio", sum(!is.na(df$FB_ratio)))
      )
    })

    output$dl_fb <- downloadHandler(
      filename = function() paste0("fb_ratio_", Sys.Date(), ".pdf"),
      content  = function(file) {
        req(fb_result())
        ggplot2::ggsave(file, fb_result()$plot, width = 8, height = 6)
      }
    )

    # ── 5. CORE MICROBIOME ────────────────────────────────────────────────
    core_result <- reactive({
      req(ps_work(), input$core_min_abund)
      tryCatch(
        identify_core_microbiome(ps_work(),
                                 prevalence_cuts = c(0.5, 0.75, 0.9),
                                 min_abundance   = input$core_min_abund,
                                 group_var       = group_var()),
        error = function(e) { message("[mod_comp] core: ", e$message); NULL }
      )
    })

    output$core_plot <- renderPlot({
      req(core_result())
      core_result()$plot
    })

    output$core_taxa_list <- renderUI({
      req(core_result())
      core_t <- core_result()$core_taxa
      if (length(core_t) == 0) {
        tags$span(class = "text-muted small", "No core taxa at 75% threshold")
      } else {
        tagList(
          tags$span(class = "badge bg-success mb-1",
                    paste(length(core_t), "taxa")),
          tags$ul(
            class = "small ps-3 mt-1",
            style = "max-height:150px; overflow-y:auto;",
            lapply(head(core_t, 15), function(t)
              tags$li(tags$em(t)))
          ),
          if (length(core_t) > 15)
            tags$p(class = "text-muted small",
                   paste0("… and ", length(core_t) - 15, " more"))
        )
      }
    })

    output$core_table <- DT::renderDataTable({
      req(core_result())
      core_result()$core_summary %>%
        mutate(prevalence_threshold = paste0(prevalence_threshold * 100, "%")) %>%
        rename(`Prevalence threshold` = prevalence_threshold,
               Group = group,
               `Core taxa` = n_core_taxa) %>%
        DT::datatable(options = list(dom = "t", pageLength = 20),
                      rownames = FALSE)
    })

    output$dl_core <- downloadHandler(
      filename = function() paste0("core_microbiome_", Sys.Date(), ".pdf"),
      content  = function(file) {
        req(core_result())
        ggplot2::ggsave(file, core_result()$plot, width = 10, height = 7)
      }
    )

    # ── 6. ABUNDANCE TABLE ────────────────────────────────────────────────
    abund_result <- reactive({
      req(ps_work())
      rank_use <- resolve_rank(input$table_rank)
      req(nchar(rank_use) > 0)
      tryCatch(
        make_abundance_table(ps_work(),
                             rank      = rank_use,
                             group_var = group_var()),
        error = function(e) { message("[mod_comp] abund table: ", e$message); NULL }
      )
    })

    output$abund_table <- DT::renderDataTable({
      req(abund_result())
      abund_result()$overall %>%
        rename(
          Taxon      = taxon,
          Lineage    = lineage,
          `Mean %`   = mean_pct,
          `SD %`     = sd_pct,
          `Median %` = median_pct,
          `Prev %`   = prev_pct,
          `N samples` = n_samples
        ) %>%
        DT::datatable(
          options = list(pageLength = 15, dom = "ftp",
                         order = list(list(2, "desc"))),
          rownames = FALSE
        ) %>%
        DT::formatRound(c("Mean %", "SD %", "Median %"), digits = 3)
    })

    output$dl_table_overall <- downloadHandler(
      filename = function() paste0("abundance_table_overall_", Sys.Date(), ".csv"),
      content  = function(file) {
        req(abund_result())
        write.csv(abund_result()$overall, file, row.names = FALSE)
      }
    )

    output$dl_table_group <- downloadHandler(
      filename = function() paste0("abundance_table_by_group_", Sys.Date(), ".csv"),
      content  = function(file) {
        req(abund_result())
        tbl <- abund_result()$by_group
        if (is.null(tbl)) {
          write.csv(data.frame(Message = "No group variable selected"), file,
                    row.names = FALSE)
        } else {
          write.csv(tbl, file, row.names = FALSE)
        }
      }
    )
  })
}
