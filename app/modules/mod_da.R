# =============================================================================
# mod_da.R — Differential Abundance Module
# =============================================================================
# Sub-tabs:
#   1. Run       — method selection, parameters, run button
#   2. Volcano   — volcano plot
#   3. Effects   — effect size lollipop plot
#   4. Heatmap   — DA taxa heatmap
#   5. Table     — full results table
# =============================================================================

# ── UI ────────────────────────────────────────────────────────────────────────
mod_da_ui <- function(id) {
  ns <- NS(id)

  navset_card_tab(
    id = ns("da_tabs"),

    # ── 1. Run ──────────────────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("play"), " Run"),
      layout_sidebar(
        sidebar = sidebar(
          width = 320, open = "always",
          tags$h6("Data preparation"),
          uiOutput(ns("rank_ui")),
          sliderInput(ns("min_prev"), "Min prevalence",
                      min = 0.05, max = 0.5, value = 0.10, step = 0.05),
          numericInput(ns("min_count"), "Min read count",
                       value = 10, min = 1, max = 100),
          hr(),
          tags$h6("Methods"),
          checkboxGroupInput(ns("methods"), NULL,
                             choices  = c("DESeq2"   = "deseq2",
                                          "ALDEx2"   = "aldex2",
                                          "ANCOM-BC" = "ancombc"),
                             selected = c("deseq2", "aldex2")),
          hr(),
          tags$h6("Parameters"),
          uiOutput(ns("reference_ui")),
          uiOutput(ns("formula_ui")),
          numericInput(ns("alpha"), "Significance threshold",
                       value = 0.05, min = 0.01, max = 0.2, step = 0.01),
          numericInput(ns("lfc_thresh"), "LFC threshold",
                       value = 1, min = 0, max = 4, step = 0.5),
          numericInput(ns("min_methods"), "Min methods for consensus",
                       value = 2, min = 1, max = 3, step = 1),
          hr(),
          actionButton(ns("run_da"), "Run analysis",
                       class = "btn btn-primary w-100",
                       icon  = icon("play"))
        ),
        uiOutput(ns("run_status"))
      )
    ),

    # ── 2. Volcano ──────────────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("chart-column"), " Volcano"),
      layout_sidebar(
        sidebar = sidebar(
          width = 280, open = "always",
          tags$h6("Parameters"),
          numericInput(ns("vol_top_n"), "Top N labels",
                       value = 15, min = 5, max = 40, step = 5),
          uiOutput(ns("vol_method_ui")),
          hr(),
          downloadButton(ns("vol_download"), "Download plot",
                         class = "btn btn-outline-secondary btn-sm w-100",
                         disabled = "disabled")
        ),
        plotOutput(ns("vol_plot"), height = "580px")
      )
    ),

    # ── 3. Effect sizes ─────────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("arrows-left-right"), " Effects"),
      layout_sidebar(
        sidebar = sidebar(
          width = 280, open = "always",
          tags$h6("Parameters"),
          numericInput(ns("eff_top_n"), "Top N taxa",
                       value = 30, min = 5, max = 60, step = 5),
          uiOutput(ns("eff_method_ui")),
          hr(),
          downloadButton(ns("eff_download"), "Download plot",
                         class = "btn btn-outline-secondary btn-sm w-100",
                         disabled = "disabled")
        ),
        plotOutput(ns("eff_plot"), height = "580px")
      )
    ),

    # ── 4. Heatmap ──────────────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("border-all"), " Heatmap"),
      layout_sidebar(
        sidebar = sidebar(
          width = 280, open = "always",
          tags$h6("Parameters"),
          selectInput(ns("heat_transform"), "Transformation",
                      choices  = c("Z-score" = "zscore",
                                   "CLR"     = "clr",
                                   "Log10"   = "log10"),
                      selected = "zscore"),
          uiOutput(ns("heat_sortby_ui")),
          checkboxInput(ns("heat_cluster_rows"), "Cluster rows (taxa)", value = TRUE),
          hr(),
          downloadButton(ns("heat_download"), "Download plot",
                         class = "btn btn-outline-secondary btn-sm w-100",
                         disabled = "disabled")
        ),
        plotOutput(ns("heat_plot"), height = "580px")
      )
    ),

    # ── 5. Table ────────────────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("table"), " Table"),
      layout_sidebar(
        sidebar = sidebar(
          width = 280, open = "always",
          tags$h6("Filter"),
          checkboxInput(ns("tab_sig_only"), "Significant only", value = FALSE),
          selectInput(ns("tab_view"), "View",
                      choices  = c("Consensus"   = "consensus",
                                   "All methods" = "all"),
                      selected = "consensus"),
          hr(),
          downloadButton(ns("tab_download"), "Download CSV",
                         class = "btn btn-outline-secondary btn-sm w-100",
                         disabled = "disabled")
        ),
        DT::DTOutput(ns("da_table"))
      )
    )
  )
}


# ── Server ────────────────────────────────────────────────────────────────────
mod_da_server <- function(id, upload_data, parent_input) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Reactive: phyloseq ───────────────────────────────────────────────────
    ps_work <- reactive({
      req(upload_data$ps)
      upload_data$ps
    })

    # ── Reactive: group variable ─────────────────────────────────────────────
    group_var <- reactive({
      gv <- parent_input$group_variable
      if (is.null(gv) || gv == "" || gv == "none") return(NULL)
      gv
    })

    # ── Dynamic UIs ──────────────────────────────────────────────────────────
    output$rank_ui <- renderUI({
      feat <- upload_data$features
      req(feat)
      ranks <- if (feat$has_taxonomy) feat$avail_ranks else "Genus"
      selectInput(ns("rank"), "Taxonomic rank",
                  choices = ranks,
                  selected = if ("Genus" %in% ranks) "Genus" else ranks[1])
    })

    output$reference_ui <- renderUI({
      gv <- group_var()
      req(gv, ps_work())
      levels_gv <- unique(as.character(sample_data(ps_work())[[gv]]))
      selectInput(ns("reference"), "Reference group",
                  choices = levels_gv, selected = levels_gv[1])
    })

    output$formula_ui <- renderUI({
      gv <- group_var()
      req(gv)
      textInput(ns("formula"), "ANCOM-BC formula",
                value = gv,
                placeholder = "e.g. disease_status + age + sex")
    })

    # Method selector for volcano/effects (populated after run)
    output$vol_method_ui <- renderUI({
      req(da_results())
      methods <- names(Filter(Negate(is.null),
                               list(deseq2  = da_results()$deseq2,
                                    aldex2  = da_results()$aldex2,
                                    ancombc = da_results()$ancombc)))
      selectInput(ns("vol_method"), "Method",
                  choices  = setNames(methods, toupper(methods)),
                  selected = methods[1])
    })

    output$eff_method_ui <- renderUI({
      req(da_results())
      methods <- names(Filter(Negate(is.null),
                               list(deseq2  = da_results()$deseq2,
                                    aldex2  = da_results()$aldex2,
                                    ancombc = da_results()$ancombc)))
      selectInput(ns("eff_method"), "Method",
                  choices  = setNames(methods, toupper(methods)),
                  selected = methods[1])
    })

    output$heat_sortby_ui <- renderUI({
      feat <- upload_data$features
      req(feat)
      choices <- c("Hierarchical clustering" = "none",
                   setNames(feat$group_vars, feat$group_vars))
      selectInput(ns("heat_sort_group"), "Sort columns by",
                  choices = choices, selected = if ("disease_status" %in% feat$group_vars) "disease_status" else choices[1])
    })

    # ── Run status UI ────────────────────────────────────────────────────────
    output$run_status <- renderUI({
      if (is.null(da_results())) {
        card(
          class = "border-0",
          card_body(
            class = "text-center py-5",
            tags$div(
              class = "text-muted",
              icon("play-circle", class = "fa-3x mb-3"),
              tags$h5("Ready to run"),
              tags$p("Configure parameters and click Run analysis.")
            )
          )
        )
      } else {
        res     <- da_results()
        cons    <- res$consensus$consensus
        n_sig   <- if (!is.null(cons)) sum(cons$consensus_da, na.rm = TRUE) else 0
        n_inc   <- if (!is.null(cons)) sum(cons$consensus_da & cons$direction == "increased", na.rm = TRUE) else 0
        n_dec   <- if (!is.null(cons)) sum(cons$consensus_da & cons$direction == "decreased", na.rm = TRUE) else 0

        tagList(
          layout_columns(
            col_widths = c(4, 4, 4),
            card(card_header("Consensus DA taxa"),
                 card_body(tags$h2(class = "text-center fw-bold text-primary", n_sig))),
            card(card_header(icon("arrow-up"), " Increased"),
                 card_body(tags$h2(class = "text-center fw-bold text-danger", n_inc))),
            card(card_header(icon("arrow-down"), " Decreased"),
                 card_body(tags$h2(class = "text-center fw-bold text-primary", n_dec)))
          ),
          tags$p(class = "text-muted small mt-2",
                 icon("circle-check", class = "text-success"),
                 " Analysis complete. View results in the tabs above.")
        )
      }
    })

    # ── Core reactive: run DA ────────────────────────────────────────────────
    da_results <- reactiveVal(NULL)

    observeEvent(input$run_da, {
      gv <- group_var()
      req(gv, ps_work())

      showNotification("Running differential abundance analysis...",
                       id = "da_running", duration = NULL, type = "message")

      tryCatch({
        rank     <- input$rank %||% "Genus"
        ps_da    <- prepare_da_data(ps_work(),
                                    rank           = rank,
                                    min_prevalence = input$min_prev,
                                    min_count      = input$min_count)

        results  <- list()
        methods  <- input$methods

        if ("deseq2" %in% methods) {
          results$deseq2 <- tryCatch(
            run_deseq2(ps_da, group_var = gv,
                       reference     = input$reference,
                       alpha         = input$alpha,
                       lfc_threshold = input$lfc_thresh),
            error = function(e) {
              message("[mod_da] DESeq2 error: ", e$message); NULL
            }
          )
        }

        if ("aldex2" %in% methods) {
          results$aldex2 <- tryCatch(
            run_aldex2(ps_da, group_var = gv, alpha = input$alpha),
            error = function(e) {
              message("[mod_da] ALDEx2 error: ", e$message); NULL
            }
          )
        }

        if ("ancombc" %in% methods) {
          results$ancombc <- tryCatch(
            run_ancombc(ps_da,
                        formula   = input$formula,
                        group_var = gv,
                        reference = input$reference,
                        alpha     = input$alpha),
            error = function(e) {
              message("[mod_da] ANCOM-BC error: ", e$message); NULL
            }
          )
        }

        results$consensus <- build_consensus(
          result_list = Filter(Negate(is.null),
                               list(ancombc = results$ancombc,
                                    deseq2  = results$deseq2,
                                    aldex2  = results$aldex2)),
          alpha       = input$alpha,
          min_methods = input$min_methods
        )

        results$ps_da <- ps_da
        results$rank  <- rank
        results$gv    <- gv

        da_results(results)
        removeNotification("da_running")
        showNotification("Analysis complete!", type = "message", duration = 3)

      }, error = function(e) {
        removeNotification("da_running")
        showNotification(paste("Error:", e$message), type = "error", duration = 8)
        message("[mod_da] run error: ", e$message)
      })
    })

    # ── Helper: get result df for chosen method ──────────────────────────────
    get_method_df <- function(method_input) {
      res <- da_results()
      req(res)
      df <- switch(method_input,
        deseq2  = res$deseq2,
        aldex2  = res$aldex2,
        ancombc = res$ancombc,
        res$consensus$all_results
      )
      req(df)
      df
    }

    # ── 2. Volcano ──────────────────────────────────────────────────────────
    vol_plot_obj <- reactive({
      req(da_results(), input$vol_method)
      df <- get_method_df(input$vol_method)
      tryCatch(
        plot_volcano(df,
                     alpha         = input$alpha,
                     lfc_threshold = input$lfc_thresh,
                     top_n_label   = input$vol_top_n,
                     title         = paste0("DA — ", toupper(input$vol_method))),
        error = function(e) { message("[mod_da] volcano error: ", e$message); NULL }
      )
    })

    output$vol_plot <- renderPlot({
      req(vol_plot_obj())
      vol_plot_obj()
    })

    output$vol_download <- downloadHandler(
      filename = function() paste0("volcano_", Sys.Date(), ".pdf"),
      content  = function(file) {
        req(vol_plot_obj())
        ggplot2::ggsave(file, vol_plot_obj(), width = 10, height = 8)
      }
    )

    # ── 3. Effect sizes ─────────────────────────────────────────────────────
    eff_plot_obj <- reactive({
      req(da_results(), input$eff_method)
      df <- get_method_df(input$eff_method)
      gv <- da_results()$gv
      tryCatch(
        plot_effect_sizes(df,
                          alpha            = input$alpha,
                          top_n            = input$eff_top_n,
                          group_comparison = gv),
        error = function(e) { message("[mod_da] effects error: ", e$message); NULL }
      )
    })

    output$eff_plot <- renderPlot({
      req(eff_plot_obj())
      eff_plot_obj()
    })

    output$eff_download <- downloadHandler(
      filename = function() paste0("effect_sizes_", Sys.Date(), ".pdf"),
      content  = function(file) {
        req(eff_plot_obj())
        ggplot2::ggsave(file, eff_plot_obj(), width = 10, height = 10)
      }
    )

    # ── 4. Heatmap ──────────────────────────────────────────────────────────
    heat_plot_obj <- reactive({
      req(da_results())
      res     <- da_results()
      cons    <- res$consensus$consensus
      req(cons)
      sig_taxa <- cons %>%
        dplyr::filter(consensus_da) %>%
        dplyr::pull(taxon)
      if (length(sig_taxa) == 0) return(NULL)
      # Use plot_abundance_heatmap (same as composition — has dendrogram)
      tryCatch({
        ps_agg <- tax_glom(res$ps_da, taxrank = res$rank, NArm = FALSE)
        taxa_names(ps_agg) <- as.character(tax_table(ps_agg)[, res$rank])
        shared <- intersect(sig_taxa, taxa_names(ps_agg))
        if (length(shared) == 0) return(NULL)
        ps_filt <- prune_taxa(shared, ps_agg)
        sort_val        <- if (!is.null(input$heat_sort_group)) input$heat_sort_group else "none"
        gv_heat         <- if (sort_val != "none") sort_val else NULL
        cluster_cols    <- is.null(gv_heat)
        plot_abundance_heatmap(
          ps_filt,
          rank            = res$rank,
          top_n           = ntaxa(ps_filt),
          group_var       = gv_heat,
          transform       = input$heat_transform,
          cluster_samples = cluster_cols,
          cluster_rows    = isTRUE(input$heat_cluster_rows)
        )
      }, error = function(e) {
        message("[mod_da] heatmap error: ", e$message)
        NULL
      })
    })
    output$heat_plot <- renderPlot({
      p <- heat_plot_obj()
      if (is.null(p)) {
        plot.new()
        text(0.5, 0.5, "No consensus DA taxa found.
Try relaxing thresholds.",
             cex = 1.2, col = "grey50")
      } else {
        print(p)
      }
    })

    output$heat_download <- downloadHandler(
      filename = function() paste0("da_heatmap_", Sys.Date(), ".pdf"),
      content  = function(file) {
        req(heat_plot_obj())
        ggplot2::ggsave(file, heat_plot_obj(), width = 12, height = 8)
      }
    )

    # ── 5. Table ────────────────────────────────────────────────────────────
    table_df <- reactive({
      req(da_results())
      res  <- da_results()
      cons <- res$consensus

      df <- if (input$tab_view == "consensus") cons$consensus else cons$all_results
      req(df)

      if (isTRUE(input$tab_sig_only)) {
        if ("consensus_da" %in% colnames(df)) df <- df[df$consensus_da == TRUE, ]
        else if ("diff_abund" %in% colnames(df)) df <- df[df$diff_abund == TRUE, ]
      }

      num_cols <- sapply(df, is.numeric)
      df[num_cols] <- lapply(df[num_cols], round, 4)
      df
    })

    output$da_table <- DT::renderDT({
      req(table_df())
      DT::datatable(table_df(),
                    options = list(pageLength = 20, scrollX = TRUE),
                    rownames = FALSE,
                    filter   = "top")
    })

    output$tab_download <- downloadHandler(
      filename = function() paste0("da_results_", Sys.Date(), ".csv"),
      content  = function(file) {
        req(table_df())
        write.csv(table_df(), file, row.names = FALSE)
      }
    )

    # ── Enable download buttons ──────────────────────────────────────────────
    observe({
      shinyjs::toggleState("vol_download",  !is.null(vol_plot_obj()))
      shinyjs::toggleState("eff_download",  !is.null(eff_plot_obj()))
      shinyjs::toggleState("heat_download", !is.null(heat_plot_obj()))
      shinyjs::toggleState("tab_download",  !is.null(da_results()))
    })
  })
}

