# =============================================================================
# mod_qc.R — Quality Control Module
# =============================================================================
# Wraps microbiome_qc.R functions:
#   qc_sequencing_depth()    → Depth sub-tab
#   qc_rarefaction_curves()  → Rarefaction sub-tab
#   qc_filter_taxa()         → Taxa filtering sub-tab
#   qc_filter_samples()      → Sample filtering sub-tab
#
# Returns: rv$ps_filtered (filtered phyloseq) for downstream modules
# =============================================================================


# ── UI ────────────────────────────────────────────────────────────────────────
mod_qc_ui <- function(id) {
  ns <- NS(id)

  tagList(
    # ── Sub-tab controls injected into the sidebar when QC is active ─────
    # (These are rendered by mod_qc_sidebar_ui and swapped in by server.R)

    navset_card_tab(
      id = ns("qc_tabs"),

      # ── 1. Sequencing depth ─────────────────────────────────────────────
      nav_panel(
        title = tagList(icon("chart-bar"), " Depth"),
        value = "depth",

        layout_columns(
          col_widths = c(3, 9),

          # Controls
          card(
            class = "border-0 bg-light",
            card_body(
              tags$h6(class = "fw-bold", "Parameters"),
              numericInput(
                ns("min_reads"),
                label = tooltip(
                  tags$span("Min reads / sample ", icon("circle-question", class = "small")),
                  "Samples below this threshold are flagged in red."
                ),
                value = 1000, min = 0, step = 100
              ),
              hr(),
              tags$h6(class = "fw-bold", "Summary"),
              uiOutput(ns("depth_summary_badges")),
              hr(),
              tags$h6(class = "fw-bold", "Flagged samples"),
              uiOutput(ns("flagged_samples_ui"))
            )
          ),

          # Plot
          card(
            card_body(
              plotOutput(ns("depth_plot"), height = "500px")
            )
          )
        )
      ),

      # ── 2. Rarefaction curves ────────────────────────────────────────────
      nav_panel(
        title = tagList(icon("chart-line"), " Rarefaction"),
        value = "rarefaction",

        layout_columns(
          col_widths = c(3, 9),

          card(
            class = "border-0 bg-light",
            card_body(
              tags$h6(class = "fw-bold", "Parameters"),
              numericInput(
                ns("rare_step"),
                label = tooltip(
                  tags$span("Step size ", icon("circle-question", class = "small")),
                  "Interval between rarefaction depths. Smaller = smoother curves, slower."
                ),
                value = 500, min = 100, step = 100
              ),
              numericInput(
                ns("rare_n_samples"),
                label = tooltip(
                  tags$span("Max samples to plot ", icon("circle-question", class = "small")),
                  "Large datasets are subsampled for speed. Leave blank to plot all."
                ),
                value = NA, min = 5
              ),
              actionButton(
                ns("run_rarefaction"),
                label = tagList(icon("play"), " Run"),
                class = "btn btn-primary btn-sm w-100 mt-2"
              ),
              hr(),
              tags$h6(class = "fw-bold", "Knee point"),
              uiOutput(ns("rare_knee_ui"))
            )
          ),

          card(
            card_body(
              plotOutput(ns("rare_plot"), height = "500px")
            )
          )
        )
      ),

      # ── 3. Taxa filtering ────────────────────────────────────────────────
      nav_panel(
        title = tagList(icon("filter"), " Taxa filter"),
        value = "taxa_filter",

        layout_columns(
          col_widths = c(3, 9),

          card(
            class = "border-0 bg-light",
            card_body(
              tags$h6(class = "fw-bold", "Filter thresholds"),
              tags$p(class = "text-muted small",
                     "These mirror the global sidebar sliders."),
              sliderInput(
                ns("taxa_prev"),
                label = "Min prevalence",
                min = 0.01, max = 0.5,
                value = 0.10, step = 0.01,
                ticks = FALSE
              ),
              sliderInput(
                ns("taxa_abund"),
                label = "Min total abundance",
                min = 1, max = 100,
                value = 10, step = 1,
                ticks = FALSE
              ),
              actionButton(
                ns("apply_taxa_filter"),
                label = tagList(icon("check"), " Apply filter"),
                class = "btn btn-success btn-sm w-100 mt-2"
              ),
              hr(),
              tags$h6(class = "fw-bold", "Result"),
              uiOutput(ns("taxa_filter_summary"))
            )
          ),

          card(
            card_body(
              plotOutput(ns("taxa_filter_plot"), height = "500px")
            )
          )
        )
      ),

      # ── 4. Sample filtering ──────────────────────────────────────────────
      nav_panel(
        title = tagList(icon("users-slash"), " Sample filter"),
        value = "sample_filter",

        layout_columns(
          col_widths = c(3, 9),

          card(
            class = "border-0 bg-light",
            card_body(
              tags$h6(class = "fw-bold", "Filter thresholds"),
              numericInput(
                ns("samp_min_reads"),
                label = "Min reads",
                value = 1000, min = 0, step = 100
              ),
              numericInput(
                ns("samp_min_taxa"),
                label = "Min observed taxa",
                value = 10, min = 1, step = 1
              ),
              actionButton(
                ns("apply_sample_filter"),
                label = tagList(icon("check"), " Apply filter"),
                class = "btn btn-success btn-sm w-100 mt-2"
              ),
              hr(),
              tags$h6(class = "fw-bold", "Result"),
              uiOutput(ns("sample_filter_summary"))
            )
          ),

          card(
            card_body(
              tags$h6(class = "fw-bold mb-3", "Removed samples"),
              DT::dataTableOutput(ns("removed_samples_table")),
              hr(),
              tags$h6(class = "fw-bold mb-3", "Retained samples"),
              DT::dataTableOutput(ns("retained_samples_table"))
            )
          )
        )
      ),

      # ── 5. Final filtered dataset ────────────────────────────────────────
      nav_panel(
        title = tagList(icon("circle-check"), " Final dataset"),
        value = "final",

        layout_columns(
          col_widths = c(6, 6),

          card(
            card_header(icon("table"), " QC summary table"),
            card_body(
              DT::dataTableOutput(ns("qc_summary_table"))
            )
          ),

          card(
            card_header(icon("chart-scatter"), " Reads vs taxon richness"),
            card_body(
              plotOutput(ns("final_plot"), height = "350px")
            )
          )
        ),

        card(
          card_header(icon("download"), " Export"),
          card_body(
            tags$p(class = "text-muted small",
                   "Download the filtered phyloseq object for use in other modules or your own analysis."),
            downloadButton(ns("download_ps"), "Download filtered phyloseq (.rds)",
                           class = "btn btn-outline-primary btn-sm")
          )
        )
      )
    )
  )
}


# ── Server ────────────────────────────────────────────────────────────────────
mod_qc_server <- function(id, upload_data, parent_input) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Convenience reactives ─────────────────────────────────────────────
    ps_raw <- reactive({
      req(upload_data$ps)
      upload_data$ps
    })

    group_var <- reactive({
      gv <- parent_input$group_variable
      if (is.null(gv) || gv == "") NULL else gv
    })

    # ── 1. SEQUENCING DEPTH — auto-reactive ───────────────────────────────
    depth_result <- reactive({
      req(ps_raw(), input$min_reads)
      tryCatch(
        qc_sequencing_depth(ps_raw(),
                            min_reads = input$min_reads,
                            group_var = group_var()),
        error = function(e) { message("[mod_qc] depth: ", e$message); NULL }
      )
    })

    output$depth_plot <- renderPlot({
      req(depth_result())
      depth_result()$plot
    })

    output$depth_summary_badges <- renderUI({
      req(depth_result())
      s <- depth_result()$summary
      tagList(
        stat_badge("Samples",      s$n_samples),
        stat_badge("Min reads",    format(s$min_reads,    big.mark = ",")),
        stat_badge("Median reads", format(s$median_reads, big.mark = ",")),
        stat_badge("Max reads",    format(s$max_reads,    big.mark = ",")),
        stat_badge("Below threshold",
                   paste0(s$n_below_threshold, " (", s$pct_below, "%)"),
                   colour = if (s$n_below_threshold > 0) "warning" else "success")
      )
    })

    output$flagged_samples_ui <- renderUI({
      req(depth_result())
      flagged <- depth_result()$flagged
      if (length(flagged) == 0) {
        tags$span(class = "text-success small", icon("check"), " No samples flagged")
      } else {
        tagList(
          tags$span(class = "text-danger small",
                    icon("triangle-exclamation"),
                    paste(length(flagged), "sample(s) below threshold:")),
          tags$ul(class = "small ps-3 mt-1", lapply(flagged, tags$li))
        )
      }
    })

    # ── 2. RAREFACTION — manual run button ───────────────────────────────
    rare_result <- eventReactive(
      list(input$run_rarefaction, ps_raw()),
      {
        req(ps_raw(), input$rare_step)
        n_samp <- if (is.null(input$rare_n_samples) ||
                      is.na(input$rare_n_samples)) NULL else as.integer(input$rare_n_samples)
        tryCatch(
          qc_rarefaction_curves(ps_raw(),
                                step      = input$rare_step,
                                group_var = group_var(),
                                n_samples = n_samp),
          error = function(e) { message("[mod_qc] rarefaction: ", e$message); NULL }
        )
      },
      ignoreNULL = TRUE
    )

    output$rare_plot <- renderPlot({
      req(rare_result())
      rare_result()$plot
    })

    output$rare_knee_ui <- renderUI({
      req(rare_result())
      median_knee <- median(rare_result()$knee$knee_depth)
      tagList(
        stat_badge("Median plateau",
                   format(round(median_knee), big.mark = ","),
                   colour = "info"),
        tags$p(class = "text-muted small mt-1",
               "Suggested rarefaction depth: ",
               tags$strong(format(round(median_knee), big.mark = ",")), " reads")
      )
    })

    # ── 3. TAXA FILTERING — auto-reactive, Apply button re-runs ──────────
    # Tracks the thresholds that were last applied
    taxa_thresholds <- reactiveValues(prev = 0.10, abund = 10)

    taxa_result <- reactive({
      req(ps_raw())
      # Depend on Apply button OR initial data load
      input$apply_taxa_filter
      isolate({
        tryCatch(
          qc_filter_taxa(ps_raw(),
                         min_prevalence = taxa_thresholds$prev,
                         min_abundance  = taxa_thresholds$abund),
          error = function(e) { message("[mod_qc] taxa filter: ", e$message); NULL }
        )
      })
    })

    observeEvent(input$apply_taxa_filter, {
      taxa_thresholds$prev  <- input$taxa_prev
      taxa_thresholds$abund <- input$taxa_abund
    })

    output$taxa_filter_plot <- renderPlot({
      req(taxa_result())
      taxa_result()$plot
    })

    output$taxa_filter_summary <- renderUI({
      req(taxa_result())
      s <- taxa_result()$summary
      tagList(
        stat_badge("Taxa before", s$n_before),
        stat_badge("Taxa after",  s$n_after, colour = "success"),
        stat_badge("Removed",
                   paste0(s$n_removed, " (", s$pct_removed, "%)"),
                   colour = if (s$n_removed > 0) "warning" else "success")
      )
    })

    # ── 4. SAMPLE FILTERING — auto-reactive, Apply button re-runs ────────
    samp_thresholds <- reactiveValues(reads = 1000, taxa = 10)

    sample_result <- reactive({
      req(ps_raw())
      input$apply_sample_filter
      isolate({
        tryCatch(
          qc_filter_samples(ps_raw(),
                            min_reads = samp_thresholds$reads,
                            min_taxa  = samp_thresholds$taxa),
          error = function(e) { message("[mod_qc] sample filter: ", e$message); NULL }
        )
      })
    })

    observeEvent(input$apply_sample_filter, {
      samp_thresholds$reads <- input$samp_min_reads
      samp_thresholds$taxa  <- input$samp_min_taxa
    })

    output$sample_filter_summary <- renderUI({
      req(sample_result())
      n_removed <- length(sample_result()$removed_samples)
      tagList(
        stat_badge("Samples before", nsamples(ps_raw())),
        stat_badge("Samples after",  nsamples(sample_result()$ps_filtered),
                   colour = "success"),
        stat_badge("Removed", n_removed,
                   colour = if (n_removed > 0) "warning" else "success")
      )
    })

    output$removed_samples_table <- DT::renderDataTable({
      req(sample_result())
      removed <- sample_result()$removed_samples
      if (length(removed) == 0)
        return(DT::datatable(data.frame(Message = "No samples removed"),
                             rownames = FALSE, options = list(dom = "t")))
      data.frame(sample_data(ps_raw()))[removed, , drop = FALSE] %>%
        tibble::rownames_to_column("Sample") %>%
        DT::datatable(options = list(pageLength = 10, dom = "tp"), rownames = FALSE)
    })

    output$retained_samples_table <- DT::renderDataTable({
      req(sample_result())
      data.frame(sample_data(sample_result()$ps_filtered)) %>%
        tibble::rownames_to_column("Sample") %>%
        DT::datatable(options = list(pageLength = 10, dom = "tp"), rownames = FALSE)
    })

    # ── 5. FINAL DATASET — derived from both filter reactives ─────────────
    ps_final <- reactive({
      req(ps_raw())
      ps_f <- ps_raw()

      if (!is.null(sample_result())) {
        keep_samps <- sample_names(sample_result()$ps_filtered)
        ps_f <- prune_samples(keep_samps, ps_f)
      }

      if (!is.null(taxa_result())) {
        keep_taxa <- intersect(taxa_names(taxa_result()$ps_filtered), taxa_names(ps_f))
        if (length(keep_taxa) > 0)
          ps_f <- prune_taxa(keep_taxa, ps_f)
      }

      ps_f
    })

    output$qc_summary_table <- DT::renderDataTable({
      req(ps_final())
      ps_f <- ps_final()
      df <- data.frame(
        Step = c("Raw data", "After sample filtering",
                 "After taxa filtering", "Final dataset"),
        Samples = c(
          nsamples(ps_raw()),
          if (!is.null(sample_result())) nsamples(sample_result()$ps_filtered) else NA,
          NA,
          nsamples(ps_f)
        ),
        Taxa = c(
          ntaxa(ps_raw()),
          NA,
          if (!is.null(taxa_result())) taxa_result()$summary$n_after else NA,
          ntaxa(ps_f)
        ),
        Total_reads = c(
          format(sum(sample_sums(ps_raw())), big.mark = ","),
          NA, NA,
          format(sum(sample_sums(ps_f)), big.mark = ",")
        ),
        Median_reads = c(
          format(round(median(sample_sums(ps_raw()))), big.mark = ","),
          NA, NA,
          format(round(median(sample_sums(ps_f))), big.mark = ",")
        )
      )
      DT::datatable(df, options = list(dom = "t"), rownames = FALSE) %>%
        DT::formatStyle("Step", fontWeight = "bold")
    })

    output$final_plot <- renderPlot({
      req(ps_final())
      ps_f  <- ps_final()
      otu_m <- as(otu_table(ps_f), "matrix")
      if (taxa_are_rows(ps_f)) otu_m <- t(otu_m)
      final_df <- data.frame(
        reads  = sample_sums(ps_f),
        n_taxa = rowSums(otu_m > 0)
      )
      ggplot(final_df, aes(x = reads, y = n_taxa)) +
        geom_point(colour = "#27ae60", alpha = 0.75, size = 3) +
        geom_smooth(method = "lm", se = TRUE, colour = "#2980b9",
                    linetype = "dashed", linewidth = 0.8, alpha = 0.15) +
        scale_x_continuous(labels = scales::label_comma()) +
        labs(
          title    = "Final dataset: reads vs taxon richness",
          subtitle = paste0(nsamples(ps_f), " samples | ", ntaxa(ps_f), " taxa | ",
                            format(sum(sample_sums(ps_f)), big.mark = ","), " total reads"),
          x = "Total reads per sample",
          y = "Observed taxa per sample"
        ) +
        theme_microbiome()
    })

    # ── Download ──────────────────────────────────────────────────────────
    output$download_ps <- downloadHandler(
      filename = function()
        paste0("MicrobialExplorer_QC_filtered_", format(Sys.Date(), "%Y%m%d"), ".rds"),
      content  = function(file) {
        req(ps_final())
        saveRDS(ps_final(), file)
      }
    )

    # ── Return filtered phyloseq to parent ────────────────────────────────
    return(ps_final)
  })
}
