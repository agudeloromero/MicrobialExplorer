# =============================================================================
# mod_alpha.R — Alpha Diversity Module
# =============================================================================

# ── UI ────────────────────────────────────────────────────────────────────────
mod_alpha_ui <- function(id) {
  ns <- NS(id)

  navset_card_tab(
    id = ns("alpha_tabs"),

    # ── 1. Boxplots ─────────────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("chart-simple"), " Boxplots"),
      value = "box",
      layout_sidebar(
        sidebar = sidebar(
          width = 320, open = "always",
          tags$h6("Parameters"),
          checkboxGroupInput(
            ns("box_metrics"),
            "Diversity metrics",
            choices  = c("Observed richness" = "observed",
                         "Chao1"             = "chao1",
                         "Shannon"           = "shannon",
                         "Simpson"           = "simpson",
                         "Pielou's evenness" = "pielou",
                         "Faith's PD"        = "faith_pd"),
            selected = c("observed", "shannon", "simpson", "pielou")
          ),
          selectInput(
            ns("box_test"), "Statistical test",
            choices = c(
              "Wilcoxon rank-sum (non-parametric, unpaired)"   = "wilcox_unpaired",
              "Wilcoxon signed-rank (non-parametric, paired)"  = "wilcox_paired",
              "Welch's t-test (parametric, unpaired)"          = "ttest_unpaired",
              "Paired t-test (parametric, paired)"             = "ttest_paired"
            ),
            selected = "wilcox_unpaired"
          ),
          uiOutput(ns("box_pair_ui")),
          selectInput(ns("box_comparisons"), "Comparisons on the plot",
                      choices  = c("Significant only" = "significant",
                                   "All pairs"        = "all",
                                   "None"             = "none"),
                      selected = "significant"),
          checkboxInput(ns("box_points"), "Show data points", value = TRUE),
          actionButton(ns("box_reset"),
                       label = tagList(icon("rotate-left"), " Reset to defaults"),
                       class = "btn btn-outline-secondary btn-sm w-100 mt-2"),
          hr(),
          downloadButton(ns("box_download"), "Download plot",
                         class = "btn btn-outline-secondary btn-sm w-100 mb-2",
                         disabled = "disabled"),
          downloadButton(ns("box_stats_download"), "Download statistics (.csv)",
                         class = "btn btn-outline-secondary btn-sm w-100")
        ),
        uiOutput(ns("alpha_depth_note")),
        plotOutput(ns("box_plot"), height = "550px"),
        hr(),
        tags$h6(class = "fw-bold", "Group comparisons"),
        DT::DTOutput(ns("box_stats")),
        hr(),

        # Help
        card(
          card_header(icon("circle-info"), " How to read this tab"),
          card_body(
            class = "small",
            layout_columns(
              col_widths = c(6, 6),
              tags$div(
                tags$h6(class = "fw-bold", "What alpha diversity is"),
                tags$ul(
                  tags$li(tags$strong("Diversity inside one sample. "),
                          "Each dot is a sample, each box a group. Alpha diversity ",
                          "asks how varied one community is; comparing whole ",
                          "communities between samples is the Beta tab."),
                  tags$li(tags$strong("The metrics measure different things. "),
                          "Observed richness counts taxa. Shannon combines richness ",
                          "with evenness. Simpson is dominated by the common taxa, so ",
                          "it is insensitive to rare ones. Pielou's evenness asks how ",
                          "equally abundance is shared. Faith's PD counts the branch ",
                          "length of the tree covered, so it needs a phylogenetic tree. ",
                          "Reporting one from each family gives a fuller picture than ",
                          "reporting several that measure the same thing."),
                  tags$li(tags$strong("Chao1 and ACE need singletons. "),
                          "They estimate unseen taxa from taxa seen once or twice. ",
                          "Pipelines that build ASVs, such as DADA2, remove singletons, ",
                          "so for ASV data these two are not interpretable. They are ",
                          "fine for OTU tables that keep singletons."),
                  tags$li(tags$strong("Everything here is computed at one depth. "),
                          "Counts are subsampled to the shallowest sample and averaged ",
                          "over ten draws, which is current best practice for comparing ",
                          "diversity between samples of unequal depth.")
                )
              ),
              tags$div(
                tags$h6(class = "fw-bold", "The tests, and how to choose"),
                tags$ul(
                  tags$li(tags$strong("Groups are compared two at a time. "),
                          "There is no overall test: each pair gets its own p-value and ",
                          "the whole set is corrected together (Benjamini-Hochberg)."),
                  tags$li(tags$strong("Non-parametric (Wilcoxon) is the default "),
                          "and is what most microbiome papers use, because diversity ",
                          "values are often skewed and groups are small. Choose ",
                          "parametric (t-test) only if you have checked that the values ",
                          "are roughly normal within each group."),
                  tags$li(tags$strong("Unpaired or paired. "),
                          "Paired means the same subject appears in both groups, such ",
                          "as before and after treatment; then choose the variable that ",
                          "identifies the subject, and only subjects present in both ",
                          "groups are used."),
                  tags$li(tags$strong("Small groups cannot reach significance. "),
                          "With three samples against three, the smallest possible ",
                          "p-value from a rank test is 0.10, so a result of ",
                          "\"not significant\" there says nothing about the biology. ",
                          "Pairs that are too small are listed in the table with a note ",
                          "instead of a p-value."),
                  tags$li(tags$strong("Many groups cost power. "),
                          "Nine groups make 36 comparisons, and correcting for all of ",
                          "them raises every p-value. If only two groups matter to your ",
                          "question, compare those by choosing a group variable with ",
                          "fewer levels."),
                  tags$li(tags$strong("Brackets show the adjusted p-value. "),
                          "By default only significant pairs are drawn; the full table ",
                          "below lists every pair and can be downloaded.")
                )
              )
            )
          )
        )
      )
    ),

    # ── 3. Diversity gradient ────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("chart-area"), " Gradient"),
      value = "gradient",
      layout_sidebar(
        sidebar = sidebar(
          width = 320, open = "always",
          tags$h6("Parameters"),
          selectInput(ns("grad_metric"), "Diversity metric",
                      choices  = c("Observed richness" = "observed",
                                   "Chao1"             = "chao1",
                                   "Shannon"           = "shannon",
                                   "Simpson"           = "simpson",
                                   "Pielou's evenness" = "pielou",
                                   "Faith's PD"        = "faith_pd"),
                      selected = "shannon"),
          uiOutput(ns("grad_xvar_ui")),
          checkboxInput(ns("grad_smooth"), "Add smoothing line", value = TRUE),
          selectInput(ns("grad_method"), "Smoothing method",
                      choices  = c("LOESS"  = "loess",
                                   "Linear" = "lm",
                                   "GAM"    = "gam"),
                      selected = "loess"),
          hr(),
          downloadButton(ns("grad_download"), "Download plot",
                         class = "btn btn-outline-secondary btn-sm w-100",
                         disabled = "disabled")
        ),
        plotOutput(ns("grad_plot"), height = "550px")
      )
    ),

    # ── 4. Diversity table ───────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("table"), " Table"),
      value = "table",
      layout_sidebar(
        sidebar = sidebar(
          width = 280, open = "always",
          tags$h6("Export"),
          downloadButton(ns("tab_download"), "Download CSV",
                         class = "btn btn-outline-secondary btn-sm w-100",
                         disabled = "disabled")
        ),
        DT::DTOutput(ns("alpha_table"))
      )
    )
  )
}


# ── Server ────────────────────────────────────────────────────────────────────
mod_alpha_server <- function(id, upload_data, parent_input) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Reactive: phyloseq object ────────────────────────────────────────────
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

    # ── Reactive: diversity data frame ───────────────────────────────────────
    diversity_df <- reactive({
      req(ps_work())
      tryCatch({
        calculate_alpha_diversity(ps_work(), rarefaction = TRUE)
      }, error = function(e) {
        message("[mod_alpha] diversity calc error: ", e$message)
        NULL
      })
    })

    # ── 1. Boxplots ──────────────────────────────────────────────────────────
    # The chosen test, split into the two things the analysis needs
    box_test_kind <- reactive({
      v <- input$box_test
      if (is.null(v)) v <- "wilcox_unpaired"
      list(
        test   = if (startsWith(v, "wilcox")) "wilcox" else "t.test",
        paired = endsWith(v, "_paired")
      )
    })

    # Subject variable, shown only for the paired tests
    output$box_pair_ui <- renderUI({
      req(box_test_kind()$paired, upload_data$ps)
      meta_df <- data.frame(sample_data(upload_data$ps), check.names = FALSE)
      gv      <- group_var()
      n_uni   <- vapply(meta_df, function(x) length(unique(stats::na.omit(x))),
                        integer(1))
      vars    <- setdiff(names(meta_df)[n_uni >= 2 & n_uni < nrow(meta_df)], gv)
      selectInput(ns("box_pair_var"), "Subject / pairing variable",
                  choices = c("Choose a variable" = "", vars), selected = "")
    })

    # Depth actually used for every metric on this tab
    output$alpha_depth_note <- renderUI({
      req(ps_work())
      depths <- sample_sums(ps_work())
      depth  <- min(depths)
      med    <- stats::median(depths)
      n_low  <- sum(depths < med / 2)
      tags$div(
        class = "alert alert-warning small",
        icon("triangle-exclamation"),
        tags$strong(" All metrics are computed at ",
                    format(round(depth), big.mark = ","), " reads"),
        ", the depth of your shallowest sample, with counts subsampled to that ",
        "depth and averaged over ten draws. ",
        if (n_low > 0) {
          tagList(
            tags$strong(n_low, if (n_low == 1) " sample has" else " samples have"),
            " fewer than half the median depth of ",
            format(round(med), big.mark = ","),
            " reads and are pulling this threshold down. Remove them in ",
            tags$strong("QC \u2192 Sample filter"),
            " and come back, or the diversity of every other sample is measured ",
            "on a fraction of its reads."
          )
        } else {
          "Sample depths are even enough that little is lost."
        }
      )
    })

    box_result <- reactive({
      req(diversity_df())
      gv <- group_var()
      req(gv, length(input$box_metrics) > 0)
      paired <- box_test_kind()$paired
      pv     <- if (paired) input$box_pair_var else NULL
      if (paired && (is.null(pv) || !nzchar(pv))) return(NULL)
      tryCatch({
        plot_alpha_diversity(diversity_df(),
                             metrics     = input$box_metrics,
                             group_var   = gv,
                             test        = box_test_kind()$test,
                             paired      = paired,
                             pair_var    = pv,
                             comparisons = input$box_comparisons,
                             show_points = input$box_points)
      }, error = function(e) {
        message("[mod_alpha] boxplot error: ", e$message)
        NULL
      })
    })

    observeEvent(input$box_reset, {
      updateCheckboxGroupInput(session, "box_metrics",
                               selected = c("observed", "shannon",
                                            "simpson", "pielou"))
      updateSelectInput(session, "box_test", selected = "wilcox_unpaired")
      updateSelectInput(session, "box_comparisons", selected = "significant")
      updateCheckboxInput(session, "box_points", value = TRUE)
    })

    box_plot_obj <- reactive({
      req(box_result())
      box_result()$plot
    })

    output$box_plot <- renderPlot({
      validate(need(!box_test_kind()$paired ||
                      (!is.null(input$box_pair_var) && nzchar(input$box_pair_var)),
                    "This test compares the same subjects in both groups. Choose the variable that identifies the subject."))
      req(box_plot_obj())
      box_plot_obj()
    })

    box_stats_df <- reactive({
      req(box_result())
      stats <- box_result()$stats
      req(stats)
      df_stats <- dplyr::bind_rows(stats)
      df_stats <- df_stats[, c(".metric", "group1", "group2", "n1", "n2",
                                "statistic", "p", "p.adj", "p.adj.signif",
                                "note")]
      names(df_stats) <- c("Metric", "Group 1", "Group 2", "n1", "n2",
                            "Statistic", "p", "p.adj", "Significance", "Note")
      num_cols <- sapply(df_stats, is.numeric)
      df_stats[num_cols] <- lapply(df_stats[num_cols], round, 4)
      df_stats
    })

    output$box_stats <- DT::renderDT({
      DT::datatable(box_stats_df(), rownames = FALSE,
                    class = "compact hover small",
                    options = list(dom = "ftp", pageLength = 10,
                                   order = list(list(7, "asc"))))
    })

    output$box_stats_download <- downloadHandler(
      filename = function() paste0("alpha_statistics_", Sys.Date(), ".csv"),
      content  = function(file) {
        req(box_result())
        out <- box_stats_df()
        out$Test <- box_result()$test_label
        write.csv(out, file, row.names = FALSE)
      }
    )

    output$box_download <- downloadHandler(
      filename = function() paste0("alpha_boxplots_", Sys.Date(), ".pdf"),
      content  = function(file) {
        req(box_plot_obj())
        ggplot2::ggsave(file, box_plot_obj(), width = 12, height = 8)
      }
    )

    # ── 2. Rarefaction curves ────────────────────────────────────────────────
    # ── 3. Diversity gradient ────────────────────────────────────────────────
    output$grad_xvar_ui <- renderUI({
      feat <- upload_data$features
      req(feat)
      cont_vars <- feat$cont_vars
      if (is.null(cont_vars) || length(cont_vars) == 0) {
        df <- tryCatch(data.frame(sample_data(ps_work())), error = function(e) NULL)
        if (!is.null(df)) {
          cont_vars <- names(df)[sapply(df, is.numeric)]
        }
      }
      if (length(cont_vars) == 0) {
        return(tags$p(class = "text-muted small",
                      "No continuous variables found in metadata."))
      }
      selectInput(ns("grad_xvar"), "X-axis variable",
                  choices = cont_vars, selected = cont_vars[1])
    })

    grad_plot_obj <- reactive({
      req(diversity_df(), input$grad_xvar, input$grad_metric)
      gv <- group_var()
      tryCatch(
        plot_diversity_gradient(diversity_df(),
                                metric        = input$grad_metric,
                                x_var         = input$grad_xvar,
                                group_var     = gv,
                                add_smooth    = input$grad_smooth,
                                smooth_method = input$grad_method),
        error = function(e) {
          message("[mod_alpha] gradient error: ", e$message)
          NULL
        }
      )
    })

    output$grad_plot <- renderPlot({
      req(grad_plot_obj())
      grad_plot_obj()
    })

    output$grad_download <- downloadHandler(
      filename = function() paste0("alpha_gradient_", Sys.Date(), ".pdf"),
      content  = function(file) {
        req(grad_plot_obj())
        ggplot2::ggsave(file, grad_plot_obj(), width = 10, height = 7)
      }
    )

    # ── 4. Diversity table ───────────────────────────────────────────────────
    output$alpha_table <- DT::renderDT({
      req(diversity_df())
      df <- diversity_df()
      num_cols <- sapply(df, is.numeric)
      df[num_cols] <- lapply(df[num_cols], round, 3)
      DT::datatable(
        df,
        options = list(pageLength = 15, scrollX = TRUE),
        rownames = FALSE,
        filter   = "top"
      )
    })

    output$tab_download <- downloadHandler(
      filename = function() paste0("alpha_diversity_", Sys.Date(), ".csv"),
      content  = function(file) {
        req(diversity_df())
        write.csv(diversity_df(), file, row.names = FALSE)
      }
    )

    # ── Enable download buttons ──────────────────────────────────────────────
    observe({
      req(identical(parent_input$top_nav, "alpha"))   # only while the tab is open
      shinyjs::toggleState("box_download",  !is.null(box_plot_obj()))
      shinyjs::toggleState("grad_download", !is.null(grad_plot_obj()))
      shinyjs::toggleState("tab_download",  !is.null(diversity_df()))
    })
  })
}

