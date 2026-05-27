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
          selectInput(ns("box_test"), "Statistical test",
                      choices  = c("Wilcoxon"       = "wilcox",
                                   "t-test"         = "t.test",
                                   "Kruskal-Wallis" = "kruskal",
                                   "ANOVA"          = "anova"),
                      selected = "wilcox"),
          checkboxInput(ns("box_points"), "Show data points", value = TRUE),
          hr(),
          downloadButton(ns("box_download"), "Download plot",
                         class = "btn btn-outline-secondary btn-sm w-100",
                         disabled = "disabled")
        ),
        plotOutput(ns("box_plot"), height = "550px"),
        hr(),
        DT::DTOutput(ns("box_stats"))
      )
    ),

    # ── 2. Rarefaction curves ────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("chart-line"), " Rarefaction"),
      layout_sidebar(
        sidebar = sidebar(
          width = 320, open = "always",
          tags$h6("Parameters"),
          numericInput(ns("rare_step"), "Step size",
                       value = 500, min = 100, max = 5000, step = 100),
          numericInput(ns("rare_n_iter"), "Iterations per step",
                       value = 3, min = 1, max = 10, step = 1),
          hr(),
          downloadButton(ns("rare_download"), "Download plot",
                         class = "btn btn-outline-secondary btn-sm w-100",
                         disabled = "disabled")
        ),
        plotOutput(ns("rare_plot"), height = "550px")
      )
    ),

    # ── 3. Diversity gradient ────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("chart-area"), " Gradient"),
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
    box_result <- reactive({
      req(diversity_df())
      gv <- group_var()
      req(gv, length(input$box_metrics) > 0)
      tryCatch({
        plot_alpha_diversity(diversity_df(),
                             metrics     = input$box_metrics,
                             group_var   = gv,
                             test        = input$box_test,
                             show_points = input$box_points)
      }, error = function(e) {
        message("[mod_alpha] boxplot error: ", e$message)
        NULL
      })
    })

    box_plot_obj <- reactive({
      req(box_result())
      box_result()$plot
    })

    output$box_plot <- renderPlot({
      req(box_plot_obj())
      box_plot_obj()
    })

    observe({
      message("[box_result debug] is null: ", is.null(box_result()))
      if (!is.null(box_result())) {
        message("[box_result debug] names: ", paste(names(box_result()), collapse=", "))
        message("[box_result debug] stats null: ", is.null(box_result()$stats))
        message("[box_result debug] stats length: ", length(box_result()$stats))
      }
    })

    output$box_stats <- DT::renderDT({
      req(box_result())
      stats <- box_result()$stats
      req(stats)
      df_stats <- dplyr::bind_rows(stats)
      df_stats <- df_stats[, c(".metric", "group1", "group2", "n1", "n2",
                                "statistic", "p", "p.adj", "p.adj.signif")]
      names(df_stats) <- c("Metric", "Group 1", "Group 2", "n1", "n2",
                            "Statistic", "p", "p.adj", "Significance")
      num_cols <- sapply(df_stats, is.numeric)
      df_stats[num_cols] <- lapply(df_stats[num_cols], round, 4)
      DT::datatable(df_stats, options = list(dom = "t", pageLength = 20),
                    rownames = FALSE)
    })

    output$box_download <- downloadHandler(
      filename = function() paste0("alpha_boxplots_", Sys.Date(), ".pdf"),
      content  = function(file) {
        req(box_plot_obj())
        ggplot2::ggsave(file, box_plot_obj(), width = 12, height = 8)
      }
    )

    # ── 2. Rarefaction curves ────────────────────────────────────────────────
    rare_plot_obj <- reactive({
      req(ps_work())
      gv <- group_var()
      tryCatch(
        plot_rarefaction_curves(ps_work(),
                                step      = input$rare_step,
                                n_iter    = input$rare_n_iter,
                                group_var = gv),
        error = function(e) {
          message("[mod_alpha] rarefaction error: ", e$message)
          NULL
        }
      )
    })

    output$rare_plot <- renderPlot({
      req(rare_plot_obj())
      rare_plot_obj()
    })

    output$rare_download <- downloadHandler(
      filename = function() paste0("rarefaction_", Sys.Date(), ".pdf"),
      content  = function(file) {
        req(rare_plot_obj())
        ggplot2::ggsave(file, rare_plot_obj(), width = 10, height = 7)
      }
    )

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
      shinyjs::toggleState("box_download",  !is.null(box_plot_obj()))
      shinyjs::toggleState("rare_download", !is.null(rare_plot_obj()))
      shinyjs::toggleState("grad_download", !is.null(grad_plot_obj()))
      shinyjs::toggleState("tab_download",  !is.null(diversity_df()))
    })
  })
}

