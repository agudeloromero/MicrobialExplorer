# =============================================================================
# mod_functional.R — Functional Prediction Analysis Module
# MicrobialExplorer
# =============================================================================
# Follows the established module pattern:
#   mod_functional_ui(id)
#   mod_functional_server(id, upload_data, parent_input)
#
# Pipeline functions sourced from R/microbiome_functional.R
# Requires PICRUSt2 output files OR falls back to simulated demo data
# =============================================================================

# --- UI -----------------------------------------------------------------------

mod_functional_ui <- function(id) {
  ns <- NS(id)

  tagList(
    # ── Header ────────────────────────────────────────────────────────────────
    fluidRow(
      column(12,
        div(class = "module-header",
          h3(icon("dna"), "Functional Prediction Analysis"),
          p("Predict the functional potential of your microbiome using PICRUSt2
            outputs. Explore KEGG pathway abundances, functional diversity, and
            differentially abundant functional features between groups.",
            class = "text-muted")
        )
      )
    ),

    # ── Input controls ────────────────────────────────────────────────────────
    fluidRow(
      column(3,
        wellPanel(
          h4(icon("sliders-h"), "Settings"),

          # Data source toggle
          div(class = "form-group",
            tags$label("Data source"),
            radioButtons(
              ns("data_source"),
              label    = NULL,
              choices  = c("Demo data (simulated)" = "demo",
                           "Upload PICRUSt2 files"  = "upload"),
              selected = "demo"
            )
          ),

          # PICRUSt2 file uploads (shown only when upload selected)
          conditionalPanel(
            condition = sprintf("input['%s'] == 'upload'", ns("data_source")),

            tags$label("KO predictions (.tsv / .tsv.gz)"),
            fileInput(ns("ko_file"), label = NULL,
                      accept = c(".tsv", ".gz")),

            tags$label("Pathway abundances (.tsv / .tsv.gz)"),
            fileInput(ns("pathway_file"), label = NULL,
                      accept = c(".tsv", ".gz")),

            tags$label("EC numbers (.tsv / .tsv.gz) — optional"),
            fileInput(ns("ec_file"), label = NULL,
                      accept = c(".tsv", ".gz")),

            tags$label("NSTI values (.tsv) — optional"),
            fileInput(ns("nsti_file"), label = NULL,
                      accept = ".tsv")
          ),

          hr(),

          # Group variable — mirrors other modules
          uiOutput(ns("group_var_ui")),

          # Analysis parameters
          numericInput(
            ns("top_n_pathways"),
            label = "Top pathways to display",
            value = 25, min = 5, max = 60, step = 5
          ),

          numericInput(
            ns("top_n_heatmap"),
            label = "Heatmap features",
            value = 40, min = 10, max = 80, step = 5
          ),

          uiOutput(ns("heat_sortby_ui")),

          numericInput(
            ns("da_alpha"),
            label = "DA significance threshold (q)",
            value = 0.05, min = 0.001, max = 0.2, step = 0.005
          ),

          numericInput(
            ns("top_n_da"),
            label = "Top DA features to plot",
            value = 25, min = 5, max = 50, step = 5
          ),

          hr(),

          actionButton(
            ns("run_analysis"),
            label = "Run analysis",
            icon  = icon("play"),
            class = "btn-primary btn-block"
          )
        )
      ),

      # ── Main output area ────────────────────────────────────────────────────
      column(9,

        # Status / progress
        uiOutput(ns("status_ui")),

        # Tabbed results
        tabsetPanel(
          id = ns("result_tabs"),

          # Tab 1: NSTI quality
          tabPanel(
            title = tagList(icon("check-circle"), "NSTI Quality"),
            value = "nsti",
            br(),
            uiOutput(ns("nsti_info")),
            withSpinner(
              plotOutput(ns("plot_nsti"), height = "500px"),
              type = 4, color = "#3498db"
            ),
            br(),
            downloadButton(ns("dl_nsti"), "Download plot", class = "btn-sm")
          ),

          # Tab 2: Pathway abundance
          tabPanel(
            title = tagList(icon("chart-bar"), "Pathway Abundance"),
            value = "pathways",
            br(),
            withSpinner(
              plotOutput(ns("plot_pathways"), height = "600px"),
              type = 4, color = "#3498db"
            ),
            br(),
            downloadButton(ns("dl_pathways"), "Download plot", class = "btn-sm"),
            downloadButton(ns("dl_pathway_table"), "Download table",
                           class = "btn-sm ml-2"),
            br(), br(),
            DT::dataTableOutput(ns("pathway_table"))
          ),

          # Tab 3: KEGG categories
          tabPanel(
            title = tagList(icon("layer-group"), "KEGG Categories"),
            value = "kegg",
            br(),
            withSpinner(
              plotOutput(ns("plot_kegg"), height = "600px"),
              type = 4, color = "#3498db"
            ),
            br(),
            downloadButton(ns("dl_kegg"), "Download plot", class = "btn-sm"),
            br(), br(),
            DT::dataTableOutput(ns("kegg_table"))
          ),

          # Tab 4: Functional diversity
          tabPanel(
            title = tagList(icon("project-diagram"), "Functional Diversity"),
            value = "diversity",
            br(),
            h5("Alpha diversity"),
            withSpinner(
              plotOutput(ns("plot_func_alpha"), height = "320px"),
              type = 4, color = "#3498db"
            ),
            br(),
            h5("Beta diversity (PCoA)"),
            withSpinner(
              plotOutput(ns("plot_func_beta"), height = "420px"),
              type = 4, color = "#3498db"
            ),
            br(),
            fluidRow(
              column(6,
                downloadButton(ns("dl_func_alpha"), "Download alpha plot",
                               class = "btn-sm")),
              column(6,
                downloadButton(ns("dl_func_beta"), "Download beta plot",
                               class = "btn-sm"))
            )
          ),

          # Tab 5: Pathway heatmap
          tabPanel(
            title = tagList(icon("th"), "Heatmap"),
            value = "heatmap",
            br(),
            withSpinner(
              plotOutput(ns("plot_heatmap"), height = "700px"),
              type = 4, color = "#3498db"
            ),
            br(),
            downloadButton(ns("dl_heatmap"), "Download plot", class = "btn-sm")
          ),

          # Tab 6: Differential abundance
          tabPanel(
            title = tagList(icon("not-equal"), "Differential Abundance"),
            value = "da",
            br(),

            tabsetPanel(
              tabPanel("Pathway DA",
                br(),
                uiOutput(ns("pathway_da_summary")),

                uiOutput(ns("pathway_da_bar_ui")),



                withSpinner(
                  plotOutput(ns("plot_pathway_da_bubble"), height = "420px"),
                  type = 4, color = "#3498db"
                ),
                br(),
                fluidRow(
                  column(4, downloadButton(ns("dl_pathway_da_bar"),
                                          "Bar chart", class = "btn-sm")),
                  column(4, downloadButton(ns("dl_pathway_da_bubble"),
                                          "Bubble plot", class = "btn-sm")),
                  column(4, downloadButton(ns("dl_pathway_da_csv"),
                                          "Results CSV", class = "btn-sm"))
                ),
                br(),
                DT::dataTableOutput(ns("pathway_da_table"))
              ),

              tabPanel("KO DA",
                br(),
                uiOutput(ns("ko_da_summary")),

                uiOutput(ns("ko_da_bar_ui")),



                withSpinner(
                  plotOutput(ns("plot_ko_da_bubble"), height = "420px"),
                  type = 4, color = "#3498db"
                ),
                br(),
                fluidRow(
                  column(4, downloadButton(ns("dl_ko_da_bar"),
                                          "Bar chart", class = "btn-sm")),
                  column(4, downloadButton(ns("dl_ko_da_bubble"),
                                          "Bubble plot", class = "btn-sm")),
                  column(4, downloadButton(ns("dl_ko_da_csv"),
                                          "Results CSV", class = "btn-sm"))
                ),
                br(),
                DT::dataTableOutput(ns("ko_da_table"))
              )
            )
          )
        ) # end tabsetPanel
      )
    )
  )
}


# --- Server -------------------------------------------------------------------

mod_functional_server <- function(id, upload_data, parent_input) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── 1. Group variable UI (mirrors other modules) ─────────────────────────

    output$group_var_ui <- renderUI({
      req(upload_data$ps)
      ps <- upload_data$ps
      req(!is.null(ps))
      meta_vars <- colnames(sample_data(ps))
      cat_vars  <- meta_vars[sapply(meta_vars, function(v) {
        x <- sample_data(ps)[[v]]
        is.factor(x) || is.character(x) ||
          (is.numeric(x) && length(unique(x)) <= 10)
      })]
      selectInput(
        ns("group_variable"),
        label   = "Grouping variable",
        choices = cat_vars,
        selected = if ("group" %in% cat_vars) "group" else cat_vars[1]
      )
    })

    # Reactive: resolve group variable (own selector or parent)
    group_var <- reactive({
      if (!is.null(input$group_variable) && nchar(input$group_variable) > 0) {
        input$group_variable
      } else {
        parent_input$group_variable
      }
    })


    # ── 2. Build functional data from phyloseq or uploaded files ──────────────

    func_data <- reactive({
      req(upload_data$ps)

      ps       <- upload_data$ps
      meta_df  <- as.data.frame(as.matrix(sample_data(ps)))

      if (input$data_source == "demo") {
        # ── Demo path: simulate KO and pathway matrices from the phyloseq ──────
        set.seed(42)
        n_samples  <- nsamples(ps)
        sample_ids <- sample_names(ps)

        n_pathways <- 150
        n_ko       <- 300

        pathway_ids <- paste0("PWY-", 1000 + seq_len(n_pathways))
        ko_ids      <- paste0("K", formatC(seq_len(n_ko), width = 5,
                                           flag = "0"))

        # Group-structured signal for demo realism
        grp     <- as.character(meta_df[[group_var()]])
        grp_lvl <- unique(grp)
        is_grp1 <- grp == grp_lvl[1]

        pathway_mat <- matrix(
          abs(rnorm(n_pathways * n_samples, mean = 500, sd = 80)),
          nrow = n_pathways, ncol = n_samples,
          dimnames = list(pathway_ids, sample_ids)
        )
        ko_mat <- matrix(
          abs(rnorm(n_ko * n_samples, mean = 300, sd = 60)),
          nrow = n_ko, ncol = n_samples,
          dimnames = list(ko_ids, sample_ids)
        )
        n_sig_pw <- round(n_pathways * 0.20)
        n_sig_ko <- round(n_ko * 0.20)
        pathway_mat[seq_len(n_sig_pw), is_grp1] <- pathway_mat[seq_len(n_sig_pw), is_grp1] * 3
        ko_mat[seq_len(n_sig_ko), is_grp1] <- ko_mat[seq_len(n_sig_ko), is_grp1] * 3
        list(
          pathway  = pathway_mat,
          ko       = ko_mat,
          metadata = meta_df,
          source   = "demo"
        )

      } else {
        # ── Upload path: parse user-provided TSV files ──────────────────────

        read_tsv_matrix <- function(file_info) {
          if (is.null(file_info)) return(NULL)
          tryCatch({
            df <- read.table(file_info$datapath, sep = "\t", header = TRUE,
                             row.names = 1, check.names = FALSE,
                             stringsAsFactors = FALSE)
            if ("description" %in% colnames(df)) {
              df <- df[, colnames(df) != "description", drop = FALSE]
            }
            as.matrix(df)
          }, error = function(e) {
            showNotification(
              paste0("Error reading ", file_info$name, ": ", e$message),
              type = "error", duration = 8
            )
            NULL
          })
        }

        pathway_mat <- read_tsv_matrix(input$pathway_file)
        ko_mat      <- read_tsv_matrix(input$ko_file)
        nsti_mat    <- read_tsv_matrix(input$nsti_file)

        # Align samples to metadata
        align_mat <- function(mat, meta) {
          if (is.null(mat)) return(NULL)
          shared <- intersect(colnames(mat), rownames(meta))
          if (length(shared) == 0) {
            showNotification(
              "No shared samples between uploaded file and metadata.",
              type = "error", duration = 8
            )
            return(NULL)
          }
          mat[, shared, drop = FALSE]
        }

        list(
          pathway  = align_mat(pathway_mat, meta_df),
          ko       = align_mat(ko_mat, meta_df),
          nsti     = nsti_mat,
          metadata = meta_df,
          source   = "upload"
        )
      }
    })


    # ── 3. Status banner ──────────────────────────────────────────────────────

    output$heat_sortby_ui <- renderUI({
      feat <- upload_data$features
      gvars <- if (!is.null(feat) && length(feat$group_vars) > 0) feat$group_vars else c("disease_status")
      choices <- c("Hierarchical clustering" = "none", setNames(gvars, gvars))
      selectInput(ns("heat_sort_group"), "Sort columns by",
                  choices = choices, selected = "none")
    })



    output$status_ui <- renderUI({
      req(upload_data$ps)
      fd <- func_data()
      req(!is.null(fd))

      has_pathway <- !is.null(fd$pathway)
      has_ko      <- !is.null(fd$ko)

      tags$div(
        class = if (fd$source == "demo") "alert alert-info" else "alert alert-success",
        icon(if (fd$source == "demo") "info-circle" else "check-circle"),
        if (fd$source == "demo") {
          paste0(" Using simulated demo data. ",
                 "Upload real PICRUSt2 files to analyse your own data.")
        } else {
          paste0(" Uploaded data loaded. Pathways: ",
                 if (has_pathway) nrow(fd$pathway) else "not provided",
                 " | KOs: ",
                 if (has_ko) nrow(fd$ko) else "not provided",
                 " | Samples: ", ncol(fd$pathway %||% fd$ko))
        }
      )
    })


    # ── 4. Run analysis on button click ───────────────────────────────────────

    results <- eventReactive(input$run_analysis, {
      req(func_data(), group_var())

      fd  <- func_data()
      gv  <- group_var()
      out <- list()

      da_alpha_val  <- input$da_alpha
      top_n_da_val  <- input$top_n_da
      withProgress(message = "Running functional analysis…", value = 0, {

        # NSTI quality
        incProgress(0.1, detail = "NSTI quality assessment")
        tryCatch({
          p2_list <- list(
            nsti     = fd$nsti,
            ko       = fd$ko,
            pathway  = fd$pathway,
            metadata = fd$metadata
          )
          out$nsti <- plot_nsti_quality(p2_list, group_var = gv)
        }, error = function(e) {
          showNotification(paste0("NSTI: ", e$message),
                           type = "warning", duration = 6)
        })

        # Pathway abundance
        if (!is.null(fd$pathway)) {
          incProgress(0.15, detail = "Pathway abundance")
          tryCatch({
            out$pathways <- plot_pathway_abundance(
              fd$pathway, fd$metadata,
              group_var = gv, top_n = input$top_n_pathways
            )
          }, error = function(e)
            showNotification(paste0("Pathway abundance: ", e$message),
                             type = "warning", duration = 6))

          incProgress(0.1, detail = "Pathway heatmap")
          tryCatch({
            sort_val <- if (!is.null(input$heat_sort_group)) input$heat_sort_group else "none"
            gv_heat  <- if (sort_val != "none") sort_val else NULL
            out$heatmap <- plot_functional_heatmap(
              fd$pathway, fd$metadata,
              group_var    = gv_heat,
              top_n        = input$top_n_heatmap,
              feature_type = "Pathway",
              cluster_cols = is.null(gv_heat)
            )

          }, error = function(e)
            showNotification(paste0("Heatmap: ", e$message), type = "warning", duration = 6))
          tryCatch({
            out$pathway_da <- test_functional_da(
              fd$pathway, fd$metadata,
              group_var    = gv,
              feature_type = "Pathway",
              alpha        = da_alpha_val,
              top_n_plot   = top_n_da_val
            )
          }, error = function(e)
            showNotification(paste0("Pathway DA: ", e$message),
                             type = "warning", duration = 6))
        }

        # KO / KEGG analyses
        if (!is.null(fd$ko)) {
          incProgress(0.1, detail = "KEGG categories")
          tryCatch({
            out$kegg <- summarise_kegg_categories(
              fd$ko, fd$metadata, group_var = gv
            )
          }, error = function(e)
            showNotification(paste0("KEGG categories: ", e$message),
                             type = "warning", duration = 6))

          incProgress(0.15, detail = "Functional diversity")
          tryCatch({
            out$func_div <- analyse_functional_diversity(
              fd$ko, fd$metadata,
              group_var    = gv,
              feature_type = "KO"
            )
          }, error = function(e)
            showNotification(paste0("Functional diversity: ", e$message),
                             type = "warning", duration = 6))

          incProgress(0.15, detail = "KO differential abundance")
          tryCatch({
            out$ko_da <- test_functional_da(
              fd$ko, fd$metadata,
              group_var    = gv,
              feature_type = "KO",
              alpha        = da_alpha_val,
              top_n_plot   = top_n_da_val
            )
          }, error = function(e)
            showNotification(paste0("KO DA: ", e$message),
                             type = "warning", duration = 6))
        }

        incProgress(0.1, detail = "Done")
      })

      out
    })


    # ── 5. NSTI plots & info ──────────────────────────────────────────────────

    output$nsti_info <- renderUI({
      r <- results()
      req(!is.null(r$nsti))
      med <- round(r$nsti$median_nsti, 4)
      quality_class <- if (med <= 0.06) "success" else
                       if (med <= 0.15) "warning" else "danger"
      tags$div(
        class = paste0("alert alert-", quality_class),
        strong("Median NSTI: "), med, " — ",
        if (med <= 0.06) "High prediction quality." else
        if (med <= 0.15) "Moderate prediction quality." else
          "Low prediction quality: some samples exceed the recommended threshold."
      )
    })

    output$plot_nsti <- renderPlot({
      req(results()$nsti)
      results()$nsti$plot
    })

    output$dl_nsti <- downloadHandler(
      filename = "nsti_quality.pdf",
      content  = function(file) {
        req(results()$nsti)
        ggsave(file, results()$nsti$plot, width = 12, height = 10)
      }
    )


    # ── 6. Pathway abundance ──────────────────────────────────────────────────

    output$plot_pathways <- renderPlot({
      req(results()$pathways)
      results()$pathways$plot
    })

    output$pathway_table <- DT::renderDataTable({
      req(results()$pathways)
      results()$pathways$summary %>%
        DT::datatable(
          rownames  = FALSE,
          filter    = "top",
          options   = list(pageLength = 15, scrollX = TRUE)
        ) %>%
        DT::formatRound(columns = c("mean_abund", "se_abund"), digits = 3)
    })

    output$dl_pathways <- downloadHandler(
      filename = "pathway_abundance.pdf",
      content  = function(file) {
        req(results()$pathways)
        ggsave(file, results()$pathways$plot, width = 12, height = 12)
      }
    )

    output$dl_pathway_table <- downloadHandler(
      filename = "pathway_abundance.csv",
      content  = function(file) {
        req(results()$pathways)
        write.csv(results()$pathways$summary, file, row.names = FALSE)
      }
    )


    # ── 7. KEGG categories ────────────────────────────────────────────────────

    output$plot_kegg <- renderPlot({
      req(results()$kegg)
      results()$kegg$plot
    })

    output$kegg_table <- DT::renderDataTable({
      req(results()$kegg)
      results()$kegg$summary %>%
        DT::datatable(
          rownames = FALSE,
          filter   = "top",
          options  = list(pageLength = 15, scrollX = TRUE)
        ) %>%
        DT::formatRound(columns = c("mean_pct", "se_pct"), digits = 2)
    })

    output$dl_kegg <- downloadHandler(
      filename = "kegg_categories.pdf",
      content  = function(file) {
        req(results()$kegg)
        ggsave(file, results()$kegg$plot, width = 14, height = 10)
      }
    )


    # ── 8. Functional diversity ───────────────────────────────────────────────

    output$plot_func_alpha <- renderPlot({
      req(results()$func_div)
      results()$func_div$p_alpha
    })

    output$plot_func_beta <- renderPlot({
      req(results()$func_div)
      results()$func_div$p_beta
    })

    output$dl_func_alpha <- downloadHandler(
      filename = "functional_alpha_diversity.pdf",
      content  = function(file) {
        req(results()$func_div)
        ggsave(file, results()$func_div$p_alpha, width = 12, height = 6)
      }
    )

    output$dl_func_beta <- downloadHandler(
      filename = "functional_beta_diversity.pdf",
      content  = function(file) {
        req(results()$func_div)
        ggsave(file, results()$func_div$p_beta, width = 9, height = 8)
      }
    )


    # ── 9. Functional heatmap ─────────────────────────────────────────────────

    output$plot_heatmap <- renderPlot({
      req(results(), func_data())
      sort_val   <- if (!is.null(input$heat_sort_group)) input$heat_sort_group else "none"
      gv_heat    <- if (sort_val != "none") sort_val else NULL
      fd         <- func_data()
      plot_functional_heatmap(
        fd$pathway, fd$metadata,
        group_var    = gv_heat,
        top_n        = input$top_n_heatmap,
        feature_type = "Pathway",
        cluster_cols = is.null(gv_heat)
      )
    })





    # ── 10. Pathway DA ────────────────────────────────────────────────────────

    output$pathway_da_summary <- renderUI({
      r <- results()$pathway_da
      req(!is.null(r))
      n    <- nrow(r$results)
      nsig <- r$n_sig
      tags$div(
        class = if (nsig > 0) "alert alert-success" else "alert alert-info",
        icon(if (nsig > 0) "check" else "info-circle"),
        sprintf(" %d pathway%s tested | %d significantly different (q < %g)",
                n, if (n == 1) "" else "s",
                nsig, input$da_alpha)
      )
    })

    output$pathway_da_bar_ui <- renderUI({
      req(results()$pathway_da)
      if (is.null(results()$pathway_da$p_bar)) {
        tags$p(class = "text-muted", icon("info-circle"), " No significantly different pathways at this threshold.")
      } else { plotOutput(ns("plot_pathway_da_bar"), height = "600px") }
    })
    output$plot_pathway_da_bar <- renderPlot({
      req(results()$pathway_da, results()$pathway_da$p_bar)
      results()$pathway_da$p_bar
    })

    output$plot_pathway_da_bubble <- renderPlot({
      req(results()$pathway_da)
      tryCatch(results()$pathway_da$p_bubble,
        error = function(e) ggplot() + annotate('text', x=0.5, y=0.5, label='Plot unavailable', size=5, colour='grey50') + theme_void())
    })


    output$pathway_da_table <- DT::renderDataTable({
      req(results()$pathway_da)
      results()$pathway_da$results %>%
        mutate(across(where(is.numeric), ~ round(., 5))) %>%
        DT::datatable(
          rownames = FALSE,
          filter   = "top",
          options  = list(pageLength = 15, scrollX = TRUE)
        ) %>%
        DT::formatStyle(
          "diff_abund",
          backgroundColor = DT::styleEqual(TRUE, "#d5f5e3")
        )
    })

    output$dl_pathway_da_bar <- downloadHandler(
      filename = "pathway_da_bar.pdf",
      content  = function(file) {
        req(results()$pathway_da$p_bar)
        ggsave(file, results()$pathway_da$p_bar, width = 12, height = 10)
      }
    )

    output$dl_pathway_da_bubble <- downloadHandler(
      filename = "pathway_da_bubble.pdf",
      content  = function(file) {
        req(results()$pathway_da)
        ggsave(file, results()$pathway_da$p_bubble, width = 10, height = 8)
      }
    )

    output$dl_pathway_da_csv <- downloadHandler(
      filename = "pathway_da_results.csv",
      content  = function(file) {
        req(results()$pathway_da)
        write.csv(results()$pathway_da$results, file, row.names = FALSE)
      }
    )


    # ── 11. KO DA ─────────────────────────────────────────────────────────────

    output$ko_da_summary <- renderUI({
      r <- results()$ko_da
      req(!is.null(r))
      n    <- nrow(r$results)
      nsig <- r$n_sig
      tags$div(
        class = if (nsig > 0) "alert alert-success" else "alert alert-info",
        icon(if (nsig > 0) "check" else "info-circle"),
        sprintf(" %d KO%s tested | %d significantly different (q < %g)",
                n, if (n == 1) "" else "s",
                nsig, input$da_alpha)
      )
    })

    output$ko_da_bar_ui <- renderUI({
      req(results()$ko_da)
      if (is.null(results()$ko_da$p_bar)) {
        tags$p(class = "text-muted", icon("info-circle"), " No significantly different KOs at this threshold.")
      } else { plotOutput(ns("plot_ko_da_bar"), height = "600px") }
    })
    output$plot_ko_da_bar <- renderPlot({
      req(results()$ko_da, results()$ko_da$p_bar)
      results()$ko_da$p_bar
    })


    output$plot_ko_da_bubble <- renderPlot({
      req(results()$ko_da)
      tryCatch(results()$ko_da$p_bubble,
        error = function(e) ggplot() + annotate('text', x=0.5, y=0.5, label='Plot unavailable', size=5, colour='grey50') + theme_void())
    })

    output$ko_da_table <- DT::renderDataTable({
      req(results()$ko_da)
      results()$ko_da$results %>%
        mutate(across(where(is.numeric), ~ round(., 5))) %>%
        DT::datatable(
          rownames = FALSE,
          filter   = "top",
          options  = list(pageLength = 15, scrollX = TRUE)
        ) %>%
        DT::formatStyle(
          "diff_abund",
          backgroundColor = DT::styleEqual(TRUE, "#d5f5e3")
        )
    })

    output$dl_ko_da_bar <- downloadHandler(
      filename = "ko_da_bar.pdf",
      content  = function(file) {
        req(results()$ko_da$p_bar)
        ggsave(file, results()$ko_da$p_bar, width = 12, height = 10)
      }
    )

    output$dl_ko_da_bubble <- downloadHandler(
      filename = "ko_da_bubble.pdf",
      content  = function(file) {
        req(results()$ko_da)
        ggsave(file, results()$ko_da$p_bubble, width = 10, height = 8)
      }
    )

    output$dl_ko_da_csv <- downloadHandler(
      filename = "ko_da_results.csv",
      content  = function(file) {
        req(results()$ko_da)
        write.csv(results()$ko_da$results, file, row.names = FALSE)
      }
    )

  }) # end moduleServer
}


# =============================================================================
# NULL-coalescing helper (base R safe)
# =============================================================================
`%||%` <- function(a, b) if (!is.null(a)) a else b
