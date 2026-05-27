# =============================================================================
# mod_beta.R — Beta Diversity Module
# =============================================================================
# Sub-tabs:
#   1. Ordination  — PCoA / NMDS with ellipses
#   2. PERMANOVA   — group difference test + pairwise table
#   3. Dispersion  — betadisper boxplot + significance
#   4. Heatmap     — hierarchically clustered distance matrix
#   5. Envfit      — metadata vectors on ordination
# =============================================================================

# ── UI ────────────────────────────────────────────────────────────────────────
mod_beta_ui <- function(id) {
  ns <- NS(id)

  navset_card_tab(
    id = ns("beta_tabs"),

    # ── 1. Ordination ────────────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("circle-dot"), " Ordination"),
      layout_sidebar(
        sidebar = sidebar(
          width = 320, open = "always",
          tags$h6("Parameters"),
          selectInput(ns("ord_distance"), "Distance metric",
                      choices  = c("Bray-Curtis"   = "bray",
                                   "Jaccard"        = "jaccard",
                                   "Aitchison"      = "aitchison",
                                   "Unweighted UniFrac" = "unifrac",
                                   "Weighted UniFrac"   = "wunifrac"),
                      selected = "bray"),
          selectInput(ns("ord_method"), "Ordination method",
                      choices  = c("PCoA" = "PCoA", "NMDS" = "NMDS"),
                      selected = "PCoA"),
          checkboxInput(ns("ord_ellipse"), "Draw 95% ellipses", value = TRUE),
          checkboxInput(ns("ord_spider"),  "Draw spider lines",  value = FALSE),
          uiOutput(ns("ord_shape_ui")),
          hr(),
          downloadButton(ns("ord_download"), "Download plot",
                         class = "btn btn-outline-secondary btn-sm w-100",
                         disabled = "disabled")
        ),
        plotOutput(ns("ord_plot"), height = "580px")
      )
    ),

    # ── 2. PERMANOVA ────────────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("vials"), " PERMANOVA"),
      layout_sidebar(
        sidebar = sidebar(
          width = 320, open = "always",
          tags$h6("Parameters"),
          selectInput(ns("perm_distance"), "Distance metric",
                      choices  = c("Bray-Curtis"   = "bray",
                                   "Jaccard"        = "jaccard",
                                   "Aitchison"      = "aitchison",
                                   "Unweighted UniFrac" = "unifrac",
                                   "Weighted UniFrac"   = "wunifrac"),
                      selected = "bray"),
          uiOutput(ns("perm_formula_ui")),
          numericInput(ns("perm_perms"), "Permutations",
                       value = 999, min = 99, max = 9999, step = 100),
          hr(),
          downloadButton(ns("perm_download"), "Download results",
                         class = "btn btn-outline-secondary btn-sm w-100",
                         disabled = "disabled")
        ),
        tags$h6("PERMANOVA result"),
        verbatimTextOutput(ns("perm_result")),
        hr(),
        tags$h6("Pairwise PERMANOVA"),
        DT::DTOutput(ns("perm_pairwise"))
      )
    ),

    # ── 3. Dispersion ───────────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("arrows-left-right"), " Dispersion"),
      layout_sidebar(
        sidebar = sidebar(
          width = 320, open = "always",
          tags$h6("Parameters"),
          selectInput(ns("disp_distance"), "Distance metric",
                      choices  = c("Bray-Curtis"   = "bray",
                                   "Jaccard"        = "jaccard",
                                   "Aitchison"      = "aitchison",
                                   "Unweighted UniFrac" = "unifrac",
                                   "Weighted UniFrac"   = "wunifrac"),
                      selected = "bray"),
          hr(),
          downloadButton(ns("disp_download"), "Download plot",
                         class = "btn btn-outline-secondary btn-sm w-100",
                         disabled = "disabled")
        ),
        plotOutput(ns("disp_plot"), height = "500px")
      )
    ),

    # ── 4. Heatmap ──────────────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("border-all"), " Heatmap"),
      layout_sidebar(
        sidebar = sidebar(
          width = 320, open = "always",
          tags$h6("Parameters"),
          selectInput(ns("heat_distance"), "Distance metric",
                      choices  = c("Bray-Curtis"   = "bray",
                                   "Jaccard"        = "jaccard",
                                   "Aitchison"      = "aitchison",
                                   "Unweighted UniFrac" = "unifrac",
                                   "Weighted UniFrac"   = "wunifrac"),
                      selected = "bray"),
          selectInput(ns("heat_clust_method"), "Clustering method",
                      choices  = c("Ward.D2" = "ward.D2", "Complete" = "complete", "Average" = "average", "Single" = "single"), selected = "ward.D2"),
          hr(),
          downloadButton(ns("heat_download"), "Download plot",
                         class = "btn btn-outline-secondary btn-sm w-100",
                         disabled = "disabled")
        ),
        plotOutput(ns("heat_plot"), height = "600px")
      )
    ),

    # ── 5. Envfit ───────────────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("arrow-up-right-dots"), " Envfit"),
      layout_sidebar(
        sidebar = sidebar(
          width = 320, open = "always",
          tags$h6("Parameters"),
          selectInput(ns("env_distance"), "Distance metric",
                      choices  = c("Bray-Curtis"   = "bray",
                                   "Jaccard"        = "jaccard",
                                   "Aitchison"      = "aitchison",
                                   "Unweighted UniFrac" = "unifrac",
                                   "Weighted UniFrac"   = "wunifrac"),
                      selected = "bray"),
          numericInput(ns("env_pthresh"), "p-value threshold",
                       value = 0.05, min = 0.01, max = 0.2, step = 0.01),
          numericInput(ns("env_perms"), "Permutations",
                       value = 999, min = 99, max = 9999, step = 100),
          hr(),
          downloadButton(ns("env_download"), "Download plot",
                         class = "btn btn-outline-secondary btn-sm w-100",
                         disabled = "disabled")
        ),
        plotOutput(ns("env_plot"), height = "580px"),
        hr(),
        tags$h6("Envfit vectors"),
        DT::DTOutput(ns("env_table"))
      )
    )
  )
}


# ── Server ────────────────────────────────────────────────────────────────────
mod_beta_server <- function(id, upload_data, parent_input) {
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

    # ── Dynamic UI: shape variable selector ─────────────────────────────────
    output$ord_shape_ui <- renderUI({
      feat <- upload_data$features
      req(feat)
      vars <- c("None" = "", feat$group_vars)
      selectInput(ns("ord_shape"), "Shape variable (optional)",
                  choices = vars, selected = "")
    })

    # ── Dynamic UI: PERMANOVA formula ───────────────────────────────────────
    output$perm_formula_ui <- renderUI({
      gv <- group_var()
      req(gv)
      textInput(ns("perm_formula"), "Formula (right-hand side)",
                value = gv,
                placeholder = "e.g. disease_status + age + sex")
    })

    # ── Helper: compute distance (cached by method) ──────────────────────────
    get_distance <- function(method) {
      req(ps_work())
      tryCatch({
        dists <- compute_distances(ps_work(),
                                   methods     = method,
                                   rarefaction = TRUE)
        dists[[method]]
      }, error = function(e) {
        message("[mod_beta] distance error (", method, "): ", e$message)
        NULL
      })
    }

    # ── 1. Ordination ────────────────────────────────────────────────────────
    ord_result <- reactive({
      req(ps_work(), input$ord_distance, input$ord_method)
      d <- get_distance(input$ord_distance)
      req(d)
      tryCatch(
        run_ordination(d, method = input$ord_method, k = 3),
        error = function(e) {
          message("[mod_beta] ordination error: ", e$message)
          NULL
        }
      )
    })

    ord_plot_obj <- reactive({
      req(ord_result())
      gv <- group_var()
      req(gv)
      shape_v <- input$ord_shape
      if (is.null(shape_v) || shape_v == "") shape_v <- NULL
      tryCatch(
        plot_ordination(ord_result(), ps_work(),
                        group_var  = gv,
                        shape_var  = shape_v,
                        ellipse    = input$ord_ellipse,
                        spider     = input$ord_spider),
        error = function(e) {
          message("[mod_beta] ord plot error: ", e$message)
          NULL
        }
      )
    })

    output$ord_plot <- renderPlot({
      req(ord_plot_obj())
      ord_plot_obj()
    })

    output$ord_download <- downloadHandler(
      filename = function() paste0("ordination_", Sys.Date(), ".pdf"),
      content  = function(file) {
        req(ord_plot_obj())
        ggplot2::ggsave(file, ord_plot_obj(), width = 9, height = 8)
      }
    )

    # ── 2. PERMANOVA ────────────────────────────────────────────────────────
    perm_result <- reactive({
      req(ps_work(), input$perm_distance, input$perm_formula)
      gv <- group_var()
      req(gv)
      d <- get_distance(input$perm_distance)
      req(d)
      tryCatch(
        run_permanova(d, ps_work(),
                      formula_rhs  = input$perm_formula,
                      permutations = input$perm_perms),
        error = function(e) {
          message("[mod_beta] permanova error: ", e$message)
          NULL
        }
      )
    })

    output$perm_result <- renderPrint({
      req(perm_result())
      print(perm_result()$permanova)
    })

    output$perm_pairwise <- DT::renderDT({
      req(perm_result())
      pw <- perm_result()$pairwise
      if (is.null(pw)) return(DT::datatable(data.frame(Note = "Only 2 groups — no pairwise needed"),
                                             options = list(dom = "t"), rownames = FALSE))
      num_cols <- sapply(pw, is.numeric)
      pw[num_cols] <- lapply(pw[num_cols], round, 4)
      DT::datatable(pw, options = list(dom = "t", pageLength = 20), rownames = FALSE)
    })

    output$perm_download <- downloadHandler(
      filename = function() paste0("permanova_", Sys.Date(), ".csv"),
      content  = function(file) {
        req(perm_result())
        pw <- perm_result()$pairwise
        if (!is.null(pw)) write.csv(pw, file, row.names = FALSE)
        else write.csv(as.data.frame(print(perm_result()$permanova)), file)
      }
    )

    # ── 3. Dispersion ───────────────────────────────────────────────────────
    disp_plot_obj <- reactive({
      req(ps_work(), input$disp_distance)
      gv <- group_var()
      req(gv)
      d <- get_distance(input$disp_distance)
      req(d)
      tryCatch(
        plot_betadisper(d, ps_work(), group_var = gv),
        error = function(e) {
          message("[mod_beta] betadisper error: ", e$message)
          NULL
        }
      )
    })

    output$disp_plot <- renderPlot({
      req(disp_plot_obj())
      disp_plot_obj()
    })

    output$disp_download <- downloadHandler(
      filename = function() paste0("betadisper_", Sys.Date(), ".pdf"),
      content  = function(file) {
        req(disp_plot_obj())
        ggplot2::ggsave(file, disp_plot_obj(), width = 12, height = 6)
      }
    )


    heat_plot_obj <- reactive({
      req(ps_work(), input$heat_distance)
      gv      <- group_var()
      req(gv)
      d       <- get_distance(input$heat_distance)
      req(d)
      clust_method <- input$heat_clust_method

      tryCatch({
        dist_mat <- as.matrix(d)
        meta_df  <- data.frame(sample_data(ps_work()))

        hc        <- hclust(as.dist(dist_mat), method = clust_method)
        col_order <- rownames(dist_mat)[hc$order]

        dist_mat <- dist_mat[col_order, col_order]

        heat_df <- as.data.frame(dist_mat) %>%
          tibble::rownames_to_column("sample1") %>%
          tidyr::pivot_longer(-sample1, names_to = "sample2", values_to = "distance") %>%
          dplyr::mutate(
            sample1 = factor(sample1, levels = col_order),
            sample2 = factor(sample2, levels = rev(col_order))
          )

        ggplot2::ggplot(heat_df,
                        ggplot2::aes(x = sample1, y = sample2, fill = distance)) +
          ggplot2::geom_tile() +
          ggplot2::scale_fill_distiller(palette = "RdYlBu", direction = -1,
                                        name = "Distance") +
          ggplot2::scale_x_discrete(expand = c(0, 0)) +
          ggplot2::scale_y_discrete(expand = c(0, 0)) +
          ggplot2::labs(
            title    = "Pairwise distance matrix",
            subtitle = paste0("Sort: ",
                              ifelse(is.null(input$heat_clust_method), "ward.D2", input$heat_clust_method),

                              " | ", nrow(dist_mat), " samples"),
            x = NULL, y = NULL
          ) +
          theme_microbiome() +
          ggplot2::theme(
            axis.text.x  = ggplot2::element_text(angle = 90, hjust = 1, size = 6),
            axis.text.y  = ggplot2::element_text(size = 6),
            panel.border = ggplot2::element_rect(colour = "grey80", fill = NA)
          )
      }, error = function(e) {
        message("[mod_beta] heatmap error: ", e$message)
        NULL
      })
    })
    output$heat_plot <- renderPlot({
      req(heat_plot_obj())
      heat_plot_obj()
    })

    output$heat_download <- downloadHandler(
      filename = function() paste0("distance_heatmap_", Sys.Date(), ".pdf"),
      content  = function(file) {
        req(heat_plot_obj())
        ggplot2::ggsave(file, heat_plot_obj(), width = 12, height = 10)
      }
    )

    # ── 5. Envfit ───────────────────────────────────────────────────────────
    env_result <- reactive({
      req(ps_work(), input$env_distance, ord_result())
      gv <- group_var()
      req(gv)
      # Recompute ordination with envfit distance if different
      d <- get_distance(input$env_distance)
      req(d)
      ord <- tryCatch(run_ordination(d, method = "PCoA", k = 3),
                      error = function(e) NULL)
      req(ord)
      tryCatch(
        plot_envfit(ord, ps_work(),
                    group_var    = gv,
                    p_threshold  = input$env_pthresh,
                    permutations = input$env_perms),
        error = function(e) {
          message("[mod_beta] envfit error: ", e$message)
          NULL
        }
      )
    })

    output$env_plot <- renderPlot({
      req(env_result())
      env_result()$plot
    })

    output$env_table <- DT::renderDT({
      req(env_result())
      vecs <- env_result()$vectors
      req(vecs)
      num_cols <- sapply(vecs, is.numeric)
      vecs[num_cols] <- lapply(vecs[num_cols], round, 4)
      DT::datatable(vecs, options = list(dom = "t", pageLength = 20),
                    rownames = FALSE)
    })

    output$env_download <- downloadHandler(
      filename = function() paste0("envfit_", Sys.Date(), ".pdf"),
      content  = function(file) {
        req(env_result())
        ggplot2::ggsave(file, env_result()$plot, width = 9, height = 8)
      }
    )

    # ── Enable download buttons ──────────────────────────────────────────────
    observe({
      shinyjs::toggleState("ord_download",  !is.null(ord_plot_obj()))
      shinyjs::toggleState("perm_download", !is.null(perm_result()))
      shinyjs::toggleState("disp_download", !is.null(disp_plot_obj()))
      shinyjs::toggleState("heat_download", !is.null(heat_plot_obj()))
      shinyjs::toggleState("env_download",  !is.null(env_result()))
    })
  })

}
