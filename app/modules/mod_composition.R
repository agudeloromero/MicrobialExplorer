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
        col_widths = c(3, 9, 12),

        card(
          class = "border-0 bg-light",
          card_body(
            tags$h6(class = "fw-bold", "Parameters"),
            uiOutput(ns("bars_rank_ui")),
            numericInput(ns("bars_top_n"), "Top N taxa",
                         value = 20, min = 5, max = 50, step = 5),
            selectInput(ns("bars_sort"), "Sort samples by",
                        choices = c("Group" = "group",
                                    "Dominant taxon" = "dominant_taxon",
                                    "None" = "none")),
            uiOutput(ns("bars_facet_ui")),
            actionButton(ns("bars_reset"),
                         label = tagList(icon("rotate-left"), " Reset to defaults"),
                         class = "btn btn-outline-secondary btn-sm w-100"),
            hr(),
            downloadButton(ns("dl_bars"), "Download plot",
                           class = "btn btn-outline-secondary btn-sm w-100")
          )
        ),

        card(
          card_body(
            plotOutput(ns("bars_plot"), height = "550px")
          )
        ),

        # Help
        card(
          card_header(icon("circle-info"), " How to read this tab"),
          card_body(
            class = "small",
            layout_columns(
              col_widths = c(6, 6),
              tags$div(
                tags$h6(class = "fw-bold", "What you are looking at"),
                tags$ul(
                  tags$li(tags$strong("The question: "),
                          "which taxa make up each sample, and do the groups look different?"),
                  tags$li(tags$strong("Each bar is one sample. "),
                          "Each colour is one taxon, and its height is the share of that ",
                          "sample's reads belonging to it. Every bar adds up to 100%."),
                  tags$li(tags$strong("The coloured bar on top "),
                          "shows the group of each sample (the group variable in the sidebar)."),
                  tags$li(tags$strong("Other (grey): "),
                          "all the taxa outside the top N, added together."),
                  tags$li(tags$strong("Unknown: "),
                          "reads that could not be given a name at this rank. A large Unknown ",
                          "share is common at genus level; try a higher rank such as family."),
                  tags$li(tags$strong("Data in use: "),
                          "the plot uses the dataset from the QC tab, including any sample or ",
                          "taxa filter you applied there.")
                ),
                tags$h6(class = "fw-bold mt-3", "What to look for"),
                tags$ul(
                  tags$li("Samples from the same group that look alike, and groups that ",
                          "look different from each other."),
                  tags$li("A sample that looks nothing like the rest of its group: check it ",
                          "for a mix-up, contamination or low reads."),
                  tags$li("One taxon filling almost a whole bar: real in some sample types, ",
                          "but also a sign of contamination in low-biomass samples.")
                )
              ),
              tags$div(
                tags$h6(class = "fw-bold", "Settings"),
                tags$ul(
                  tags$li(tags$strong("Taxonomic rank: "),
                          "the level the reads are grouped at. Phylum gives the overview; ",
                          "family or genus give detail but more Unknown."),
                  tags$li(tags$strong("Top N taxa: "),
                          "how many of the most abundant taxa get their own colour. 10 to 20 ",
                          "is readable. Above 20 the colours start to repeat."),
                  tags$li(tags$strong("Sort samples by: "),
                          "Group puts samples of the same group side by side. Dominant taxon ",
                          "puts samples with the same most abundant taxon together, which shows ",
                          "patterns that cut across groups. None keeps the original order."),
                  tags$li(tags$strong("Facet by: "),
                          "splits the plot into one panel per level of a variable. Only ",
                          "variables with 12 levels or fewer are offered."),
                  tags$li(tags$strong("Reset to defaults: "),
                          "puts rank, Top N, sorting and facet back to their starting values.")
                ),
                tags$h6(class = "fw-bold mt-3", "Keep in mind"),
                tags$ul(
                  tags$li(tags$strong("These are proportions, not amounts. "),
                          "If one taxon goes up, the others must go down in the plot even if ",
                          "their real numbers did not change."),
                  tags$li(tags$strong("This plot describes; it does not test. "),
                          "Use the DA tab to test whether a taxon differs between groups, and ",
                          "the Beta tab to test whether overall composition differs.")
                )
              )
            )
          )
        )
      )
    ),

    # ── 2. Mean composition ──────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("chart-column"), " Mean composition"),
      value = "mean_comp",

      layout_columns(
        col_widths = c(3, 9, 12),

        card(
          class = "border-0 bg-light",
          card_body(
            tags$h6(class = "fw-bold", "Parameters"),
            uiOutput(ns("mean_rank_ui")),
            numericInput(ns("mean_top_n"), "Top N taxa",
                         value = 10, min = 3, max = 20, step = 1),
            selectInput(ns("mean_error"), "Error bars",
                        choices = c("Standard error (SE)"      = "se",
                                    "Standard deviation (SD)"  = "sd",
                                    "95% confidence interval"  = "ci95",
                                    "None"                     = "none"),
                        selected = "se"),
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
        ),

        # Help
        card(
          card_header(icon("circle-info"), " How to read this tab"),
          card_body(
            class = "small",
            layout_columns(
              col_widths = c(6, 6),
              tags$div(
                tags$h6(class = "fw-bold", "What you are looking at"),
                tags$ul(
                  tags$li(tags$strong("The question: "),
                          "on average, how abundant is each taxon in each group?"),
                  tags$li(tags$strong("Each cluster of bars is one group "),
                          "(the group variable in the sidebar). The number of samples in the ",
                          "group is shown under its name."),
                  tags$li(tags$strong("Each bar is one taxon. "),
                          "Its height is the mean share of reads from that taxon across the ",
                          "samples of the group."),
                  tags$li(tags$strong("Other (grey): "),
                          "all the taxa outside the top N, added together. With Other ",
                          "included, the bars of a group add up to 100%."),
                  tags$li(tags$strong("Unknown: "),
                          "reads that could not be given a name at this rank.")
                ),
                tags$h6(class = "fw-bold mt-3", "Why this is useful"),
                tags$ul(
                  tags$li("It condenses the stacked bars into one profile per group, which is ",
                          "easier to compare and to report."),
                  tags$li("It shows which taxa characterise each group, and the error bars ",
                          "show how consistent that is among the samples of the group."),
                  tags$li("It gives a first idea of which taxa may differ between groups, to ",
                          "follow up in the DA tab.")
                )
              ),
              tags$div(
                tags$h6(class = "fw-bold", "Error bars"),
                tags$ul(
                  tags$li(tags$strong("Standard error (SE): "),
                          "how precisely the group mean is estimated. It gets smaller as the ",
                          "number of samples grows. Use it to compare group means."),
                  tags$li(tags$strong("Standard deviation (SD): "),
                          "how much the samples of a group differ from each other. It does not ",
                          "shrink with more samples. Use it to show variability."),
                  tags$li(tags$strong("95% confidence interval: "),
                          "the range likely to contain the true group mean. It is wide when a ",
                          "group has few samples, which makes small groups visible."),
                  tags$li(tags$strong("None: "),
                          "bars only."),
                  tags$li("Bars stop at 0%, because a proportion cannot be negative. No error ",
                          "bar is drawn for a group with a single sample.")
                ),
                tags$h6(class = "fw-bold mt-3", "Keep in mind"),
                tags$ul(
                  tags$li(tags$strong("A mean hides the individual samples. "),
                          "A tall bar with a long error bar can come from one sample. ",
                          "Check the Stacked bars tab."),
                  tags$li(tags$strong("With 2 or 3 samples per group "),
                          "any error bar is a rough guide only. Always state the number of ",
                          "samples and which error bar you used."),
                  tags$li(tags$strong("Error bars are not a test. "),
                          "Bars that overlap or do not overlap do not tell you whether groups ",
                          "differ. Use the DA tab for that."),
                  tags$li(tags$strong("Data in use: "),
                          "the dataset from the QC tab, including any filter you applied there.")
                )
              )
            )
          )
        )
      )
    ),

    # ── 3. Heatmap ───────────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("grid"), " Heatmap"),
      value = "heatmap",

      layout_columns(
        col_widths = c(3, 9, 12),

        card(
          class = "border-0 bg-light",
          card_body(
            tags$h6(class = "fw-bold", "Parameters"),
            uiOutput(ns("heat_rank_ui")),
            numericInput(ns("heat_top_n"), "Top N taxa",
                         value = 30, min = 5, max = 60, step = 5),
            selectInput(ns("heat_transform"), "Transformation",
                        choices = c("CLR" = "clr",
                                    "Log10" = "log10",
                                    "Relative abundance" = "relative"),
                        selected = "clr"),
            checkboxInput(ns("heat_scale"), "Scale each taxon (row z-score)",
                          value = FALSE),
            uiOutput(ns("heat_group_ui")),
            selectInput(ns("heat_cols"), "Columns (samples)",
                        choices = c("Block by group, cluster within each" = "block",
                                    "Cluster all samples together"        = "cluster",
                                    "Keep the order of the data"          = "none"),
                        selected = "block"),
            checkboxInput(ns("heat_cluster_rows"), "Cluster rows (taxa)",
                          value = TRUE),
            actionButton(ns("heat_reset"),
                         label = tagList(icon("rotate-left"), " Reset to defaults"),
                         class = "btn btn-outline-secondary btn-sm w-100"),
            hr(),
            downloadButton(ns("dl_heat"), "Download plot",
                           class = "btn btn-outline-secondary btn-sm w-100")
          )
        ),

        card(
          card_body(
            plotOutput(ns("heat_plot"), height = "600px")
          )
        ),

        # Help
        card(
          card_header(icon("circle-info"), " How to read this tab"),
          card_body(
            class = "small",
            layout_columns(
              col_widths = c(6, 6),
              tags$div(
                tags$h6(class = "fw-bold", "What you are looking at"),
                tags$ul(
                  tags$li(tags$strong("The question: "),
                          "which taxa are high or low in which samples, and do samples and ",
                          "taxa fall into groups with similar patterns?"),
                  tags$li(tags$strong("Each column is one sample, each row one taxon. "),
                          "The colour of a cell is the abundance of that taxon in that sample: ",
                          "red is high, blue is low, white is in between."),
                  tags$li(tags$strong("The coloured bar on top "),
                          "shows the group of each sample."),
                  tags$li(tags$strong("The tree on the left "),
                          "joins taxa with similar patterns across samples. Taxa on nearby ",
                          "branches rise and fall together."),
                  tags$li(tags$strong("The tree on top "),
                          "(only when all samples are clustered together) joins samples with ",
                          "similar composition. The lower two branches join, the more alike ",
                          "they are.")
                ),
                tags$h6(class = "fw-bold mt-3", "Why this is useful"),
                tags$ul(
                  tags$li("It shows many taxa and all samples at once, where the bar plots ",
                          "show only the most abundant taxa."),
                  tags$li("Blocks of red or blue reveal sets of taxa that mark a group of ",
                          "samples."),
                  tags$li("Clustering all samples together is an unbiased check: if samples ",
                          "from the same group end up side by side without being told to, the ",
                          "groups really differ in composition.")
                )
              ),
              tags$div(
                tags$h6(class = "fw-bold", "Settings"),
                tags$ul(
                  tags$li(tags$strong("Taxonomic rank and Top N taxa: "),
                          "the level and how many of the most abundant taxa to show."),
                  tags$li(tags$strong("Transformation: "),
                          "CLR compares each taxon with the average taxon in the same sample ",
                          "and is the recommended choice. Log10 shows amounts on a compressed ",
                          "scale. Relative abundance shows plain proportions, which a few ",
                          "dominant taxa tend to overwhelm."),
                  tags$li(tags$strong("Scale each taxon: "),
                          "colours then show whether a sample is above or below the average ",
                          "for that taxon (z-score). Use it to see patterns in rare taxa; ",
                          "you can no longer compare one taxon with another."),
                  tags$li(tags$strong("Group: "),
                          "the variable used for the colour bar and for the blocks."),
                  tags$li(tags$strong("Columns, block by group: "),
                          "samples of each group are kept together in their own panel, and ",
                          "ordered by similarity inside it. Use it to compare known groups."),
                  tags$li(tags$strong("Columns, cluster all samples together: "),
                          "samples are ordered only by similarity, with the tree on top. ",
                          "Use it to see whether the groups separate without help."),
                  tags$li(tags$strong("Columns, keep the order of the data: "),
                          "no reordering."),
                  tags$li(tags$strong("Cluster rows: "),
                          "orders taxa by similarity and draws the tree on the left. ",
                          "Untick to keep taxa ordered by abundance."),
                  tags$li(tags$strong("Reset to defaults: "),
                          "puts every setting back to its starting value.")
                ),
                tags$h6(class = "fw-bold mt-3", "Keep in mind"),
                tags$ul(
                  tags$li("Clustering uses Euclidean distance and complete linkage on the ",
                          "values shown. With CLR this is the Aitchison distance."),
                  tags$li("Clustering always produces a tree, even from noise. Treat it as a ",
                          "description; use the Beta tab to test whether groups differ."),
                  tags$li("Sample names are hidden when there are more than 60 samples."),
                  tags$li(tags$strong("Data in use: "),
                          "the dataset from the QC tab, including any filter you applied there.")
                )
              )
            )
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
      # A facet with many levels gives panels too narrow to read
      meta_df <- data.frame(sample_data(upload_data$ps), check.names = FALSE)
      facet_choices <- Filter(function(v) {
        v %in% colnames(meta_df) && dplyr::n_distinct(meta_df[[v]]) <= 12
      }, facet_choices)
      all_vars <- c("None" = "", facet_choices)
      selectInput(ns("bars_facet"), "Facet by",
                  choices = all_vars, selected = "")
    })

    # Reset the stacked-bar parameters to their starting values
    observeEvent(input$bars_reset, {
      feat <- upload_data$features
      preferred <- intersect(c("Genus", "Family", "Order", "Phylum"), feat$avail_ranks)
      updateSelectInput(session, "bars_rank",
                        selected = if (length(preferred) > 0) preferred[1] else feat$avail_ranks[1])
      updateNumericInput(session, "bars_top_n", value = 20)
      updateSelectInput(session, "bars_sort",  selected = "group")
      updateSelectInput(session, "bars_facet", selected = "")
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
    output$bars_rank_ui  <- renderUI({ .rank_selector("bars_rank", "Taxonomic rank") })
    output$mean_rank_ui  <- renderUI({ .rank_selector("mean_rank", "Taxonomic rank") })

    # Heatmap rank and group selectors
    output$heat_rank_ui <- renderUI({ .rank_selector("heat_rank", "Taxonomic rank") })

    output$heat_group_ui <- renderUI({
      feat <- upload_data$features
      req(feat)
      choices <- c("None" = "", feat$group_vars)
      gv      <- parent_input$group_variable
      selectInput(ns("heat_group"), "Group (colour bar and blocks)",
                  choices  = choices,
                  selected = if (!is.null(gv) && gv %in% feat$group_vars) gv else "")
    })

    # Reset the heatmap parameters to their starting values
    observeEvent(input$heat_reset, {
      feat <- upload_data$features
      preferred <- intersect(c("Genus", "Family", "Order", "Phylum"), feat$avail_ranks)
      gv        <- parent_input$group_variable
      updateSelectInput(session, "heat_rank",
                        selected = if (length(preferred) > 0) preferred[1] else feat$avail_ranks[1])
      updateNumericInput(session, "heat_top_n", value = 30)
      updateSelectInput(session, "heat_transform", selected = "clr")
      updateCheckboxInput(session, "heat_scale", value = FALSE)
      updateSelectInput(session, "heat_group",
                        selected = if (!is.null(gv) && gv %in% feat$group_vars) gv else "")
      updateSelectInput(session, "heat_cols", selected = "block")
      updateCheckboxInput(session, "heat_cluster_rows", value = TRUE)
    })

    # Helper: resolve rank
    # each sub-tab has its own rank selector and falls back to the sidebar
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
      rank_use <- resolve_rank(input$bars_rank)
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
      rank_use <- resolve_rank(input$bars_rank)
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
      rank_use <- resolve_rank(input$mean_rank)
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
      rank_use <- resolve_rank(input$mean_rank)
      tryCatch(
        plot_mean_composition(ps_mean(),
                              rank      = rank_use,
                              group_var = group_var(),
                              top_n     = input$mean_top_n,
                              error_bar = input$mean_error %||% "se"),
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
      rank_use <- resolve_rank(input$heat_rank)
      # heat_group gives the colour bar; heat_cols decides how samples are ordered
      hg <- if (is.null(input$heat_group) || input$heat_group == "") NULL
            else input$heat_group
      cols_mode <- input$heat_cols %||% "block"
      tryCatch(
        plot_abundance_heatmap(ps_work(),
                               rank            = rank_use,
                               top_n           = input$heat_top_n,
                               group_var       = hg,
                               transform       = input$heat_transform,
                               cluster_samples = cols_mode %in% c("block", "cluster"),
                               cluster_rows    = isTRUE(input$heat_cluster_rows),
                               block_by_group  = cols_mode == "block" && !is.null(hg),
                               scale_rows      = isTRUE(input$heat_scale)),
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
