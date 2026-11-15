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
          col_widths = c(3, 9, 12),

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
                    tags$li(tags$strong("Top plot: "),
                            "one bar per sample, sorted from fewest to most reads. ",
                            "The red dashed line is your threshold; samples under it are flagged."),
                    tags$li(tags$strong("Bottom plot: "),
                            "how many samples fall in each read-depth range. ",
                            "Red line = threshold, green line = median."),
                    tags$li(tags$strong("Summary: "),
                            "Samples, Min, Median and Max reads describe the data and do not ",
                            "change with the threshold. Only \"Below threshold\" and ",
                            "\"Flagged samples\" do."),
                    tags$li(tags$strong("Nothing is removed here. "),
                            "This tab only flags samples. To remove them, use the ",
                            "Sample filter tab.")
                  )
                ),
                tags$div(
                  tags$h6(class = "fw-bold", "Choosing a minimum"),
                  tags$ul(
                    tags$li("There is no universal value: it depends on the sequencing run ",
                            "and the sample type."),
                    tags$li("For 16S amplicon data, 1,000 reads is a commonly used lower limit, ",
                            "and 5,000 to 10,000 is often used for well-sequenced samples such as stool."),
                    tags$li("Low-biomass samples (for example milk, airway or skin) often have ",
                            "fewer reads; a lower threshold may be needed, and negative controls ",
                            "become more important."),
                    tags$li("Use the top plot: set the threshold in the gap between the few ",
                            "samples that sit far below the rest and the main group."),
                    tags$li("Check the Rarefaction tab: the threshold should be at a depth ",
                            "where the curves have flattened."),
                    tags$li("A higher threshold gives more reliable samples but fewer of them.")
                  )
                )
              )
            )
          )
        )
      ),

      # ── 2. Rarefaction curves ────────────────────────────────────────────
      nav_panel(
        title = tagList(icon("chart-line"), " Rarefaction"),
        value = "rarefaction",

        layout_columns(
          col_widths = c(3, 9, 12),

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
                            "were the samples sequenced deeply enough to capture their diversity?"),
                    tags$li(tags$strong("Each line is one sample. "),
                            "It shows how many ASVs would be found if only that many reads ",
                            "had been sequenced."),
                    tags$li(tags$strong("A flat line is good: "),
                            "more reads would add few new ASVs. A line still climbing at the ",
                            "right edge means deeper sequencing would find more."),
                    tags$li(tags$strong("The curves stop at the shallowest sample. "),
                            "One very shallow sample cuts all curves short."),
                    tags$li(tags$strong("Orange dotted line (knee point): "),
                            "for each sample, the depth where the curve gains about 1 new ASV ",
                            "per 100 reads; the line is the median across samples. It is a guide ",
                            "only: if the curves are still rising there, they have not levelled off."),
                    tags$li(tags$strong("Nothing is changed here. "),
                            "This tab only draws the curves; the data are not rarefied or filtered.")
                  )
                ),
                tags$div(
                  tags$h6(class = "fw-bold", "Settings and what to do next"),
                  tags$ul(
                    tags$li(tags$strong("Step size: "),
                            "the gap in reads between calculated points. 500 suits most data. ",
                            "Use 100 if your shallowest sample has only a few thousand reads; ",
                            "use a larger step to make deep datasets run faster."),
                    tags$li(tags$strong("Max samples to plot: "),
                            "leave blank to draw every sample. With hundreds of samples, ",
                            "enter a number such as 50 to draw a random subset and save time."),
                    tags$li(tags$strong("If most curves are flat: "),
                            "the sequencing depth is sufficient. Note the depth where they flatten."),
                    tags$li(tags$strong("If a few samples are much shallower than the rest: "),
                            "use that depth as the minimum reads in the Sample filter tab, ",
                            "so those samples are removed."),
                    tags$li(tags$strong("If most curves are still climbing: "),
                            "richness is underestimated. Compare richness between groups with ",
                            "caution, or use evenness-based measures such as Shannon, which are ",
                            "less sensitive to depth."),
                    tags$li("Denoised ASV data without singletons tends to flatten sooner than ",
                            "OTU data, so flat curves are less informative there.")
                  )
                )
              )
            )
          )
        )
      ),

      # ── 3. Sample filtering ──────────────────────────────────────────────
      nav_panel(
        title = tagList(icon("users-slash"), " Sample filter"),
        value = "sample_filter",

        layout_columns(
          col_widths = c(3, 9, 12),

          card(
            class = "border-0 bg-light",
            card_body(
              tags$h6(class = "fw-bold", "Filter thresholds"),
              uiOutput(ns("samp_min_reads_ui")),
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
              actionButton(
                ns("remove_sample_filter"),
                label = tagList(icon("rotate-left"), " Remove filter"),
                class = "btn btn-outline-secondary btn-sm w-100 mt-2"
              ),
              uiOutput(ns("sample_filter_status")),
              hr(),
              tags$h6(class = "fw-bold", "Preview with these values"),
              uiOutput(ns("sample_filter_summary"))
            )
          ),

          card(
            card_body(
              tags$h6(class = "fw-bold mb-3", "Samples that would be removed"),
              DT::dataTableOutput(ns("removed_samples_table")),
              hr(),
              tags$h6(class = "fw-bold mb-3", "Samples that would be kept"),
              DT::dataTableOutput(ns("retained_samples_table"))
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
                            "which samples have too little data to be analysed reliably?"),
                    tags$li(tags$strong("Min reads: "),
                            "samples with fewer total reads are removed. Choose the value from ",
                            "the Depth and Rarefaction tabs. It is not used when the data are ",
                            "not raw counts."),
                    tags$li(tags$strong("Min observed taxa: "),
                            "samples with fewer different taxa are removed. This catches failed ",
                            "or nearly empty samples."),
                    tags$li(tags$strong("Preview: "),
                            "the counts and the two tables show what the current values would do. ",
                            "Check which samples would be lost and from which group.")
                  ),
                  tags$h6(class = "fw-bold mt-3", "You are in control"),
                  tags$ul(
                    tags$li(tags$strong("Nothing is removed until you click Apply filter. "),
                            "Changing the values only changes the preview."),
                    tags$li(tags$strong("Apply filter: "),
                            "the analysis will use only the samples that pass. ",
                            "The status under the buttons shows what is applied."),
                    tags$li(tags$strong("Remove filter: "),
                            "goes back to using all samples."),
                    tags$li("Loading a new dataset starts again with no filter.")
                  )
                ),
                tags$div(
                  tags$h6(class = "fw-bold", "How this affects the rest of the analysis"),
                  tags$ul(
                    tags$li(tags$strong("Rarefaction depth: "),
                            "alpha and beta diversity rarefy every sample to the depth of the ",
                            "shallowest sample that remains. Removing a few shallow samples lets ",
                            "all the others keep more reads."),
                    tags$li(tags$strong("Statistical power: "),
                            "every removed sample is one fewer in its group. Check that no group ",
                            "becomes too small to compare."),
                    tags$li(tags$strong("Balance: "),
                            "if most removed samples come from one group, ask why before removing ",
                            "them. Low reads can be a real feature of that group (for example ",
                            "low-biomass samples)."),
                    tags$li(tags$strong("Taxa filter: "),
                            "the sample filter is used first, so prevalence in the Taxa filter tab ",
                            "is calculated on the samples that remain."),
                    tags$li(tags$strong("A stricter filter "),
                            "gives more reliable samples but fewer of them. ",
                            "Report the thresholds you used and how many samples were removed.")
                  )
                )
              )
            )
          )
        )
      ),

      # ── 4. Taxa filtering ────────────────────────────────────────────────
      nav_panel(
        title = tagList(icon("filter"), " Taxa filter"),
        value = "taxa_filter",

        layout_columns(
          col_widths = c(3, 9, 12),

          card(
            class = "border-0 bg-light",
            card_body(
              tags$h6(class = "fw-bold", "Filter thresholds"),
              sliderInput(
                ns("taxa_prev"),
                label = "Min prevalence",
                min = 0.01, max = 0.5,
                value = 0.10, step = 0.01,
                ticks = FALSE
              ),
              uiOutput(ns("taxa_abund_ui")),
              actionButton(
                ns("apply_taxa_filter"),
                label = tagList(icon("check"), " Apply filter"),
                class = "btn btn-success btn-sm w-100 mt-2"
              ),
              actionButton(
                ns("remove_taxa_filter"),
                label = tagList(icon("rotate-left"), " Remove filter"),
                class = "btn btn-outline-secondary btn-sm w-100 mt-2"
              ),
              uiOutput(ns("taxa_filter_status")),
              hr(),
              tags$h6(class = "fw-bold", "Preview with these values"),
              uiOutput(ns("taxa_filter_summary"))
            )
          ),

          card(
            card_body(
              plotOutput(ns("taxa_filter_plot"), height = "500px")
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
                            "which taxa are too rare to be analysed reliably?"),
                    tags$li(tags$strong("Each point is one taxon. "),
                            "Left to right: the fraction of samples it is found in (prevalence). ",
                            "Bottom to top: its total reads across all samples, on a log scale ",
                            "(2 = 100 reads, 4 = 10,000 reads)."),
                    tags$li(tags$strong("Red dashed lines: "),
                            "vertical = minimum prevalence, horizontal = minimum total abundance. ",
                            "A taxon is kept only if it passes both, so the kept taxa are in the ",
                            "top-right area. Crosses are removed taxa."),
                    tags$li(tags$strong("Preview: "),
                            "the plot and the counts show what the current slider values would do. ",
                            "They update as you move the sliders."),
                    tags$li(tags$strong("Taxa are removed, not samples. "),
                            "Rare taxa are usually many in number but hold few reads, so a large ",
                            "share of taxa can be removed while most reads are kept.")
                  ),
                  tags$h6(class = "fw-bold mt-3", "You are in control"),
                  tags$ul(
                    tags$li(tags$strong("Nothing is removed until you click Apply filter. "),
                            "Moving the sliders only changes the preview."),
                    tags$li(tags$strong("Apply filter: "),
                            "the analysis will use the dataset filtered with these values. ",
                            "The status under the buttons shows what is applied."),
                    tags$li(tags$strong("Remove filter: "),
                            "goes back to using all taxa."),
                    tags$li(tags$strong("Try different values: "),
                            "change the sliders, click Apply filter again and look at how your ",
                            "results change. Conclusions that hold across reasonable values are ",
                            "the ones to trust."),
                    tags$li("If a sample filter is applied, it is used first: prevalence is ",
                            "calculated on the samples that remain."),
                    tags$li("Loading a new dataset starts again with no filter.")
                  ),
                  tags$h6(class = "fw-bold mt-3", "Choosing values"),
                  tags$ul(
                    tags$li(tags$strong("Min prevalence: "),
                            "0.05 to 0.10 (present in 5 to 10% of samples) is commonly used. ",
                            "Think in samples: 0.10 of 26 samples means at least 3 samples."),
                    tags$li("If your groups are very different (for example body sites or ",
                            "environments), a taxon can be real but present in one group only. ",
                            "Use a lower prevalence so group-specific taxa are not lost."),
                    tags$li(tags$strong("Min total abundance: "),
                            "10 reads removes taxa seen only a handful of times. ",
                            "Raise it for very deeply sequenced data. It is not used when the ",
                            "data are not raw counts.")
                  )
                ),
                tags$div(
                  tags$h6(class = "fw-bold", "How this affects the rest of the analysis"),
                  tags$ul(
                    tags$li(tags$strong("Differential abundance: "),
                            "benefits most. Fewer taxa are tested, so the correction for ",
                            "multiple testing is less severe, and results driven by taxa that ",
                            "are almost all zeros are avoided."),
                    tags$li(tags$strong("Networks and correlations: "),
                            "rare taxa produce unreliable correlations; filtering makes these ",
                            "analyses more stable and faster."),
                    tags$li(tags$strong("Beta diversity: "),
                            "abundance-based distances (Bray-Curtis, weighted UniFrac) change little. ",
                            "Presence/absence distances (Jaccard, unweighted UniFrac) are more ",
                            "sensitive to rare taxa and change more."),
                    tags$li(tags$strong("Alpha diversity: "),
                            "richness (Observed, Chao1) drops when rare taxa are removed, so alpha ",
                            "diversity is normally calculated on unfiltered data. Shannon and ",
                            "Simpson change little."),
                    tags$li(tags$strong("Composition: "),
                            "barplots at phylum or genus level barely change, because the removed ",
                            "taxa hold few reads."),
                    tags$li(tags$strong("A stricter filter "),
                            "gives cleaner statistics but can remove real low-abundance organisms. ",
                            "Report the thresholds you used.")
                  ),
                  tags$h6(class = "fw-bold mt-3", "Where the filtered data are used"),
                  tags$ul(
                    tags$li("The filtered dataset is used by the Composition tab and by the ",
                            "download on the Final dataset tab."),
                    tags$li("The Alpha, Beta and DA tabs currently use the data as uploaded.")
                  )
                )
              )
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
        ),

        # Help
        card(
          card_header(icon("circle-info"), " How to read this tab"),
          card_body(
            class = "small",
            layout_columns(
              col_widths = c(6, 6),
              tags$div(
                tags$h6(class = "fw-bold", "QC summary table"),
                tags$ul(
                  tags$li(tags$strong("Uploaded data: "),
                          "the dataset as you loaded it."),
                  tags$li(tags$strong("Sample filter / Taxa filter: "),
                          "the dataset after each filter. A row marked \"not applied\" ",
                          "means that filter is not in use, so its numbers repeat the row above."),
                  tags$li(tags$strong("Final dataset (in use): "),
                          "what the analysis works with after the filters you applied."),
                  tags$li(tags$strong("What to check: "),
                          "how many samples and reads you keep. Losing many taxa but few reads ",
                          "is normal. Losing a large share of the reads or of the samples means ",
                          "the thresholds are probably too strict for this dataset.")
                ),
                tags$h6(class = "fw-bold mt-3", "Download"),
                tags$ul(
                  tags$li("The file contains the final dataset exactly as shown in the last row, ",
                          "as an R phyloseq object. Open it in R with readRDS().")
                )
              ),
              tags$div(
                tags$h6(class = "fw-bold", "Reads vs taxon richness"),
                tags$ul(
                  tags$li(tags$strong("The question: "),
                          "does the number of taxa found in a sample depend on how deeply it ",
                          "was sequenced?"),
                  tags$li(tags$strong("Each point is one sample "),
                          "in the final dataset. Left to right: its total reads. Bottom to top: ",
                          "the number of different taxa found in it. The dashed line is the ",
                          "overall trend and the grey band its uncertainty."),
                  tags$li(tags$strong("Good: "),
                          "a flat or gently sloping line with points scattered around it. ",
                          "Richness then reflects the samples, not the sequencing depth."),
                  tags$li(tags$strong("A warning sign: "),
                          "a clear upward slope with points close to the line. Samples with more ",
                          "reads show more taxa simply because they were sequenced more, so ",
                          "richness must be compared at equal depth (rarefaction)."),
                  tags$li(tags$strong("Also look for: "),
                          "points alone in the bottom-left corner (few reads and few taxa). ",
                          "These are candidates for the Sample filter."),
                  tags$li(tags$strong("Keep in mind: "),
                          "very different sample types have very different richness, which ",
                          "spreads the points vertically. That is biology, not a quality problem."),
                  tags$li(tags$strong("Use it to compare settings: "),
                          "apply a filter and come back. Removing rare taxa usually flattens ",
                          "the trend, because rare taxa are the ones most dependent on depth.")
                )
              )
            )
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

    # ── FILTER STATE ──────────────────────────────────────────────────────
    # Nothing is filtered until the user clicks Apply.
    # NULL = no filter applied; otherwise the thresholds in use.
    samp_applied <- reactiveVal(NULL)
    taxa_applied <- reactiveVal(NULL)

    # A new dataset starts with no filter
    observeEvent(upload_data$ps, {
      samp_applied(NULL)
      taxa_applied(NULL)
    }, ignoreNULL = FALSE)

    # Thresholds in reads only make sense for raw counts
    is_counts <- reactive({
      identical(upload_data$features$data_type, "counts")
    })

    not_counts_note <- function(label) {
      tags$p(class = "text-muted small",
             tags$strong(label), tags$br(),
             "Not used: the data are not raw counts.")
    }

    filter_status <- function(applied_txt, changed, n_txt) {
      if (is.null(applied_txt)) {
        tags$div(class = "alert alert-secondary p-2 small mt-2 mb-0",
                 icon("circle-info"), " ",
                 tags$strong("Not applied. "), n_txt)
      } else {
        tagList(
          tags$div(class = "alert alert-success p-2 small mt-2 mb-0",
                   icon("circle-check"), " ",
                   tags$strong("Applied: "), paste0(applied_txt, ". ", n_txt)),
          if (changed)
            tags$div(class = "alert alert-warning p-2 small mt-2 mb-0",
                     icon("triangle-exclamation"),
                     " Values changed. Click Apply filter to use them.")
        )
      }
    }

    # ── 3. SAMPLE FILTERING — live preview, used only after Apply ─────────
    output$samp_min_reads_ui <- renderUI({
      req(ps_raw())
      if (is_counts()) {
        numericInput(ns("samp_min_reads"), label = "Min reads",
                     value = 1000, min = 0, step = 100)
      } else {
        not_counts_note("Min reads")
      }
    })

    samp_values <- reactive({
      reads <- if (is_counts()) input$samp_min_reads else 0
      taxa  <- input$samp_min_taxa
      req(!is.null(reads), !is.na(reads), !is.null(taxa), !is.na(taxa))
      list(reads = reads, taxa = taxa)
    })
    samp_values_d <- debounce(samp_values, 500)

    # Metadata rows for a set of samples, with the sample name as first column.
    # The metadata may already have a column called "Sample".
    sample_meta_table <- function(ps, ids) {
      df <- data.frame(sample_data(ps), check.names = FALSE)[ids, , drop = FALSE]
      cbind(data.frame(`Sample name` = rownames(df), check.names = FALSE),
            df, row.names = NULL)
    }

    run_sample_filter <- function(ps, th) {
      tryCatch(
        qc_filter_samples(ps, min_reads = th$reads, min_taxa = th$taxa),
        error = function(e) { message("[mod_qc] sample filter: ", e$message); NULL }
      )
    }

    sample_result <- reactive({
      req(ps_raw())
      run_sample_filter(ps_raw(), samp_values_d())
    })

    observeEvent(input$apply_sample_filter, {
      th  <- samp_values()
      res <- run_sample_filter(ps_raw(), th)
      if (is.null(res) || is.null(res$ps_filtered)) {
        showNotification("These values would remove every sample. Filter not applied.",
                         type = "error")
        return()
      }
      samp_applied(th)
    })

    observeEvent(input$remove_sample_filter, {
      samp_applied(NULL)
    })

    # Dataset after the applied sample filter (all samples if none applied)
    ps_after_samples <- reactive({
      req(ps_raw())
      th <- samp_applied()
      if (is.null(th)) return(ps_raw())
      res <- run_sample_filter(ps_raw(), th)
      if (is.null(res) || is.null(res$ps_filtered)) ps_raw() else res$ps_filtered
    })

    output$sample_filter_status <- renderUI({
      req(ps_raw())
      th <- samp_applied()
      n_txt <- paste0(nsamples(ps_after_samples()), " of ", nsamples(ps_raw()),
                      " samples are in use.")
      if (is.null(th)) return(filter_status(NULL, FALSE, n_txt))
      cur <- samp_values()
      applied_txt <- paste0(
        if (is_counts()) paste0("min reads ", format(th$reads, big.mark = ","), ", ") else "",
        "min taxa ", th$taxa)
      filter_status(applied_txt, !identical(cur, th), n_txt)
    })

    output$sample_filter_summary <- renderUI({
      req(sample_result())
      res <- sample_result()
      if (is.null(res$ps_filtered))
        return(tags$span(class = "text-danger small",
                         icon("triangle-exclamation"),
                         " These values would remove every sample."))
      n_removed <- length(res$removed_samples)
      tagList(
        stat_badge("Samples before", nsamples(ps_raw())),
        stat_badge("Samples after",  nsamples(res$ps_filtered),
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
      sample_meta_table(ps_raw(), removed) %>%
        DT::datatable(options = list(pageLength = 10, dom = "tp"), rownames = FALSE)
    })

    output$retained_samples_table <- DT::renderDataTable({
      req(sample_result())
      if (is.null(sample_result()$ps_filtered))
        return(DT::datatable(data.frame(Message = "No samples kept"),
                             rownames = FALSE, options = list(dom = "t")))
      sample_meta_table(ps_raw(), sample_names(sample_result()$ps_filtered)) %>%
        DT::datatable(options = list(pageLength = 10, dom = "tp"), rownames = FALSE)
    })

    # ── 4. TAXA FILTERING — live preview, used only after Apply ───────────
    # Runs on the samples kept by the applied sample filter.
    output$taxa_abund_ui <- renderUI({
      req(ps_raw())
      if (is_counts()) {
        sliderInput(ns("taxa_abund"), label = "Min total abundance",
                    min = 1, max = 100, value = 10, step = 1, ticks = FALSE)
      } else {
        not_counts_note("Min total abundance")
      }
    })

    taxa_values <- reactive({
      prev  <- input$taxa_prev
      abund <- if (is_counts()) input$taxa_abund else 0
      req(!is.null(prev), !is.null(abund))
      list(prev = prev, abund = abund)
    })
    taxa_values_d <- debounce(taxa_values, 500)

    run_taxa_filter <- function(ps, th, plot_filter = TRUE) {
      tryCatch(
        qc_filter_taxa(ps,
                       min_prevalence = th$prev,
                       min_abundance  = th$abund,
                       plot_filter    = plot_filter),
        error = function(e) { message("[mod_qc] taxa filter: ", e$message); NULL }
      )
    }

    taxa_result <- reactive({
      req(ps_after_samples())
      run_taxa_filter(ps_after_samples(), taxa_values_d())
    })

    observeEvent(input$apply_taxa_filter, {
      th  <- taxa_values()
      res <- run_taxa_filter(ps_after_samples(), th, plot_filter = FALSE)
      if (is.null(res)) {
        showNotification("These values would remove every taxon. Filter not applied.",
                         type = "error")
        return()
      }
      taxa_applied(th)
    })

    observeEvent(input$remove_taxa_filter, {
      taxa_applied(NULL)
    })

    output$taxa_filter_plot <- renderPlot({
      req(taxa_result())
      taxa_result()$plot
    })

    output$taxa_filter_status <- renderUI({
      req(ps_final())
      th <- taxa_applied()
      n_txt <- paste0(format(ntaxa(ps_final()), big.mark = ","), " of ",
                      format(ntaxa(ps_raw()), big.mark = ","), " taxa are in use.")
      if (is.null(th)) return(filter_status(NULL, FALSE, n_txt))
      cur <- taxa_values()
      applied_txt <- paste0(
        "prevalence ", th$prev,
        if (is_counts()) paste0(", abundance ", th$abund) else "")
      filter_status(applied_txt, !identical(cur, th), n_txt)
    })

    output$taxa_filter_summary <- renderUI({
      if (is.null(taxa_result()))
        return(tags$span(class = "text-danger small",
                         icon("triangle-exclamation"),
                         " These values would remove every taxon."))
      s <- taxa_result()$summary
      tagList(
        stat_badge("Taxa before", s$n_before),
        stat_badge("Taxa after",  s$n_after, colour = "success"),
        stat_badge("Removed",
                   paste0(s$n_removed, " (", s$pct_removed, "%)"),
                   colour = if (s$n_removed > 0) "warning" else "success")
      )
    })

    # ── 5. FINAL DATASET — applied sample filter, then applied taxa filter ─
    ps_final <- reactive({
      req(ps_after_samples())
      ps_f <- ps_after_samples()

      th <- taxa_applied()
      if (!is.null(th)) {
        res <- run_taxa_filter(ps_f, th, plot_filter = FALSE)
        if (!is.null(res)) ps_f <- res$ps_filtered
      }

      ps_f
    })

    output$qc_summary_table <- DT::renderDataTable({
      req(ps_final())
      ps_s <- ps_after_samples()
      ps_f <- ps_final()
      stage_row <- function(ps) {
        c(Samples      = nsamples(ps),
          Taxa         = format(ntaxa(ps), big.mark = ","),
          Total_reads  = format(round(sum(sample_sums(ps))), big.mark = ","),
          Median_reads = format(round(median(sample_sums(ps))), big.mark = ","))
      }
      df <- data.frame(
        Step = c("Uploaded data",
                 if (is.null(samp_applied())) "Sample filter: not applied"
                 else "After sample filter",
                 if (is.null(taxa_applied())) "Taxa filter: not applied"
                 else "After taxa filter",
                 "Final dataset (in use)"),
        rbind(stage_row(ps_raw()), stage_row(ps_s), stage_row(ps_f), stage_row(ps_f)),
        stringsAsFactors = FALSE
      )
      DT::datatable(df, options = list(dom = "t", ordering = FALSE), rownames = FALSE) %>%
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
