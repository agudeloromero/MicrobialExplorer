# =============================================================================
# ui.R — MicrobialExplorer Shiny App
# =============================================================================

ui <- page_navbar(
  title = tags$span(
    tags$img(src = "logo.png", height = "28px", style = "margin-right:8px;"),
    APP_TITLE,
    tags$small(class = "text-muted ms-2 fs-6 fw-normal", APP_VERSION)
  ),
  window_title = APP_TITLE,
  theme = bs_theme(
    version    = 5,
    bootswatch = "flatly",
    primary    = "#2c3e50",
    success    = "#1abc9c",
    danger     = "#e74c3c",
    base_font  = font_google("Inter"),
    code_font  = font_google("Fira Code")
  ),
  fillable = FALSE,

  # ── Sidebar ──────────────────────────────────────────────────────────────
  sidebar = sidebar(
    width = 300,
    open  = "open",

    # Upload module UI (defined in mod_upload.R)
    mod_upload_ui("upload"),

    hr(),

    # ── Data summary card — populated by server once data is loaded ──────
    conditionalPanel(
      condition = "output.data_loaded",

      card(
        class = "border-0 bg-light p-2",
        card_body(
          class = "p-2",
          tags$h6(class = "fw-bold mb-2",
                  icon("database"), " Dataset summary"),
          uiOutput("sidebar_summary")
        )
      ),

      hr(),

      # ── Global selectors ─────────────────────────────────────────────
      tags$h6(class = "fw-bold mb-2 text-muted text-uppercase small",
              "Analysis parameters"),

      selectInput(
        "group_variable",
        label = tooltip(
          tags$span("Group variable ", icon("circle-question")),
          "Categorical metadata column used for group comparisons across all modules."
        ),
        choices  = NULL,
        selected = NULL
      ),

      selectInput(
        "tax_rank",
        label = tooltip(
          tags$span("Taxonomic rank ", icon("circle-question")),
          "Rank used for composition and differential abundance analyses."
        ),
        choices  = NULL,
        selected = NULL
      ),

      sliderInput(
        "prev_threshold",
        label = tooltip(
          tags$span("Prevalence threshold ", icon("circle-question")),
          "Minimum proportion of samples a taxon must appear in (used for filtering and networks)."
        ),
        min   = 0.01, max = 0.5,
        value = 0.10, step = 0.01,
        ticks = FALSE
      ),

      sliderInput(
        "abund_threshold",
        label = tooltip(
          tags$span("Abundance threshold ", icon("circle-question")),
          "Minimum read count a taxon must have in at least one sample."
        ),
        min   = 1, max = 50,
        value = 5, step = 1,
        ticks = FALSE
      ),

      # ── Warnings panel ────────────────────────────────────────────────
      uiOutput("sidebar_warnings")
    )
  ),

  # ── Nav panels ───────────────────────────────────────────────────────────

  # 1. Welcome / upload landing
  nav_panel(
    title = tagList(icon("house"), " Home"),
    value = "home",
    mod_upload_landing_ui("upload")
  ),

  nav_spacer(),

  # 2. Dashboard
  nav_panel(
    title = tagList(icon("chart-pie"), " Dashboard"),
    value = "dashboard",
    uiOutput("dashboard_ui")
  ),

  # 3. QC
  nav_panel(
    title = tagList(icon("microscope"), " QC"),
    value = "qc",
    uiOutput("qc_ui")
  ),

  # 4. Composition
  nav_panel(
    title = tagList(icon("layer-group"), " Composition"),
    value = "composition",
    uiOutput("composition_ui")
  ),

  # 5. Alpha Diversity
  nav_panel(
    title = tagList(icon("chart-bar"), " Alpha"),
    value = "alpha",
    uiOutput("alpha_ui")
  ),

  # 6. Beta Diversity
  nav_panel(
    title = tagList(icon("circle-dot"), " Beta"),
    value = "beta",
    uiOutput("beta_ui")
  ),

  # 7. Differential Abundance
  nav_panel(
    title = tagList(icon("not-equal"), " DA"),
    value = "da",
    uiOutput("da_ui")
  ),

  # 8. Functional
  nav_panel(
    title = tagList(icon("dna"), " Functional"),
    value = "functional",
    uiOutput("functional_ui")
  ),

  # 9. Network
  nav_panel(
    title = tagList(icon("circle-nodes"), " Network"),
    value = "network",
    uiOutput("network_ui")
  ),

  # 10. ML
  nav_panel(
    title = tagList(icon("robot"), " ML"),
    value = "ml",
    uiOutput("ml_ui")
  ),

  # 11. Correlation
  nav_panel(
    title = tagList(icon("arrows-left-right"), " Correlation"),
    value = "correlation",
    uiOutput("correlation_ui")
  ),

  # 12. Longitudinal
  nav_panel(
    title = tagList(icon("timeline"), " Longitudinal"),
    value = "longitudinal",
    uiOutput("longitudinal_ui")
  ),

  nav_spacer(),

  # Export
  nav_panel(
    title = tagList(icon("download"), " Export"),
    value = "export",
    uiOutput("export_ui")
  ),

  # ── Footer nav items ─────────────────────────────────────────────────────
  nav_item(
    tags$a(
      href   = "https://github.com/agudeloromero/MicrobialExplorer",
      target = "_blank",
      icon("github"), " GitHub"
    )
  ),

  # useShinyjs must be called inside the UI
  header = tagList(
    useShinyjs(),
    tags$head(
      tags$link(rel = "stylesheet", href = "custom.css")
    )
  )
)
