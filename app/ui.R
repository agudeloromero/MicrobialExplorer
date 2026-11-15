# =============================================================================
# ui.R — MicrobialExplorer Shiny App
# =============================================================================

ui <- page_fluid(
  title = APP_TITLE,
  theme = bs_theme(
    version    = 5,
    bootswatch = "flatly",
    primary    = "#2c3e50",
    success    = "#1abc9c",
    danger     = "#e74c3c",
    base_font  = font_google("Inter"),
    code_font  = font_google("Fira Code")
  ),
  class = "me-page",

  # useShinyjs must be called inside the UI
  useShinyjs(),
  # styles are placed in the page itself, so a stylesheet kept by the
  # browser from an earlier session can never be used
  tags$head(
    includeCSS("www/custom.css")
  ),

  # ── Top banner: title on the left, housekeeping links on the right ───────
  tags$div(
    class = "me-banner",
    tags$span(
      class = "me-brand",
      if (file.exists("www/logo.png"))
        tags$img(src = "logo.png", height = "28px", style = "margin-right:8px;"),
      APP_TITLE,
      tags$small(class = "me-version ms-2 fs-6 fw-normal", APP_VERSION)
    ),
    tags$div(
      class = "me-banner-links",
      actionLink("nav_home",   label = tagList(icon("house"), " Home")),
      actionLink("nav_export", label = tagList(icon("download"), " Export")),
      tags$a(
        href   = "https://github.com/agudeloromero/MicrobialExplorer",
        target = "_blank",
        icon("github"), " GitHub"
      )
    )
  ),

  # ── Strip: group variable and data in use (once data is loaded) ──────────
  conditionalPanel(
    condition = "output.data_loaded",
    tags$div(
      class = "me-strip",
      tags$div(
        class = "me-strip-group",
        tooltip(
          tags$span(class = "me-strip-label", "Group variable ",
                    icon("circle-question")),
          "Categorical metadata column used for group comparisons across all modules."
        ),
        selectInput("group_variable", label = NULL,
                    choices = NULL, selected = NULL, width = "240px")
      ),
      uiOutput("data_in_use", class = "me-strip-data")
    )
  ),

  # ── Analysis menu (left) and pages ───────────────────────────────────────
  tags$div(
    class = "me-main",
    navset_pill_list(
      id       = "top_nav",
      selected = "home",
      widths   = c(2, 10),
      well     = FALSE,

  # 1. Welcome / upload landing (opened from the banner)
  nav_panel_hidden(
    value = "home",
    mod_upload_landing_ui("upload")
  ),

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

  # Export (opened from the banner)
  nav_panel_hidden(
    value = "export",
    uiOutput("export_ui")
  )

    )
  )
)
