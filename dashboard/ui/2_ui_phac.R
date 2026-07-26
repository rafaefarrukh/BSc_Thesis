###################################################################################################
# UI: phaC
###################################################################################################

page_phac <- function() {
  layout_sidebar(
    sidebar = sidebar(
      title = "Structure Controls",
      width = 300,

      tags$h6("Representation"),
      selectInput(
        "phac_style_select", "Style",
        choices  = c("Cartoon", "Stick", "Sphere", "Line"),
        selected = "Cartoon"
      ),
      selectInput(
        "phac_color_select", "Color Scheme",
        choices  = c("Uniform (grey)", "Spectrum (rainbow)"),
        selected = "Uniform (grey)"
      ),
      selectInput(
        "phac_bg_select", "Background",
        choices  = c("White", "Black", "Grey"),
        selected = "White"
      ),
      checkboxInput("phac_spin_toggle", "Spin structure", value = FALSE),

      hr(),
      tags$h6("View"),
      layout_columns(
        col_widths = c(6, 6),
        actionButton("phac_btn_zoom_in",  "Zoom In",  icon = icon("magnifying-glass-plus"),  width = "100%"),
        actionButton("phac_btn_zoom_out", "Zoom Out", icon = icon("magnifying-glass-minus"), width = "100%")
      ),
      actionButton(
        "phac_btn_reset", "Reset View",
        icon = icon("rotate-left"), class = "btn-outline-secondary", width = "100%"
      ),

      hr(),
      tags$h6("Catalytic Triad"),
      actionButton("phac_btn_full",  "Show Full Structure",       class = "btn-secondary",        width = "100%"),
      actionButton("phac_btn_triad", "Highlight Catalytic Triad", class = "btn-primary",          width = "100%"),
      actionButton("phac_btn_zoom",  "Zoom to Triad",             class = "btn-outline-primary",  width = "100%"),
      tags$small(
        tags$strong("Catalytic triad"), tags$br(),
        "Cys133 (nucleophile)", tags$br(),
        "His318 (general base)", tags$br(),
        "Asp289 (orienting residue)"
      )
    ),

    card(
      card_header("phaC — Predicted Structure (ColabFold / AlphaFold2, rank_001)"),
      r3dmolOutput("phac_structure_viewer", height = "650px")
    )
  )
}

###################################################################################################
