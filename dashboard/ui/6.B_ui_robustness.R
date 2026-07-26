page_robustness <- function() {
  layout_sidebar(
    sidebar = sidebar(
      title = "Sweep Parameters",
      selectInput(
        "robust_nutrient", "Select Nutrient (1D Sweep):",
        choices  = c("Oxygen" = "O2", "Glucose" = "Glc", "Maltose" = "Malt", "Ammonia (Nitrogen)" = "NH4"),
        selected = "O2"
      ),
      selectInput(
        "robust_pair", "Select Nutrient Pair (2D Heatmap):",
        choices  = c("Glucose \u00D7 Oxygen" = "glc_o2", "Glucose \u00D7 Ammonia" = "glc_nh4", "Maltose \u00D7 Ammonia" = "malt_nh4"),
        selected = "glc_o2"
      )
    ),
    navset_card_tab(
      nav_panel("Single-Nutrient Sweeps", card(full_screen = TRUE, card_header("1D Nutrient Sweep"), plotlyOutput("robust_line_plot", height = "450px"))),
      nav_panel("Paired-Nutrient Heatmaps", card(full_screen = TRUE, card_header("2D Nutrient Grid Response"), plotlyOutput("robust_heatmap_plot", height = "450px")))
    )
  )
}
