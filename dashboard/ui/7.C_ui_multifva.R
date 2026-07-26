page_multifva <- function() {
  layout_sidebar(
    sidebar = sidebar(
      title = "Filter Controls",
      radioButtons("multifva_method", "FVA Calculation Method:", choices = c("Loopless FVA" = "loopless", "Standard FVA" = "standard"), selected = "loopless"),
      checkboxGroupInput("multifva_fractions", "PHB Optimality Fractions:", choices = c("100% of PHB optimum" = "100", "90% of PHB optimum" = "90", "50% of PHB optimum" = "50"), selected = c("100", "90", "50")),
      sliderInput("multifva_max_distance", "Max Distance from PHB Synthesis:", min = 0, max = 4, value = 4, step = 1),
      sliderInput("multifva_min_range", "Minimum Flux Range (mmol/gDW/h):", min = 0, max = 40, value = 0, step = 1),
      searchInput("multifva_search", "Search Reaction / Subsystem:", placeholder = "e.g., ACACT", btnSearch = icon("search"), btnReset = icon("remove"), width = "100%")
    ),
    layout_column_wrap(
      width = 1,
      card(full_screen = TRUE, card_header("PHB Pathway Flux Ranges Across PHB-Optimality Fractions"), plotlyOutput("multifva_bar_plot", height = "400px")),
      card(full_screen = TRUE, card_header("Multi-FVA Detailed Results"), DTOutput("multifva_table"))
    )
  )
}
