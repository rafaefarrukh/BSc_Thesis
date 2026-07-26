page_phpp <- function() {
  layout_sidebar(
    sidebar = sidebar(
      title = "PhPP Controls",
      selectInput(
        "phpp_carbon_source", "Carbon Source:",
        choices  = c("Glucose", "Maltose"),
        selected = "Glucose"
      ),
      card(
        class = "bg-light mt-3",
        card_header("PhPP Guide"),
        p(
          class = "small",
          "Phenotype Phase Planes map optimal growth rate across a grid of carbon-source and ",
          "oxygen uptake fluxes. Markers indicate candidate phase-transition points identified ",
          "via curvature-based detection (top of the growth surface's second-derivative magnitude)."
        )
      )
    ),
    layout_column_wrap(
      width = 1,
      card(full_screen = TRUE, card_header("Phenotype Phase Plane"), plotlyOutput("phpp_contour_plot", height = "500px")),
      card(card_header("Candidate Phase Transition Points"), DTOutput("phpp_region_table"))
    )
  )
}
