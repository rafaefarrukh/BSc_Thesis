page_dfba <- function() {
  layout_sidebar(
    sidebar = sidebar(
      title = "dFBA Controls",
      radioButtons(
        "dfba_carbon", "Carbon Source:",
        choices = c("Glucose" = "glucose", "Maltose" = "maltose"),
        selected = "glucose"
      ),
      radioButtons(
        "dfba_condition", "Feeding Scenario:",
        choices = c("Standard Batch" = "standard", "Nitrogen-Limited Batch" = "nlimited", "Two-Stage Fed-Batch" = "twostage"),
        selected = "standard"
      )
    ),
    layout_column_wrap(
      width = 1,
      card(
        class = "bg-warning-subtle border-warning mb-3",
        card_body(
          h5("Physiological Boundary Warning"),
          p(
            "Dynamic integration past ~80 hours exhibits non-physiological accumulation artifacts ",
            "in the Standard and Nitrogen-Limited scenarios, due to substrate exhaustion combined ",
            "with explicit-Euler integration under a purely growth-maximizing objective. The ",
            "Two-Stage scenario switches its objective to PHB maximization at t = 24h and does not ",
            "exhibit this artifact."
          )
        )
      ),
      card(full_screen = TRUE, card_header("Dynamic Trajectories"), plotlyOutput("dfba_timeseries_plot", height = "500px")),
      card(card_header("Simulation Data"), DTOutput("dfba_table"))
    )
  )
}
