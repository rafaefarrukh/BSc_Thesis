###################################################################################################
# UI: PHB Pathway Flexibility
###################################################################################################

page_flexibility <- function() {
  layout_sidebar(
    sidebar = sidebar(
      title = "Filter Controls",
      width = 320,
      
      sliderInput(
        inputId = "flexibility_max_distance",
        label   = "Max Distance from PHB Synthesis:",
        min     = 0,
        max     = 4,
        value   = 4,
        step    = 1
      ),
      sliderInput(
        inputId = "flexibility_min_range",
        label   = "Minimum Flux Range (mmol/gDW/h):",
        min     = 0,
        max     = 30,
        value   = 0,
        step    = 0.5
      ),
      searchInput(
        inputId     = "flexibility_search",
        label       = "Search Reaction / Subsystem:",
        placeholder = "e.g., ACACT",
        btnSearch   = icon("search"),
        btnReset    = icon("remove"),
        width       = "100%"
      )
    ),
    
    layout_column_wrap(
      width = 1,
      
      # Upper Card: Bar Chart
      card(
        full_screen = TRUE,
        card_header("Feasible Flux Range by Reaction"),
        plotlyOutput("flexibility_bar_plot", height = "400px")
      ),
      
      # Lower Card: Data Table
      card(
        full_screen = TRUE,
        card_header("Pathway Reactions & Flexibility Details"),
        DTOutput("flexibility_table")
      )
    )
  )
}

###################################################################################################