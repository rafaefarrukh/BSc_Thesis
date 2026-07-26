###################################################################################################
# UI: FBA & Production Envelope
###################################################################################################

page_fba <- function() {
  layout_columns(
    col_widths = c(6, 6),
    
    # Left Column: 2 stacked cards
    layout_column_wrap(
      width = 1,
      
      # Upper Card: FBA Scenarios Description
      card(
        card_header("FBA Scenarios Overview"),
        DTOutput("fba_scenarios_table")
      ),
      
      # Bottom Card: Production Envelope
      card(
        full_screen = TRUE,
        card_header("PHB Production Envelope"),
        plotlyOutput("fba_production_envelope_plot")
      )
    ),
    
    # Right Column: 1 card containing pFBA table
    card(
      full_screen = TRUE,
      card_header("pFBA Flux Distribution"),
      DTOutput("fba_pfba_table")
    )
  )
}

###################################################################################################