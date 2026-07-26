###################################################################################################
# UI: Environmental Screening
###################################################################################################

page_environmental_screening <- function() {
  
  tagList(
    navset_card_tab(
      
      title = "Environmental Screening",
      
      # TAB 1: MEDIUM SCREENING (Step 7.06) #########################
      nav_panel(
        title = "Medium Screening",
        
        layout_columns(
          col_widths = c(7, 5),
          
          card(
            card_header("Growth Rate & Growth-Coupled PHB Flux Across Candidate Media"),
            plotlyOutput("env_medium_comparison_plot", height = "420px")
          ),
          
          card(
            card_header("Minimal Medium — Required Exchange Reactions"),
            p(
              class = "text-muted small",
              "Minimum set of active exchange reactions required to sustain 50% of ",
              "maximal aerobic, nutrient-replete growth."
            ),
            DTOutput("env_minimal_medium_table")
          )
        ),
        
        card(
          card_header("Medium Screening — Summary Table"),
          DTOutput("env_medium_table")
        )
      ),
      
      # TAB 2: AERATION x CARBON SOURCE (Step 7.09) #################
      nav_panel(
        title = "Aeration \u00D7 Carbon Source",
        
        layout_columns(
          col_widths = c(6, 6),
          
          card(
            card_header("Growth Rate — Aeration \u00D7 Carbon Source Heatmap"),
            plotlyOutput("env_aeration_growth_heatmap", height = "480px")
          ),
          
          card(
            card_header("Growth-Coupled PHB Flux — Aeration \u00D7 Carbon Source Heatmap"),
            plotlyOutput("env_aeration_phb_heatmap", height = "480px")
          )
        ),
        
        card(
          card_header("Condition Matrix — Detail Table"),
          p(
            class = "text-muted small",
            "Full 4 (aeration levels) \u00D7 12 (carbon sources) results grid. ",
            "Select a carbon source below to filter."
          ),
          uiOutput("env_carbon_source_selector"),
          DTOutput("env_aeration_table")
        )
      )
    )
  )
}

###################################################################################################
