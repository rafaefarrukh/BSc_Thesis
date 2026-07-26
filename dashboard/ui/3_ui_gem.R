###################################################################################################
# UI: GEM
###################################################################################################

page_gem <- function() {

  layout_sidebar(

    sidebar = sidebar(
      width = 300,

      h5("Metabolic Network Presets"),
      actionButton(
        inputId = "gem_preset_phb",
        label   = "PHB Synthesis",
        icon    = icon("diagram-project"),
        class   = "btn-primary w-100 mb-2"
      ),
      actionButton(
        inputId = "gem_preset_complete",
        label   = "Complete Network",
        icon    = icon("circle-nodes"),
        class   = "btn-outline-secondary w-100 mb-2"
      ),
      div(
        class = "text-muted small mb-3",
        "The PHB Synthesis preset draws the reactions in the curated PHB pathway distance table up to
         the chosen hop distance, gradient-colored by topological distance from PHB synthase. The
         Complete Network preset shows every reaction in that table regardless of distance."
      ),

      hr(),

      h5("Network Options"),
      sliderInput(
        inputId = "gem_max_distance",
        label   = "Max distance from PHB synthase (hops)",
        min     = 0,
        max     = 9,
        value   = 3,
        step    = 1,
        width   = "100%"
      ),
      checkboxInput(
        inputId = "gem_hide_currency",
        label   = "Hide currency metabolites (ATP, NAD(P)H, H2O, ...)",
        value   = TRUE
      ),
      checkboxInput(
        inputId = "gem_show_metabolites",
        label   = "Show metabolite nodes",
        value   = TRUE
      ),
      actionButton(
        inputId = "gem_toggle_physics",
        label   = "Freeze Layout",
        icon    = icon("lock"),
        class   = "btn-outline-secondary w-100 mt-2"
      )
    ),

    layout_columns(
      col_widths = c(9, 3),

      card(
        card_header("Network"),
        full_screen = TRUE,
        visNetworkOutput("gem_network_plot", height = "650px")
      ),

      card(
        card_header("Selected Node"),
        uiOutput("gem_node_details")
      )

    ),

    layout_columns(
      col_widths = c(12),
      card(
        card_header("Reactions in View"),
        full_screen = TRUE,
        DTOutput("gem_reaction_table")
      )
    )

  )

}

###################################################################################################
