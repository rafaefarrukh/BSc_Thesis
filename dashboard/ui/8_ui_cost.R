###################################################################################################
# UI: Cost
###################################################################################################

page_cost <- function() {
  
  layout_sidebar(
    sidebar = sidebar(
      title = "Controls & Settings",
      width = 320,
      
      # Medium & Recovery Selection
      selectInput("cost_select_medium",   "Medium Name",      choices = NULL),
      selectInput("cost_select_recovery", "Recovery Method", choices = NULL),
      hr(),
      
      # Table Actions
      tags$label("Add Table Entries", class = "fw-bold mb-2"),
      div(
        class = "d-grid gap-2 mb-3",
        actionButton("cost_add_component_row", "Add Component",       icon = icon("plus"), class = "btn-sm btn-outline-primary"),
        actionButton("cost_add_medium_col",     "Add Medium",          icon = icon("plus"), class = "btn-sm btn-outline-primary"),
        actionButton("cost_add_recovery_row",   "Add Recovery Method", icon = icon("plus"), class = "btn-sm btn-outline-primary")
      ),
      hr(),
      
      # Import Data
      tags$label("Import Data (CSV)", class = "fw-bold mb-2"),
      fileInput("cost_import_media",    "Import Media CSV",      accept = c(".csv")),
      fileInput("cost_import_recovery", "Import Downstream CSV", accept = c(".csv")),
      fileInput("cost_import_yield",    "Import PHB Yield CSV",  accept = c(".csv")),
      hr(),
      
      # Export Data
      tags$label("Export Data", class = "fw-bold mb-2"),
      div(
        class = "d-grid gap-2",
        downloadButton("cost_export_media",    "Export Media CSV",      class = "btn-sm"),
        downloadButton("cost_export_recovery", "Export Recovery CSV",   class = "btn-sm"),
        downloadButton("cost_export_yield",    "Export PHB Yield CSV",  class = "btn-sm")
      )
    ),
    
    # Main Area
    layout_columns(
      col_widths = c(12),
      
      # Tables
      card(
        card_header("Cost & Production Data"),
        navset_card_tab(
          
          nav_panel(
            title = "Overview",
            card(
              card_header("Cost Calculation Overview"),
              markdown("
              ### Medium Composition Cost
              * Describes the cost for various media.
              * All costs are in **PKR/kg**.
              * The composition of each medium is in **g/L**.
              * New components and media can be added using the sidebar.
              
              ---
              
              ### Downstream Recovery Cost
              * Describes the cost for various downstream recovery methods.
              * All costs are in **PKR/1L medium processed**.
              * New methods can be added using the sidebar.
              
              ---
              
              ### PHB Yield
              * Describes the PHB yield (**g/L**) from various media.
              * Yields can be modified directly in the table.
              
              ---
              
              ### Import & Export
              * Media, Downstream Recovery Methods, and PHB Yields can be exported and imported as CSV files for ease.
                       ")
            )
          ),
          
          nav_panel(
            title = "Medium Composition Cost",
            DTOutput("cost_medium_table")
          ),
          
          nav_panel(
            title = "Downstream Recovery Cost",
            DTOutput("cost_recovery_table")
          ),
          
          nav_panel(
            title = "PHB Yield",
            DTOutput("cost_yield_table")
          ),
          
          nav_panel(
            title = "Cost Comparison",
            layout_sidebar(
              sidebar = sidebar(
                title = "Filter Media",
                position = "right",
                width = 220,
                open = TRUE,
                checkboxGroupInput(
                  "cost_select_plot_media",
                  "Display Media:",
                  choices = NULL
                ),
                div(
                  class = "d-grid gap-2",
                  actionButton("cost_plot_select_all",   "Select All",   class = "btn-sm btn-outline-secondary"),
                  actionButton("cost_plot_deselect_all", "Deselect All", class = "btn-sm btn-outline-secondary")
                )
              ),
              plotlyOutput("cost_comparison_plot")
            )
          )
          
        )
      ),
      
      # Summary Value Boxes
      card(
        card_header("Estimated PHB Production Cost"),
        layout_column_wrap(
          width = 1/3,
          value_box(
            title = "Medium Cost",
            value = textOutput("cost_medium_cost_out"),
            showcase = icon("flask")
          ),
          value_box(
            title = "Downstream Recovery Cost",
            value = textOutput("cost_recovery_cost_out"),
            showcase = icon("filter")
          ),
          value_box(
            title = "Total Cost (per kg PHB)",
            value = textOutput("cost_total_cost_out"),
            showcase = icon("coins")
          )
        )
      )
    )
  )
  
}

###################################################################################################