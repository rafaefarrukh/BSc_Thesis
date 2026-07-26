###################################################################################################
# UI: Strain Design
###################################################################################################

page_strain_design <- function() {
  
  tagList(
    navset_card_tab(
      
      title = "Strain Design",
      
      # TAB 1: ESSENTIAL REACTIONS (Step 7.16) ########################
      nav_panel(
        title = "Essential Reactions",
        
        layout_columns(
          col_widths = c(4, 8),
          
          card(
            card_header("Knockout Classification"),
            p(class = "small text-muted", "Growth fraction relative to wild type:"),
            tags$ul(
              class = "small",
              tags$li(tags$b("Essential"), " \u2014 \u2264 1% of wild-type growth"),
              tags$li(tags$b("Critical"), " \u2014 1\u201350% of wild-type growth"),
              tags$li(tags$b("Non-essential"), " \u2014 \u2265 50% of wild-type growth")
            ),
            hr(),
            DTOutput("strain_essential_summary_table")
          ),
          
          card(
            card_header("Top 30 Essential / Critical Knockouts by Resulting Growth Fraction"),
            plotlyOutput("strain_essential_plot", height = "560px")
          )
        )
      ),
      
      # TAB 2: GROWTH-COUPLED KNOCKOUTS (Step 7.17) ###################
      nav_panel(
        title = "Growth-Coupled Knockouts",
        
        card(
          class = "border-danger",
          card_header(
            class = "bg-danger-subtle",
            tags$strong("Result: No Growth-Coupling Achieved")
          ),
          p(
            "Every single- and double-reaction knockout candidate examined yielded a forced ",
            "minimum PHB flux of exactly ", tags$b("0.0000 mmol/gDW/h"), " at maximum growth ",
            "\u2014 identical to the wild-type baseline. ", tags$b("No single- or double-reaction "),
            tags$b("knockout was found that couples growth to obligate PHB production"),
            " in this model."
          ),
          p(
            class = "mb-0 small text-muted",
            "This negative result is consistent with the production envelope, which showed a ",
            "minimum feasible PHB flux of zero at every tested growth rate: sufficiently many ",
            "alternate routes exist to balance acetyl-CoA/redox flux without invoking PHB ",
            "synthesis that blocking one or two reactions cannot eliminate the ",
            "\u201Cgrow-without-producing\u201D phenotype."
          )
        ),
        
        card(
          card_header("Ranked Knockout Candidates"),
          radioButtons(
            "strain_knockout_type",
            label = NULL,
            choices = c("Single Knockouts" = "single", "Double Knockouts" = "double"),
            selected = "single",
            inline = TRUE
          ),
          DTOutput("strain_knockout_table")
        )
      )
    )
  )
}

###################################################################################################
