page_diagnostics <- function() {
  nav_panel("Shadow Prices & Reduced Costs", layout_columns(col_widths = c(6, 6), card(card_header("Growth Optimum Profile"), DTOutput("diag_growth_opt_table")), card(card_header("PHB Optimum Profile"), DTOutput("diag_phb_opt_table"))))
}
