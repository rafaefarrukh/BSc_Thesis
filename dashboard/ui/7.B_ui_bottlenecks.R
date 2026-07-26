page_bottlenecks <- function() {
    nav_panel("Pathway Bottlenecks", layout_column_wrap(width = 1, card(full_screen = TRUE, card_header("Ranked Bottlenecks"), plotlyOutput("diag_bottleneck_plot", height = "400px")), card(card_header("Bottleneck Details"), DTOutput("diag_bottleneck_table"))))
}
