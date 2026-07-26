###################################################################################################
# UI
###################################################################################################

ui <- page_navbar(
  
  title = "BSc Thesis",
  navbar_options = navbar_options(underline = TRUE),
  
  nav_panel(title = "Overview", page_overview()),
  
  nav_panel(title = "Genome Annotation", page_genome()),
  nav_panel(title = "phaC Structure", page_phac()),
  nav_panel(title = "GEM Network", page_gem()),
  nav_panel(title = "FBA", page_fba()),
  nav_panel(title = "Flexibility", page_flexibility()),
  
  nav_menu(
    title = "Sensitivity & Perturbations",
    nav_panel("Environmental Screening", page_environmental_screening()),
    nav_panel("Robustness", page_robustness()),
    nav_panel("PhPP", page_phpp()),
    nav_panel("Dynamic FBA", page_dfba()),
  ),
  
  nav_menu(
    title = "System Diagnostics",
    nav_panel("Diagnostics", page_diagnostics()),
    nav_panel("Bottlenecks", page_bottlenecks()),
    nav_panel("Multi-FVA", page_multifva()),
    nav_panel("Strain Design", page_strain_design()),
  ),
  
  nav_panel(title = "Cost Analysis", page_cost()),
  
  nav_spacer()
)

###################################################################################################