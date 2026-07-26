###################################################################################################
# PREP: GEM Network Data
###################################################################################################

library(tidyverse)
library(jsonlite)

# CONFIG #########################################

model_path     <- "models/model_5.json"
distances_path <- "7_analyses/step_7_results/step_7.11_phb_pathway_distances.csv"
out_path       <- "dashboard/data/gem_network.rds"

# metabolites so ubiquitous ("currency" metabolites) that including them in every
# reaction they touch turns the graph into an unreadable hairball. Flagged (not
# dropped) so the app can toggle them on/off.
currency_metabolites <- c(
  "h_c", "h_e", "h_p", "h2o_c", "h2o_e", "h2o_p",
  "atp_c", "adp_c", "amp_c", "pi_c", "ppi_c",
  "nad_c", "nadh_c", "nadp_c", "nadph_c",
  "coa_c", "co2_c", "o2_c", "nh4_c", "na1_c",
  "ACP_c"
)

# LOAD ###########################################

model     <- fromJSON(model_path, simplifyVector = FALSE)
distances <- read_csv(distances_path, show_col_types = FALSE)

reactions   <- model$reactions
metabolites <- model$metabolites

# lookups by id
rxn_by_id <- set_names(reactions, map_chr(reactions, "id"))
met_by_id <- set_names(metabolites, map_chr(metabolites, "id"))

# REACTION NODE TABLE (PHB preset scope only) ####

phb_reaction_ids <- distances$reaction_id

reaction_nodes <- distances %>%
  transmute(
    id         = reaction_id,
    node_id    = paste0("rxn_", reaction_id),
    label      = reaction_id,
    name       = if_else(is.na(reaction_name) | reaction_name == "", reaction_id, reaction_name),
    subsystem  = subsystem,
    distance   = as.integer(distance_from_PHB),
    node_type  = "reaction"
  )

# pull bounds / gpr / equation string from the model for each reaction in scope
safe_scalar <- function(x, default) if (is.null(x)) default else x

reaction_meta <- map_dfr(phb_reaction_ids, function(rid) {
  rxn <- rxn_by_id[[rid]]
  if (is.null(rxn)) {
    return(tibble(
      id = rid, lower_bound = NA_real_, upper_bound = NA_real_,
      gene_reaction_rule = NA_character_, equation = NA_character_,
      n_metabolites = 0L
    ))
  }
  mets <- rxn$metabolites
  coefs <- unlist(mets)
  substrate_ids <- names(mets)[coefs < 0]
  product_ids   <- names(mets)[coefs > 0]
  eqn <- paste(
    paste(substrate_ids, collapse = " + "),
    "\u2192",
    paste(product_ids, collapse = " + ")
  )
  tibble(
    id = rid,
    lower_bound = as.numeric(safe_scalar(rxn$lower_bound, NA_real_)),
    upper_bound = as.numeric(safe_scalar(rxn$upper_bound, NA_real_)),
    gene_reaction_rule = as.character(safe_scalar(rxn$gene_reaction_rule, "")),
    equation = eqn,
    n_metabolites = length(mets)
  )
})

reaction_nodes <- reaction_nodes %>% left_join(reaction_meta, by = "id")

# METABOLITE NODE TABLE (only mets touched by in-scope reactions) ####

met_ids_in_scope <- reactions %>%
  keep(~ .x$id %in% phb_reaction_ids) %>%
  map(~ names(.x$metabolites)) %>%
  unlist() %>%
  unique()

metabolite_nodes <- map_dfr(met_ids_in_scope, function(mid) {
  met <- met_by_id[[mid]]
  met_name <- safe_scalar(met$name, mid)
  if (is.null(met_name) || met_name == "") met_name <- mid
  tibble(
    id          = mid,
    node_id     = paste0("met_", mid),
    label       = met_name,
    name        = met_name,
    compartment = as.character(safe_scalar(met$compartment, NA_character_)),
    formula     = as.character(safe_scalar(met$formula, NA_character_)),
    node_type   = "metabolite",
    is_currency = mid %in% currency_metabolites
  )
})

# EDGE TABLE (reaction <-> metabolite, bipartite) ####

edges <- reactions %>%
  keep(~ .x$id %in% phb_reaction_ids) %>%
  map_dfr(function(rxn) {
    mets <- rxn$metabolites
    tibble(
      reaction_id   = rxn$id,
      metabolite_id = names(mets),
      coefficient   = unlist(mets),
      role          = if_else(unlist(mets) < 0, "substrate", "product")
    )
  }) %>%
  mutate(
    from = if_else(role == "substrate", paste0("met_", metabolite_id), paste0("rxn_", reaction_id)),
    to   = if_else(role == "substrate", paste0("rxn_", reaction_id), paste0("met_", metabolite_id))
  )

# SAVE ###########################################

gem_network <- list(
  reaction_nodes   = reaction_nodes,
  metabolite_nodes = metabolite_nodes,
  edges            = edges,
  distance_range   = range(reaction_nodes$distance, na.rm = TRUE)
)

write_rds(gem_network, out_path, compress = "xz")

cat(sprintf(
  "Saved %s: %d reaction nodes, %d metabolite nodes (%d currency), %d edges. Distance range: %d-%d\n",
  out_path,
  nrow(reaction_nodes),
  nrow(metabolite_nodes),
  sum(metabolite_nodes$is_currency),
  nrow(edges),
  gem_network$distance_range[1],
  gem_network$distance_range[2]
))

###################################################################################################
