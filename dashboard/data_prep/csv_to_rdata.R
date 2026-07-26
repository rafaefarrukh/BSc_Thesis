###################################################################################################
# DATA PREPARATION: COMPILE CSV INTO RDATA
###################################################################################################

# import all csv files
csv_list <- lapply(list.files(path = "7_analyses/step_7_results", pattern = "\\.csv$", full.names = TRUE), read.csv)

# rename all dataframes
names(csv_list) <- paste0("df_", seq_along(csv_list))

# rename to match with server req -----------------------------------------------------------------

# FBA / pFBA / Production Envelope (page 4)
names(csv_list)[1]  <- "fba_results"
names(csv_list)[2]  <- "pfba_flux_distribution"
names(csv_list)[3]  <- "production_envelope"

# Phenotype Phase Plane (page 6.C)
names(csv_list)[4]  <- "phpp_glucose_grid"
names(csv_list)[5]  <- "phpp_glucose_transitions"
names(csv_list)[6]  <- "phpp_maltose_grid"
names(csv_list)[7]  <- "phpp_maltose_transitions"

# Dynamic FBA (page 6.D)
names(csv_list)[8]  <- "dfba_maltose_nlimited"
names(csv_list)[9]  <- "dfba_glucose_nlimited"
names(csv_list)[10] <- "dfba_maltose_standard"
names(csv_list)[11] <- "dfba_glucose_standard"
names(csv_list)[12] <- "dfba_maltose_twostage"
names(csv_list)[13] <- "dfba_glucose_twostage"

# Environmental Screening - medium screen (page 6.A)
names(csv_list)[14] <- "medium_screen_results"
names(csv_list)[15] <- "minimal_medium_exchanges"

# Robustness - 1D & 2D sweeps (page 6.B)
names(csv_list)[16] <- "robust_1d_glucose"
names(csv_list)[17] <- "robust_1d_maltose"
names(csv_list)[18] <- "robust_1d_nitrogen"
names(csv_list)[19] <- "robust_1d_oxygen"
names(csv_list)[20] <- "robust_2d_glc_nh4"
names(csv_list)[21] <- "robust_2d_malt_nh4"
names(csv_list)[22] <- "robust_2d_glc_o2"

# Environmental Screening - aeration x carbon source (page 6.A)
names(csv_list)[23] <- "aeration_carbon_results"

# System Diagnostics - shadow prices & reduced costs (page 7.A)
names(csv_list)[24] <- "diag_reduced_cost_growth"
names(csv_list)[25] <- "diag_reduced_cost_phb"
names(csv_list)[26] <- "diag_shadow_price_growth"
names(csv_list)[27] <- "diag_shadow_price_phb"

# Full PHB pathway BFS distance table (source table for GEM network page; not directly
# consumed by any 6.A-7.D page, but renamed here for clarity/consistency)
names(csv_list)[28] <- "phb_pathway_distance"

# Pathway Bottleneck Analysis (page 7.B)
names(csv_list)[29] <- "bottleneck_results"

# PHB Pathway Flexibility (page 5) - also reused as reaction metadata source for Multi-FVA (7.C)
names(csv_list)[30] <- "pathway_flexibility"

# Flux Shift Analysis (growth-optimal vs PHB-optimal states); not consumed by any current
# page, renamed for clarity in case a future page is added
names(csv_list)[31] <- "flux_shift_results"

# Multi-fraction / Loopless FVA (page 7.C)
# Fraction/loopless assignment confirmed against multifva_summary's n_blocked counts:
#   100%: loopless -> 103 blocked, standard -> 92 blocked
#    90%: loopless ->  13 blocked, standard ->  9 blocked
#    50%: loopless ->  13 blocked, standard ->  9 blocked
names(csv_list)[32] <- "multifva_100_loopless"
names(csv_list)[33] <- "multifva_100_standard"
names(csv_list)[34] <- "multifva_50_loopless"
names(csv_list)[35] <- "multifva_50_standard"
names(csv_list)[36] <- "multifva_90_loopless"
names(csv_list)[37] <- "multifva_90_standard"
names(csv_list)[38] <- "multifva_summary"

# Strain Design / Knockout Screening (page 7.D)
names(csv_list)[39] <- "knockout_essentiality"
names(csv_list)[40] <- "knockout_double"
names(csv_list)[41] <- "knockout_single"
names(csv_list)[42] <- "knockout_cutsets"

# convert to data.frames
csv_list <- lapply(csv_list, as.data.frame)

# add to global envr
list2env(csv_list, envir = .GlobalEnv)

# export
save(list = names(csv_list), file = "dashboard/data/step_7_results.RData")

###################################################################################################
