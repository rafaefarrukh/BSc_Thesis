###################################################################################################
# Step 6.1: MODEL MODIFICATIONS
###################################################################################################

# libraries
library(tidyverse)
library(ape)
library(cobrar)

# work dir
getwd()
# should be ../GEM

# import model 1
model <- readSBMLmod("models/model_1.xml")

# EXTRACT INFO ##################################

# metabolites
metabolites <- data.frame(
  id = model@met_id,
  name = model@met_name,
  compartment = model@met_comp
)

# reactions
reactions <- data.frame(
  id = model@react_id,
  name = model@react_name,
  reac = printReaction(model, 1:react_num(model), T),
  reaction = printReaction(model, 1:react_num(model)),
  lb = model@lowbnd,
  ub = model@uppbnd
)

# fba
fba(model)@obj

# CHECK PHB PRODUCING ENZYMES ###################

View(metabolites[str_detect(metabolites$name, regex("3-hydroxy", ignore_case = TRUE)), ])
# bhb_c & bhb_p & bhb_e: the monomer form of PHB
# 3hbcoa_c: direct CoA-bound monomer, but (R) isomer is required.
# 3hptcoa_c: precursor for 3HV units in copolymers.
# 3hhcoa_c, 3hocoa_c, 3hdcoa_c, 3hddcoa_c, 3htdcoa_c, 3hhdcoa_c, 3hodcoa_c: These are fatty acid beta-oxidation which channel into mcl pathways

View(metabolites[str_detect(metabolites$name, regex("acetoacetyl", ignore_case = TRUE)), ])
# aacoa_c: acetoacetyl-CpA

View(metabolites[str_detect(metabolites$name, regex("acetyl-coA", ignore_case = TRUE)), ])
# accoa_c: acetyl-CoA

View(reactions[str_detect(reactions$reac, regex("bhb", ignore_case = TRUE)), ])
# BDH: converts monomer into acetoacetate
# BHBtpp: moves bhb from periplasm into cytoplasm
# RHA40tex: allows diffusion of bhb between periplasm and extracellular space
# EX_bhb_e: exchange reaction which allows GEM to produce/secrete bhb

View(reactions[str_detect(reactions$reac, regex("aacoa_c", ignore_case = TRUE)), ])
# ACACCT: irreversible shortcut to produce Acetoacetyl-CoA
# ACACT1r: first step of the canonical PHB synthesis pathway (~phaA)
# HACD1 / HACD1_1 (Forward): Reduces aacoa_c to 3hbcoa_c / 3hbycoa_c (~phaB)
# HACD1i: oxidation of 3hbcoa_c back to Acetoacetyl-CoA during polymer mobilization or fatty acid beta-oxidation

View(reactions[str_detect(reactions$reac, regex("accoa_c", ignore_case = TRUE)), ])
# PDH / PDHbr / POR_syn: generate acetyl-CoA
# ACS / ACS_1: converts free acetate into acetyl-CoA (organic acids as substrate)
# CS / ACCOAC / KAS2 / KAS7 / KAS8 / KAS13 / KAS17: competing reactions for acetyl CoA

# missing reactions
## PHB Synthase (PhaC): 3hbcoa_c (or R-enantiomer) --> phb_c+coa_c
## PHB Depolymerase (PhaZ): phb_c+h2o_c --> bhb_c
## Epimerase: 3hbcoa_c  <--> R-3hbcoa_c
## R-specific Reductase: aacoa_c+h_c+nadph_c⇌R-3hbcoa_c+nadp_c

# ADD MISSING INFO ##############################

# missing metabolites
model <- addMetabolite(model,
                     id = c("3hbcoa__R", "phb_c"),
                     comp = c("C_c", "C_c"),
                     name = c("(R)-3-hydroxybutyryl-CoA", "Poly-beta-hydroxybutyrate"),
                     chemicalFormula = c("C25H38N7O18P3S", "C4H6O2"), 
                     charge = c(-4, 0)
)

# phaB
model <- addReact(model, 
                id = "AACOAR_syn", 
                met = c("aacoa_c", "h_c", "nadph_c", "3hbcoa__R", "nadp_c"),
                Scoef = c(-1, -1, -1, 1, 1), 
                lb = -10, 
                ub = 10,
                reactName = "Acetoacetyl-CoA reductase (R-specific)",
)

# phaC
model <- addReact(model, 
                id = "PHBS_syn_1", 
                met = c("3hbcoa__R", "phb_c", "coa_c"),
                Scoef = c(-1, 1, 1), 
                lb = 0, 
                ub = 10,
                reactName = "PHB synthetase"
)

# demand reaction for optimization objective
model <- addReact(model, 
                id = "SK_phb_c", 
                met = "phb_c",
                Scoef = -1, 
                lb = 0, 
                ub = 10,
                reactName = "Demand reaction for PHB accumulation"
)

# SUGAR UTILIZATION #############################

# from a semi-quantitative test, we know that the bacteria utilizes the following sugars
# Glucose > Sucrose > Maltose > Fructose > Mannose > Mannitol
# from QTS24, we we know that the bacteria utilizes the following sugars
# Arabinose, Rhamnose, Melibiose
# from both, we know it doesnt utilize the following sugars
# Sorbitol, Lactose, Adonitol, Raffinose

# checking the non-utilized sugars
View(metabolites[str_detect(metabolites$name, regex("sorbitol|lactose|adonitol|raffinose", ignore_case = TRUE)), ])
met <- metabolites$id[str_detect(metabolites$name, regex("sorbitol|lactose|adonitol|raffinose", ignore_case = TRUE))][9:15]

# checking reactions
View(reactions[str_detect(reactions$reac, regex(paste0(met, collapse = "|"), ignore_case = TRUE)), ])
react <- reactions$id[str_detect(reactions$reac, regex(paste0(met, collapse = "|"), ignore_case = TRUE))]

# remove the reactions and metabolites from the model
model <- rmReact(model, react, TRUE)

# VERIFICATION AND BOUNDS

# bounds

model@lowbnd[model@lowbnd < 0] <- -10
model@lowbnd[model@uppbnd < 0] <- 10

model <- changeBounds(model, react = "EX_o2_e", lb = -20.0, ub = 10)
model <- changeBounds(model, react = "ATPM", lb = 8.39, ub = 10)
model <- changeBounds(model, react = "EX_nh4_e", lb = -15.0, ub = 10)
model <- changeBounds(model, react = "EX_pi_e", lb = -5.0, ub = 10)

model <- changeBounds(model, react = "EX_malt_e", lb = -10.0, ub = 10)
model <- changeBounds(model, react = "EX_sucr_e", lb = -9.0, ub = 10)
model <- changeBounds(model, react = "EX_glc__D_e", lb = -8.0, ub = 10)
model <- changeBounds(model, react = "EX_fru_e", lb = -5.0, ub = 10)
model <- changeBounds(model, react = "EX_man_e", lb = -1.0, ub = 10)
model <- changeBounds(model, react = "EX_mnl_e", lb = -1.0, ub = 10)
model <- changeBounds(model, react = "EX_melib_e", lb = -5.0, ub = 10)
model <- changeBounds(model, react = "EX_arab__L_e", lb = -5.0, ub = 10)

# fba
fba(model)@obj

# metabolites
metabolites <- data.frame(
  id = model@met_id,
  name = model@met_name,
  compartment = model@met_comp
)

# reactions
reactions <- data.frame(
  id = model@react_id,
  name = model@react_name,
  reac = printReaction(model, 1:react_num(model), T),
  reaction = printReaction(model, 1:react_num(model)),
  lb = model@lowbnd,
  ub = model@uppbnd
)

# EXPORT MODEL ##################################

writeSBMLmod(model, "models/model_2.xml")

###################################################################################################