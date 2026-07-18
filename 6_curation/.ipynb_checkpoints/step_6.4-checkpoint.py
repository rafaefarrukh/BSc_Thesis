import cobra
from cobra.flux_analysis import pfba, flux_variability_analysis

# 1. Load the metabolic model (replace with your actual file)
model = cobra.io.read_sbml_model("../models/model_4.xml")

# 2. Define the list of target reactions
target_reactions = [
    "13PPDH","13PPDH2","4ABUTD","5DGLCNR","5DKGR","ACKr","ACONT","ACONTa","ACONTb","ACTD",
    "ACTD_1","ACTDa","ADAPAT","AH6PI","ALATA_D","ALATA_L","ALAtex","ALCD19","ALCD19y",
    "ALCD2ir","ALCD2x","ALCD2y","ASPTA","ASPtpp","BUTKr","CO2t","CO2tex","CO2tpp",
    "COBALT2tex","COBALTt5","DAPDA","DURIPP","DURIPP_1","ECOAH1","ECOAH1_1","ECOAH3",
    "ETOHtex","ETOHtrpp","FE2tex","FE3PYOVDL2","FFSD1r","G6PDH2r","G6PI","GALKr",
    "GALpts","GAPD","GDH1","GK1","GK2","GLUDxi","GLUDy","GLUDy_1","GLUR","GLYBtex",
    "GLYCtex","GLYCtpp","GLYOX_1","GTHOr","H2Ot","H2Otpp","HACD1","HACD1_1","HACD1i",
    "HACD3","HACD8","HACD8i","HG2tex","HSDxi","HSTPTr","Htex","INOSR","Kt1","LEUTA",
    "LEUtex","LGTHL","LLEUDr","METB1","MNtex","MTAM_1","MTHFC","MTHFD","MYCTR",
    "NADDPp_1","NADKd","NADTRHD","NDPK8","NH3c","NH4t","NH4tex","NH4tpp","OCBT_1",
    "ORNTA","P5CR","P5CRx","PGK","PGM","PGMT","PGMT_2","PGMT_B","PRO1x","PTHPS_1",
    "PYNP1","PYNP1_1","PYNP2r","RECOAH3","RHACOAR80","RU5PS","SDPTA","SHSL2","SHSL2r",
    "SMIA1","SMIB1","SSCOARy","SUCOAS","TARTRtpp","THPAT","UAGPT2_1","VALDHr","VALTA",
    "VALtex","VPAMTr","YUMPS","r2465_1","EX_23dhbzs3_e","EX_4abut_e","EX_cit_e",
    "EX_cmcbtt_e","EX_cpgn_un_e","EX_enter_e","EX_fe2_e","EX_fe3_e","EX_fe3pyovd_kt_e",
    "EX_feoxam_un_e","EX_glyald_e","EX_gthrd_e","EX_h2o_e","EX_h_e","EX_istfrnA_e",
    "EX_istfrnB_e","EX_mcbtt_e","EX_nmn_e","EX_o2_e","EX_salchs2_e","EX_salchs4_e",
    "EX_stfrnA_e","EX_stfrnB_e"
]

# 3. Apply physiological bounds to target reactions
# Default bounds of +/- 1000 often cause unbound flux loops.
for rxn_id in target_reactions:
    if rxn_id in model.reactions:
        reaction = model.reactions.get_by_id(rxn_id)
        # Constrain to a realistic biological maximum (e.g., +/- 50 mmol/gDW/h)
        new_lower = max(reaction.lower_bound, -50.0)
        new_upper = min(reaction.upper_bound, 50.0)
        reaction.bounds = (new_lower, new_upper)

# 4. Filter list for reactions that actually exist in the model
valid_reactions = [rxn for rxn in target_reactions if rxn in model.reactions]

# 5. Run loopless Flux Variability Analysis (FVA)
# The loopless=True flag ensures thermodynamically infeasible loops are ignored.
fva_results = flux_variability_analysis(model, valid_reactions, loopless=True)
print("FVA Results (Loopless):")
print(fva_results)

# 6. Run Parsimonious FBA (pFBA)
# pFBA minimizes the sum of all absolute fluxes, eliminating futile cycles.
try:
    pfba_solution = pfba(model)
    print("\npFBA Solution Objective Value:", pfba_solution.objective_value)
    
    # Extract resolved fluxes for the target list
    print("\nResolved Fluxes for Target Reactions (pFBA):")
    for rxn_id in valid_reactions:
        print(f"{rxn_id}: {pfba_solution.fluxes[rxn_id]}")
except Exception as e:
    print("\npFBA failed:", e)