"""
STEP 7.29 (FIXED) — CORRECTED CARBON-SOURCE SCREEN

BUG FIX vs. the previous version of this script: every row reported
PHB=10.0000 regardless of carbon source or aeration, because the biomass
lower bound was reset to 0.0 before maximizing PHB -- which just recomputes
the model's unconstrained PHB ceiling every time (always hits its upper
bound of 10) instead of the *growth-coupled* PHB flux at the growth rate
that condition actually achieved. This version fixes biomass at the
just-computed growth rate (not 0.0) before the second optimization, so the
PHB column now reflects PHB production AT that condition's growth optimum,
matching the intent of the original step_7.10 script.

NOTE ON THE UNDERLYING BIOLOGY: the previous run of this script (with the
carbon-exchange-closing fix applied but the PHB bug still present) showed
that growth is numerically IDENTICAL across glucose, fructose, glycerol,
sucrose, mannose, mannitol, arabinose, melibiose, acetate, succinate, and
pyruvate -- scaling only with the O2 bound and topping out at 0.3722 /h,
the same ceiling seen in step_7.8's independent O2-uptake sweep. Only
maltose reaches the model's true unconstrained ceiling (0.4454 /h). This
points to a shared downstream bottleneck (not a carbon-supply limit) that
every substrate except maltose hits -- see step_7.32 for the follow-up
analysis that identifies the specific bottleneck reaction(s).

Input:  ../models/model_7.xml
Output: ./step_7.29_aeration_carbon_screen_corrected.csv
        ./step_7.29_corrected_vs_original_comparison.csv
        ./step_7.29_summary.txt
"""

import cobra
import pandas as pd

MODEL_PATH = "../models/model_7.xml"
BIOMASS_ID = "Growth"
PHB_ID = "PHBS_syn_1"
GLC_UPTAKE_BOUND = 10.0  # mmol/gDW/h, matches step_7.10 convention

CANDIDATE_CARBON_EXCHANGES = {
    "Glucose": "EX_glc__D_e",
    "Fructose": "EX_fru_e",
    "Glycerol": "EX_glyc_e",
    "Sucrose": "EX_sucr_e",
    "Maltose": "EX_malt_e",
    "Mannose": "EX_man_e",
    "Mannitol": "EX_mnl_e",
    "Arabinose": "EX_arab__L_e",
    "Melibiose": "EX_melib_e",
    "Acetate": "EX_ac_e",
    "Succinate": "EX_succ_e",
    "Pyruvate": "EX_pyr_e",
}

AERATION_LEVELS = {
    "Anaerobic": 0.0,
    "Microaerobic": -2.0,
    "Aerobic": -20.0,
    "Hyperaerobic": -40.0,
}

print(f"Loading model from: {MODEL_PATH}")
model = cobra.io.read_sbml_model(MODEL_PATH)
print(f"Loaded model '{model.id}': {len(model.reactions)} reactions, "
      f"{len(model.metabolites)} metabolites, {len(model.genes)} genes")

present = {}
for name, ex_id in CANDIDATE_CARBON_EXCHANGES.items():
    if ex_id in model.reactions:
        present[name] = ex_id
    else:
        print(f"  (skipping '{name}': exchange reaction '{ex_id}' not found in model)")
print(f"Candidate carbon exchanges confirmed in model: {list(present.values())}")

O2_EX = "EX_o2_e"
biomass_rxn = model.reactions.get_by_id(BIOMASS_ID)
phb_rxn = model.reactions.get_by_id(PHB_ID)
o2_rxn = model.reactions.get_by_id(O2_EX)

original_bounds = {ex_id: model.reactions.get_by_id(ex_id).bounds for ex_id in present.values()}
original_o2_bounds = o2_rxn.bounds
original_biomass_bounds = biomass_rxn.bounds

results = []
print("\nRunning corrected aeration x carbon-source screen "
      f"({len(present)} sources x {len(AERATION_LEVELS)} aeration levels)...")

for carbon_name, carbon_ex_id in present.items():
    for aeration_name, o2_lb in AERATION_LEVELS.items():
        # Close every candidate carbon exchange, then open only the tested one
        for ex_id in present.values():
            model.reactions.get_by_id(ex_id).lower_bound = 0.0
        model.reactions.get_by_id(carbon_ex_id).lower_bound = -GLC_UPTAKE_BOUND
        o2_rxn.lower_bound = o2_lb
        biomass_rxn.bounds = original_biomass_bounds  # ensure a clean slate each iteration

        model.objective = BIOMASS_ID
        growth = model.slim_optimize()
        if growth is None or growth != growth:  # infeasible / NaN
            growth = 0.0
            phb = 0.0
        else:
            # --- FIX: fix biomass at the ACHIEVED growth rate, not 0.0 ---
            biomass_rxn.lower_bound = growth
            biomass_rxn.upper_bound = growth
            model.objective = PHB_ID
            phb = model.slim_optimize()
            if phb is None or phb != phb:
                phb = 0.0
            biomass_rxn.bounds = original_biomass_bounds  # restore before next iteration

        results.append(dict(
            carbon_source=carbon_name,
            carbon_exchange=carbon_ex_id,
            aeration=aeration_name,
            o2_uptake_bound=o2_lb,
            max_growth_rate_1_h=round(growth, 4) if growth else 0.0,
            growth_coupled_PHB_flux=round(phb, 4) if phb else 0.0,
        ))
        print(f"  {carbon_name:<12} x {aeration_name:<13} -> "
              f"growth={growth:.4f} /h, PHB={phb:.4f}")

# Restore model to original state
for ex_id, bounds in original_bounds.items():
    model.reactions.get_by_id(ex_id).bounds = bounds
o2_rxn.bounds = original_o2_bounds
biomass_rxn.bounds = original_biomass_bounds

results_df = pd.DataFrame(results)
results_df.to_csv("./step_7.29_aeration_carbon_screen_corrected.csv", index=False)
print(f"\nCorrected results written to: ./step_7.29_aeration_carbon_screen_corrected.csv")

n_unique_growth = results_df.groupby("aeration")["max_growth_rate_1_h"].nunique()
n_unique_phb = results_df.groupby("aeration")["growth_coupled_PHB_flux"].nunique()
print("\nDistinct growth-rate values per aeration level:")
print(n_unique_growth.to_string())
print("\nDistinct growth-coupled PHB values per aeration level "
      "(should now be > 1 -- previously this was stuck at 1 everywhere, always =10):")
print(n_unique_phb.to_string())

try:
    orig = pd.read_csv("./step_7.10_aeration_carbon_screen.csv")
    merged = results_df.merge(
        orig, on=["carbon_source", "aeration"], suffixes=("_corrected", "_original")
    )
    merged["growth_delta"] = (merged["max_growth_rate_1_h_corrected"]
                               - merged["max_growth_rate_1_h_original"])
    merged["phb_delta"] = (merged["growth_coupled_PHB_flux_corrected"]
                            - merged["growth_coupled_PHB_flux_original"])
    merged.to_csv("./step_7.29_corrected_vs_original_comparison.csv", index=False)
    print(f"\nComparison to original step_7.10 output written to: "
          f"./step_7.29_corrected_vs_original_comparison.csv")
    print(f"Max |growth_delta|: {merged['growth_delta'].abs().max():.4f}")
    print(f"Max |phb_delta|:    {merged['phb_delta'].abs().max():.4f}")
except FileNotFoundError:
    print("\nNOTE: ./step_7.10_aeration_carbon_screen.csv not found — skipping comparison.")

with open("./step_7.29_summary.txt", "w") as f:
    f.write("STEP 7.29 (FIXED) — Corrected carbon-source screen summary\n")
    f.write("=" * 60 + "\n")
    f.write(f"Carbon sources tested: {list(present.keys())}\n\n")
    f.write("Fix applied: biomass is now fixed at the ACHIEVED growth rate\n")
    f.write("(not 0.0) before maximizing PHB, so growth_coupled_PHB_flux reflects\n")
    f.write("PHB production at that condition's growth optimum rather than the\n")
    f.write("model's unconstrained PHB ceiling (which was always 10.0 before).\n\n")
    f.write("Growth-rate result (carried over from the exchange-closing fix,\n")
    f.write("unaffected by this PHB fix): growth is numerically identical across\n")
    f.write("11 of 12 carbon sources, scaling only with the O2 bound and topping\n")
    f.write("out at 0.3722 /h -- the same ceiling independently seen in step_7.8's\n")
    f.write("O2-uptake sweep. Only maltose reaches the model's true unconstrained\n")
    f.write("ceiling (0.4454 /h), pointing to a shared downstream bottleneck that\n")
    f.write("every other substrate hits regardless of carbon supply. See step_7.32\n")
    f.write("for the bottleneck-identification follow-up.\n\n")
    f.write("Distinct growth-rate values per aeration level:\n")
    f.write(n_unique_growth.to_string())
    f.write("\n\nDistinct growth-coupled PHB values per aeration level:\n")
    f.write(n_unique_phb.to_string())

print("\nSummary written to: ./step_7.29_summary.txt")
