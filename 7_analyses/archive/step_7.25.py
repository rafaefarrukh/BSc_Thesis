"""
STEP 7.25 — LOOPLESS RE-VERIFICATION OF KEY FLUX RESULTS

MEMOTE flagged 134 reactions (6.83%) participating in stoichiometrically
balanced cycles (SBCs) and 120 non-blocked reactions (17.2%) capable of
carrying unbounded flux. Several step_7 analyses (7.2, 7.3, 7.13, 7.16, 7.23)
were run with standard (non-loopless) FBA/pFBA, so their absolute flux
magnitudes may be inflated by thermodynamically infeasible internal cycles.

This script re-solves the growth-optimal and PHB-optimal states with loopless
constraints (cobra.flux_analysis.loopless.loopless_solution) and reports how
much each flux value shifts relative to the original (non-loopless) pFBA
results, so you can see exactly which reactions were cycle-inflated.

Input:  ../models/model_7.xml
        ./step_7.3_pfba_flux_distribution.csv   (non-loopless PHB-optimal pFBA)
Output: ./step_7.25_loopless_growth_optimal.csv
        ./step_7.25_loopless_phb_optimal.csv
        ./step_7.25_loopless_vs_standard_delta.csv
        ./step_7.25_summary.txt
"""

import cobra
import pandas as pd
import numpy as np
from cobra.flux_analysis.loopless import loopless_solution
from cobra.flux_analysis import pfba

MODEL_PATH = "../models/model_7.xml"
BIOMASS_ID = "Growth"
PHB_ID = "PHBS_syn_1"
NON_LOOPLESS_PFBA_CSV = "./step_7.3_pfba_flux_distribution.csv"
DELTA_THRESHOLD = 1e-6  # ignore near-zero shifts when reporting "changed" reactions

print(f"Loading model from: {MODEL_PATH}")
model = cobra.io.read_sbml_model(MODEL_PATH)
print(f"Loaded model '{model.id}': {len(model.reactions)} reactions, "
      f"{len(model.metabolites)} metabolites, {len(model.genes)} genes")

biomass_rxn = model.reactions.get_by_id(BIOMASS_ID)
phb_rxn = model.reactions.get_by_id(PHB_ID)
print(f"Biomass reaction: {biomass_rxn.id} | bounds={biomass_rxn.bounds}")
print(f"PHB reaction:     {phb_rxn.id} | bounds={phb_rxn.bounds}")

print("=" * 70)
print("STEP 7.25 — LOOPLESS RE-VERIFICATION")
print("=" * 70)

# ---- Growth-optimal, loopless ----
model.objective = BIOMASS_ID
sol_growth = model.optimize()
print(f"\n[growth_optimum] standard objective value = {sol_growth.objective_value:.6f}")
try:
    ll_growth = loopless_solution(model, fluxes=sol_growth.fluxes)
    growth_df = ll_growth.fluxes.rename("loopless_flux").to_frame()
    growth_df["standard_flux"] = sol_growth.fluxes
    growth_df["delta"] = growth_df["loopless_flux"] - growth_df["standard_flux"]
    growth_df.index.name = "reaction_id"
    growth_df.to_csv("./step_7.25_loopless_growth_optimal.csv")
    n_changed = (growth_df["delta"].abs() > DELTA_THRESHOLD).sum()
    print(f"  loopless growth flux (Growth) = {ll_growth.fluxes[BIOMASS_ID]:.6f}")
    print(f"  reactions with non-trivial loopless correction: {n_changed}")
except Exception as e:
    print(f"  WARNING: loopless_solution failed at growth optimum: {e}")
    growth_df = None

# ---- PHB-optimal, loopless ----
model.objective = PHB_ID
sol_phb = model.optimize()
print(f"\n[PHB_optimum] standard objective value = {sol_phb.objective_value:.6f}")
try:
    ll_phb = loopless_solution(model, fluxes=sol_phb.fluxes)
    phb_df = ll_phb.fluxes.rename("loopless_flux").to_frame()
    phb_df["standard_flux"] = sol_phb.fluxes
    phb_df["delta"] = phb_df["loopless_flux"] - phb_df["standard_flux"]
    phb_df.index.name = "reaction_id"
    phb_df.to_csv("./step_7.25_loopless_phb_optimal.csv")
    n_changed_phb = (phb_df["delta"].abs() > DELTA_THRESHOLD).sum()
    print(f"  loopless PHB flux ({PHB_ID}) = {ll_phb.fluxes[PHB_ID]:.6f}")
    print(f"  reactions with non-trivial loopless correction: {n_changed_phb}")
except Exception as e:
    print(f"  WARNING: loopless_solution failed at PHB optimum: {e}")
    phb_df = None

# ---- Compare against the original (non-loopless) pFBA result from step_7.3 ----
try:
    orig = pd.read_csv(NON_LOOPLESS_PFBA_CSV).set_index("reaction_id")
    if phb_df is not None:
        cmp = phb_df[["loopless_flux"]].join(orig[["flux"]].rename(columns={"flux": "step_7.3_flux"}))
        cmp["delta_vs_step_7.3"] = cmp["loopless_flux"] - cmp["step_7.3_flux"]
        cmp = cmp.reindex(cmp["delta_vs_step_7.3"].abs().sort_values(ascending=False).index)
        cmp.to_csv("./step_7.25_loopless_vs_standard_delta.csv")
        print("\nTop 15 reactions most affected by loopless correction (vs. step_7.3 pFBA):")
        print(cmp.head(15).to_string())
except FileNotFoundError:
    print(f"\nNOTE: {NON_LOOPLESS_PFBA_CSV} not found — skipping comparison to step_7.3.")
    cmp = None

with open("./step_7.25_summary.txt", "w") as f:
    f.write("STEP 7.25 — Loopless re-verification summary\n")
    f.write("=" * 60 + "\n")
    f.write(f"Growth-optimal objective (standard): {sol_growth.objective_value:.6f}\n")
    f.write(f"PHB-optimal objective (standard):    {sol_phb.objective_value:.6f}\n")
    if growth_df is not None:
        f.write(f"Reactions changed by loopless correction (growth optimum): "
                f"{(growth_df['delta'].abs() > DELTA_THRESHOLD).sum()}\n")
    if phb_df is not None:
        f.write(f"Reactions changed by loopless correction (PHB optimum): "
                f"{(phb_df['delta'].abs() > DELTA_THRESHOLD).sum()}\n")
    if cmp is not None:
        f.write("\nInterpretation: large |delta_vs_step_7.3| values indicate reactions\n"
                "whose step_7.3 flux was likely inflated by a stoichiometrically\n"
                "balanced cycle (SBC) rather than reflecting a real metabolic route.\n"
                "Treat step_7.3/7.13/7.16/7.23 flux magnitudes for these reactions\n"
                "with caution until reconciled against this loopless result.\n")

print("\nSummary written to: ./step_7.25_summary.txt")
