"""
STEP 7.32 — BOTTLENECK IDENTIFICATION: WHY MALTOSE DIFFERS FROM EVERY OTHER
CARBON SOURCE

step_7.29 (fixed) confirmed that glucose, fructose, glycerol, sucrose,
mannose, mannitol, arabinose, melibiose, acetate, succinate, and pyruvate
all cap out at an identical growth rate (0.3722 /h aerobic) that matches
the O2-uptake ceiling from step_7.8 -- while maltose alone reaches the
model's true unconstrained maximum (0.4454 /h). That pattern means most
substrates are hitting a SHARED downstream bottleneck that has nothing to
do with carbon supply, and maltose's pathway apparently bypasses it.

This script finds that bottleneck by:
  1. Running pFBA for every candidate carbon source (aerobic, single
     substrate open, matching step_7.29's convention).
  2. Flagging every reaction whose flux magnitude is saturated at its
     bound (the classic definition of a binding/limiting reaction).
  3. Reporting which reactions are saturated in EVERY non-maltose
     substrate run but NOT saturated (or not even active) in the maltose
     run -- these are the strongest bottleneck candidates.
  4. Directly comparing the bounds and pFBA flux on the substrate-specific
     transport/catabolic reactions themselves (e.g., MALT/MALTabc/MALTpts
     for maltose vs. GLCpts/GLCabc/HEX1 for glucose) in case the
     difference is simply a higher capacity bound on the maltose pathway.

Input:  ../models/model_7.xml
Output: ./step_7.32_saturated_reactions_by_substrate.csv
        ./step_7.32_shared_bottleneck_candidates.csv
        ./step_7.32_substrate_pathway_bounds.csv
        ./step_7.32_summary.txt
"""

import cobra
import pandas as pd
from cobra.flux_analysis import pfba

MODEL_PATH = "../models/model_7.xml"
BIOMASS_ID = "Growth"
GLC_UPTAKE_BOUND = 10.0
SATURATION_TOL = 1e-6

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

# First reaction(s) each substrate passes through after uptake, for a direct
# bounds/flux comparison. IDs are checked against the model before use.
SUBSTRATE_PATHWAY_REACTIONS = {
    "Glucose": ["EX_glc__D_e", "GLCpts", "GLCabc", "GLCptspp", "HEX1", "PGI"],
    "Maltose": ["EX_malt_e", "MALTabc", "MALTpts", "MALTt2", "MALT"],
    "Fructose": ["EX_fru_e", "FRUpts", "FRUpts2", "HEX7"],
    "Sucrose": ["EX_sucr_e", "SUCpts", "SUCR", "SUCRt2"],
}

print(f"Loading model from: {MODEL_PATH}")
model = cobra.io.read_sbml_model(MODEL_PATH)
print(f"Loaded model '{model.id}': {len(model.reactions)} reactions, "
      f"{len(model.metabolites)} metabolites, {len(model.genes)} genes")

present = {name: ex_id for name, ex_id in CANDIDATE_CARBON_EXCHANGES.items()
           if ex_id in model.reactions}
missing = set(CANDIDATE_CARBON_EXCHANGES) - set(present)
if missing:
    print(f"  (skipping, not found in model: {sorted(missing)})")

o2_rxn = model.reactions.get_by_id("EX_o2_e")
original_bounds = {ex_id: model.reactions.get_by_id(ex_id).bounds for ex_id in present.values()}
original_o2_bounds = o2_rxn.bounds

print("\n" + "=" * 70)
print("STEP 7.32 — SUBSTRATE-SPECIFIC BOTTLENECK IDENTIFICATION")
print("=" * 70)

all_saturated = {}   # substrate -> set of saturated reaction ids
growth_by_substrate = {}
flux_by_substrate = {}

for carbon_name, carbon_ex_id in present.items():
    for ex_id in present.values():
        model.reactions.get_by_id(ex_id).lower_bound = 0.0
    model.reactions.get_by_id(carbon_ex_id).lower_bound = -GLC_UPTAKE_BOUND
    o2_rxn.lower_bound = -20.0  # "Aerobic" level, matches step_7.29

    model.objective = BIOMASS_ID
    growth = model.slim_optimize()
    growth_by_substrate[carbon_name] = growth

    if growth is None or growth != growth or growth < 1e-6:
        all_saturated[carbon_name] = set()
        flux_by_substrate[carbon_name] = {}
        print(f"  {carbon_name:<12} infeasible/zero growth — skipping saturation scan")
        continue

    sol = pfba(model)
    flux_by_substrate[carbon_name] = sol.fluxes.to_dict()

    saturated = set()
    for rxn in model.reactions:
        flux = sol.fluxes[rxn.id]
        lb, ub = rxn.bounds
        if ub > 0 and abs(flux - ub) < SATURATION_TOL and abs(ub) > SATURATION_TOL:
            saturated.add(rxn.id)
        elif lb < 0 and abs(flux - lb) < SATURATION_TOL and abs(lb) > SATURATION_TOL:
            saturated.add(rxn.id)
    all_saturated[carbon_name] = saturated
    print(f"  {carbon_name:<12} growth={growth:.4f} /h, "
          f"{len(saturated)} reactions saturated at their bound")

for ex_id, bounds in original_bounds.items():
    model.reactions.get_by_id(ex_id).bounds = bounds
o2_rxn.bounds = original_o2_bounds

# ---- Long-format saturated-reactions table ----
rows = []
for substrate, rxns in all_saturated.items():
    for rid in rxns:
        rows.append(dict(carbon_source=substrate, reaction_id=rid,
                          flux=flux_by_substrate[substrate].get(rid)))
sat_df = pd.DataFrame(rows)
sat_df.to_csv("./step_7.32_saturated_reactions_by_substrate.csv", index=False)

# ---- Find reactions saturated in EVERY non-maltose substrate but not in maltose ----
non_maltose = [s for s in present if s != "Maltose" and growth_by_substrate.get(s, 0)]
cand_df = pd.DataFrame()
if non_maltose:
    common_non_maltose = set.intersection(*(all_saturated[s] for s in non_maltose))
    maltose_saturated = all_saturated.get("Maltose", set())
    bottleneck_candidates = common_non_maltose - maltose_saturated

    cand_rows = []
    for rid in bottleneck_candidates:
        rxn = model.reactions.get_by_id(rid)
        cand_rows.append(dict(
            reaction_id=rid,
            reaction_name=rxn.name,
            bounds=str(rxn.bounds),
            flux_glucose=flux_by_substrate.get("Glucose", {}).get(rid),
            flux_maltose=flux_by_substrate.get("Maltose", {}).get(rid),
        ))
    cand_df = pd.DataFrame(cand_rows)
    if len(cand_df):
        cand_df = cand_df.sort_values("reaction_id")
    cand_df.to_csv("./step_7.32_shared_bottleneck_candidates.csv", index=False)

    print(f"\nReactions saturated in ALL {len(non_maltose)} non-maltose substrates "
          f"but NOT in maltose: {len(bottleneck_candidates)}")
    if len(cand_df):
        print(cand_df.to_string(index=False))
else:
    print("\nNo non-maltose substrates produced positive growth — cannot compare.")

# ---- Direct bounds/flux comparison for each substrate's own uptake pathway ----
path_rows = []
for substrate, rxn_ids in SUBSTRATE_PATHWAY_REACTIONS.items():
    for rid in rxn_ids:
        if rid in model.reactions:
            rxn = model.reactions.get_by_id(rid)
            path_rows.append(dict(
                substrate=substrate, reaction_id=rid, reaction_name=rxn.name,
                bounds=str(rxn.bounds),
                flux_in_own_scenario=flux_by_substrate.get(substrate, {}).get(rid),
            ))
path_df = pd.DataFrame(path_rows)
path_df.to_csv("./step_7.32_substrate_pathway_bounds.csv", index=False)
print("\nSubstrate-specific pathway bounds/flux:")
print(path_df.to_string(index=False))

with open("./step_7.32_summary.txt", "w") as f:
    f.write("STEP 7.32 — Bottleneck identification summary\n")
    f.write("=" * 60 + "\n")
    f.write("Growth rate by substrate (aerobic, single substrate open):\n")
    for s, g in growth_by_substrate.items():
        f.write(f"  {s:<12} {g:.4f} /h\n")
    f.write(f"\nReactions saturated in every non-maltose substrate but not in "
            f"maltose: {len(cand_df) if len(cand_df) else 0}\n")
    if len(cand_df):
        f.write("Top candidates (see step_7.32_shared_bottleneck_candidates.csv "
                "for full list and bounds):\n")
        f.write(cand_df.head(15).to_string(index=False))
        f.write("\n\nInterpretation: these reactions cap flux at an identical value "
                "regardless of how much carbon or oxygen is supplied for any "
                "substrate except maltose. Check whether maltose's uptake/"
                "catabolic route (MALTabc/MALTpts/MALT) bypasses one of these "
                "reactions entirely, or whether it simply has a larger bound "
                "capacity -- see step_7.32_substrate_pathway_bounds.csv.\n")
    else:
        f.write("No common bottleneck reaction found via this method -- "
                "consider widening the saturation tolerance or checking FVA "
                "ranges instead of a single pFBA solution.\n")

print("\nSummary written to: ./step_7.32_summary.txt")
