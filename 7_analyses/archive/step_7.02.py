"""
step_7.2_fba_growth_and_phb.py

Flux Balance Analysis (FBA) for:
    1. Maximum growth rate (biomass objective), and the resulting PHB flux
    2. Maximum PHB production (PHB objective), and the resulting growth rate
    3. Growth-coupled PHB production: growth constrained to a fraction of
       its maximum (default 90%), then PHB maximized subject to that —
       a more "physiologically realistic" production scenario than
       unconstrained PHB maximization, since real cells rarely divert
       100% of resources away from growth.

Outputs:
    step_7_outputs/step_7.2_fba_results.csv
    step_7_outputs/step_7.2_fba_results.txt
"""

import csv
import os

from step_7_utils import (
    MODEL_PATH,
    OUTPUT_DIR,
    describe_reaction,
    ensure_output_dir,
    find_biomass_reaction,
    find_phb_reaction,
    load_model,
)

# --------------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------------
# If auto-detection picks the wrong PHB reaction, set its ID here directly,
# e.g. PHB_REACTION_ID = "EX_phb_e"
PHB_REACTION_ID = "PHBS_syn_1"

# Fraction of maximum growth to hold fixed for the growth-coupled scenario
GROWTH_COUPLED_FRACTION = 0.90


def get_reactions(model):
    biomass = find_biomass_reaction(model)
    phb = model.reactions.get_by_id(PHB_REACTION_ID) if PHB_REACTION_ID else find_phb_reaction(model)
    if phb is None:
        raise ValueError(
            "Could not auto-detect a PHB reaction. Set PHB_REACTION_ID at the "
            "top of this script to the correct reaction ID (e.g. 'EX_phb_e')."
        )
    return biomass, phb


def main():
    ensure_output_dir(OUTPUT_DIR)
    model = load_model(MODEL_PATH)
    biomass, phb = get_reactions(model)

    describe_reaction(biomass, "Biomass reaction")
    describe_reaction(phb, "PHB reaction")

    results = []

    # --- Scenario 1: maximize growth, report resulting PHB flux -----------
    with model as m:
        m.objective = biomass.id
        sol = m.optimize()
        growth_max = sol.objective_value
        phb_at_growth_max = sol.fluxes[phb.id]
    results.append({
        "scenario": "max_growth",
        "status": sol.status,
        "growth_rate_1_h": growth_max,
        "phb_flux": phb_at_growth_max,
    })

    # --- Scenario 2: maximize PHB, report resulting growth rate -----------
    with model as m:
        m.objective = phb.id
        sol = m.optimize()
        phb_max = sol.objective_value
        growth_at_phb_max = sol.fluxes[biomass.id]
    results.append({
        "scenario": "max_PHB",
        "status": sol.status,
        "growth_rate_1_h": growth_at_phb_max,
        "phb_flux": phb_max,
    })

    # --- Scenario 3: growth-coupled PHB production -------------------------
    with model as m:
        m.objective = biomass.id
        max_growth = m.slim_optimize()
        fixed_growth = GROWTH_COUPLED_FRACTION * max_growth
        biomass_rxn_m = m.reactions.get_by_id(biomass.id)
        biomass_rxn_m.lower_bound = fixed_growth
        biomass_rxn_m.upper_bound = fixed_growth
        m.objective = phb.id
        sol = m.optimize()
        phb_at_fixed_growth = sol.objective_value
    results.append({
        "scenario": f"growth_coupled_{int(GROWTH_COUPLED_FRACTION * 100)}pct",
        "status": sol.status,
        "growth_rate_1_h": fixed_growth,
        "phb_flux": phb_at_fixed_growth,
    })

    # --- Write outputs -------------------------------------------------------
    csv_path = os.path.join(OUTPUT_DIR, "step_7.2_fba_results.csv")
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["scenario", "status", "growth_rate_1_h", "phb_flux"])
        writer.writeheader()
        writer.writerows(results)

    txt_lines = ["=" * 70, "STEP 7.2 — FBA: GROWTH AND PHB PRODUCTION", "=" * 70]
    txt_lines.append(f"Biomass reaction: {biomass.id}")
    txt_lines.append(f"PHB reaction:     {phb.id}")
    txt_lines.append("")
    for r in results:
        txt_lines.append(
            f"[{r['scenario']}] status={r['status']} "
            f"growth={r['growth_rate_1_h']:.6f} /h  PHB_flux={r['phb_flux']:.6f}"
        )
    txt_path = os.path.join(OUTPUT_DIR, "step_7.2_fba_results.txt")
    with open(txt_path, "w") as fh:
        fh.write("\n".join(txt_lines))

    print("\n".join(txt_lines))
    print(f"\nResults written to:\n  {csv_path}\n  {txt_path}")


if __name__ == "__main__":
    main()
