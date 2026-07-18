"""
Flux Balance Analysis (FBA) for:
    1. Maximum growth rate (biomass objective), and the resulting PHB flux
    2. Maximum PHB production (PHB objective), and the resulting growth rate
    3. Growth-coupled PHB production: growth constrained to a fraction of its maximum (default 90%), then PHB maximized subject to that a more "physiologically realistic" production scenario than unconstrained PHB maximization, since real cells rarely divert 100% of resources away from growth.

Outputs:
    step_7.01_fba_results.csv
    step_7.01_fba_results.txt
"""

import csv
import os
import cobra
from step_7_config import (
    MODEL_PATH,
    OUTPUT_DIR,
    BIOMASS_ID,
    PHB_REACTION_ID,
    growth_coupled_phb
)

# fraction of maximum growth to hold fixed for the growth-coupled scenario
GROWTH_COUPLED_FRACTION = 0.90

def main():
    model = cobra.io.read_sbml_model(MODEL_PATH)
    biomass_reaction = model.reactions.get_by_id(BIOMASS_ID)
    phb_reaction = model.reactions.get_by_id(PHB_REACTION_ID)

    results = []

    # --- Scenario 1: maximize growth, report resulting PHB flux -----------
    with model as m:
        m.objective = BIOMASS_ID
        sol = m.optimize()
        growth_max = sol.objective_value
        phb_at_growth_max = sol.fluxes[PHB_REACTION_ID]
    results.append({
        "scenario": "max_growth",
        "status": sol.status,
        "growth_rate_1_h": growth_max,
        "phb_flux": phb_at_growth_max,
    })

    # --- Scenario 2: maximize PHB, report resulting growth rate -----------
    with model as m:
        m.objective = PHB_REACTION_ID
        sol = m.optimize()
        phb_max = sol.objective_value
        growth_at_phb_max = sol.fluxes[BIOMASS_ID]
    results.append({
        "scenario": "max_PHB",
        "status": sol.status,
        "growth_rate_1_h": growth_at_phb_max,
        "phb_flux": phb_max,
    })

    # --- Scenario 3: growth-coupled PHB production -------------------------
    # Utilizing the optimized helper function from step_7_config
    fixed_growth, phb_at_fixed_growth = growth_coupled_phb(
        model,
        biomass=biomass_reaction,
        phb=phb_reaction,
        fraction=GROWTH_COUPLED_FRACTION
    )

    # Run a quick optimize inside context just to capture the solver status for logs
    with model as m:
        m.objective = BIOMASS_ID
        sol = m.optimize()

    results.append({
        "scenario": f"growth_coupled_{int(GROWTH_COUPLED_FRACTION * 100)}pct",
        "status": sol.status,
        "growth_rate_1_h": fixed_growth,
        "phb_flux": phb_at_fixed_growth,
    })

    # --- Write outputs -------------------------------------------------------
    csv_path = os.path.join(OUTPUT_DIR, "step_7.01_fba_results.csv")
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["scenario", "status", "growth_rate_1_h", "phb_flux"])
        writer.writeheader()
        writer.writerows(results)

    txt_lines = ["=" * 70, "STEP 7.01 — FBA: GROWTH AND PHB PRODUCTION", "=" * 70]
    txt_lines.append(f"Biomass reaction: {BIOMASS_ID}")
    txt_lines.append(f"PHB reaction:     {PHB_REACTION_ID}")
    txt_lines.append("")
    for r in results:
        txt_lines.append(
            f"[{r['scenario']}] status={r['status']} "
            f"growth={r['growth_rate_1_h']:.6f} /h  PHB_flux={r['phb_flux']:.6f}"
        )
    txt_path = os.path.join(OUTPUT_DIR, "step_7.01_fba_results.txt")
    with open(txt_path, "w") as fh:
        fh.write("\n".join(txt_lines))

    print("\n".join(txt_lines))
    print(f"\nResults written to:\n  {csv_path}\n  {txt_path}")


if __name__ == "__main__":
    main()
