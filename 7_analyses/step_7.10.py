"""
step_7.10.py

Shadow price and reduced cost analysis, used to identify which nutrients
and reactions are actually limiting the optimal solution.

Definitions (standard LP duality, as exposed by cobrapy's Solution object):
    - shadow_prices: the dual value of each metabolite's mass-balance
      constraint. A strongly negative shadow price on an exchanged
      metabolite (e.g. glucose, O2, NH4) means the objective would
      improve substantially if more of that metabolite were made
      available — i.e. it is a *limiting nutrient*.
    - reduced_costs: the dual value associated with each reaction's flux
      bounds. A nonzero reduced cost on a reaction sitting at one of its
      bounds means that bound is *binding* / actively constraining the
      optimum; relaxing it would change the objective.

This script computes both at two optima for comparison:
    (a) the growth optimum (objective = biomass)
    (b) the PHB optimum (objective = PHB, growth unconstrained)

Outputs:
    step_7_outputs/step_7.10_shadow_prices_growth_optimum.csv
    step_7_outputs/step_7.10_shadow_prices_PHB_optimum.csv
    step_7_outputs/step_7.10_reduced_costs_growth_optimum.csv
    step_7_outputs/step_7.10_reduced_costs_PHB_optimum.csv
    step_7_outputs/step_7.10_limiting_nutrients_summary.txt
    step_7_outputs/step_7.10_shadow_price_reduced_cost.png
"""

import os
import pandas as pd
import cobra
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from step_7_config import (
    MODEL_PATH,
    OUTPUT_DIR,
    BIOMASS_ID,
    PHB_REACTION_ID,
    describe_reaction,
)

TOP_N = 15  # how many top-ranked metabolites/reactions to show per panel


def solve_and_get_duals(model, objective_id):
    with model as m:
        m.objective = objective_id
        sol = m.optimize()
        if sol.status != "optimal":
            raise RuntimeError(f"Solver status '{sol.status}' when optimizing {objective_id}")
        # copy out of the context before it's reverted
        shadow_prices = sol.shadow_prices.copy()
        reduced_costs = sol.reduced_costs.copy()
        objective_value = sol.objective_value
    return objective_value, shadow_prices, reduced_costs


def exchange_metabolite_table(model, shadow_prices):
    """Restrict shadow prices to metabolites involved in exchange
    reactions (the ones a person could actually feed/withhold), since
    these are the directly actionable 'limiting nutrients'."""
    exch_mets = set()
    for rxn in model.exchanges:
        exch_mets.update(m.id for m in rxn.metabolites)

    rows = []
    for met_id, price in shadow_prices.items():
        if met_id in exch_mets:
            met = model.metabolites.get_by_id(met_id)
            rows.append({"metabolite_id": met_id, "metabolite_name": met.name, "shadow_price": price})
    df = pd.DataFrame(rows).sort_values("shadow_price")
    return df


def exchange_reaction_table(model, reduced_costs):
    rows = []
    for rxn in model.exchanges:
        rc = reduced_costs.get(rxn.id, 0.0)
        rows.append({
            "reaction_id": rxn.id, "reaction_name": rxn.name,
            "reduced_cost": rc, "lower_bound": rxn.lower_bound, "upper_bound": rxn.upper_bound,
        })
    df = pd.DataFrame(rows)
    df["abs_reduced_cost"] = df["reduced_cost"].abs()
    return df.sort_values("abs_reduced_cost", ascending=False).drop(columns="abs_reduced_cost")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    model = cobra.io.read_sbml_model(MODEL_PATH)

    biomass = model.reactions.get_by_id(BIOMASS_ID)
    phb = model.reactions.get_by_id(PHB_REACTION_ID)

    describe_reaction(biomass, "Biomass reaction")
    describe_reaction(phb, "PHB reaction")

    results = {}
    for label, obj_id in [("growth_optimum", biomass.id), ("PHB_optimum", phb.id)]:
        print(f"\nSolving at {label} (objective = {obj_id})...")
        obj_val, shadow_prices, reduced_costs = solve_and_get_duals(model, obj_id)
        print(f"  objective value = {obj_val:.6f}")

        exch_met_df = exchange_metabolite_table(model, shadow_prices)
        exch_rxn_df = exchange_reaction_table(model, reduced_costs)

        sp_csv = os.path.join(OUTPUT_DIR, f"step_7.10_shadow_prices_{label}.csv")
        rc_csv = os.path.join(OUTPUT_DIR, f"step_7.10_reduced_costs_{label}.csv")
        exch_met_df.to_csv(sp_csv, index=False)
        exch_rxn_df.to_csv(rc_csv, index=False)
        print(f"  shadow prices  -> {sp_csv}")
        print(f"  reduced costs  -> {rc_csv}")

        results[label] = {
            "objective_value": obj_val,
            "shadow_prices": exch_met_df,
            "reduced_costs": exch_rxn_df,
        }

    # --- Text summary of the most limiting nutrients/reactions -------------
    summary_lines = ["STEP 7.10 — Shadow Price & Reduced Cost Summary", "=" * 70]
    for label in results:
        sp = results[label]["shadow_prices"]
        rc = results[label]["reduced_costs"]
        summary_lines.append(f"\n[{label}] objective value = {results[label]['objective_value']:.6f}")

        most_limiting = sp[sp["shadow_price"] < -1e-9].head(TOP_N)
        if len(most_limiting):
            summary_lines.append("  Most limiting exchanged metabolites (most negative shadow price):")
            for _, row in most_limiting.iterrows():
                summary_lines.append(
                    f"    {row['metabolite_id']:15s} {row['shadow_price']:12.4f}  {row['metabolite_name']}"
                )
        else:
            summary_lines.append("  No metabolites with a meaningfully negative shadow price found.")

        binding = rc[rc["reduced_cost"].abs() > 1e-9].head(TOP_N)
        if len(binding):
            summary_lines.append("  Reactions with a binding (nonzero reduced-cost) bound:")
            for _, row in binding.iterrows():
                summary_lines.append(
                    f"    {row['reaction_id']:15s} {row['reduced_cost']:12.4f}  {row['reaction_name']}"
                )
        else:
            summary_lines.append("  No reactions with a nonzero reduced cost found.")

    summary_path = os.path.join(OUTPUT_DIR, "step_7.10_limiting_nutrients_summary.txt")
    with open(summary_path, "w") as fh:
        fh.write("\n".join(summary_lines))
    print("\n" + "\n".join(summary_lines))
    print(f"\nSummary written to: {summary_path}")

    # --- Figure: top shadow prices and reduced costs, both optima ----------
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for col, label in enumerate(results):
        sp = results[label]["shadow_prices"]
        rc = results[label]["reduced_costs"]

        sp_top = pd.concat([sp.head(TOP_N // 2), sp.tail(TOP_N // 2)]).drop_duplicates()
        ax_sp = axes[0, col]
        colors = ["tab:red" if v < 0 else "tab:blue" for v in sp_top["shadow_price"]]
        ax_sp.barh(sp_top["metabolite_id"], sp_top["shadow_price"], color=colors)
        ax_sp.set_title(f"Shadow prices — {label}")
        ax_sp.axvline(0, color="black", lw=0.8)

        rc_top = rc.head(TOP_N)
        ax_rc = axes[1, col]
        ax_rc.barh(rc_top["reaction_id"], rc_top["reduced_cost"], color="tab:purple")
        ax_rc.set_title(f"Top |reduced cost| exchanges — {label}")
        ax_rc.axvline(0, color="black", lw=0.8)

    fig.suptitle("Shadow Price & Reduced Cost Analysis", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    png_path = os.path.join(OUTPUT_DIR, "step_7.10_shadow_price_reduced_cost.png")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    print(f"Figure written to: {png_path}")


if __name__ == "__main__":
    main()
