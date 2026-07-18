"""
step_7.14.py

Compares the top-flux exchange reactions and top-flux internal reactions
between two pFBA-derived states:
    - growth-optimal (objective = biomass, unconstrained PHB)
    - PHB-optimal     (objective = PHB, growth-coupled at 90% of its own max)

This highlights which exchanges/pathways the cell relies on most heavily
in each state, and how resource allocation shifts when the objective
moves from growth to PHB production.

Outputs:
    step_7_outputs/step_7.14_flux_comparison_full.csv        (all reactions, both states)
    step_7_outputs/step_7.14_top_exchange_flux_comparison.png
    step_7_outputs/step_7.14_top_internal_flux_comparison.png
"""

import os
import pandas as pd
import cobra
from cobra.flux_analysis import pfba
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

GROWTH_COUPLED_FRACTION = 0.90
TOP_N = 15


def get_pfba_solution(model, objective_id, growth_coupled=False, biomass_id=None):
    with model as m:
        if growth_coupled:
            m.objective = biomass_id
            max_growth = m.slim_optimize()
            fixed_growth = GROWTH_COUPLED_FRACTION * max_growth
            biomass_rxn = m.reactions.get_by_id(biomass_id)
            biomass_rxn.lower_bound = fixed_growth
            biomass_rxn.upper_bound = fixed_growth
        m.objective = objective_id
        sol = pfba(m)
        fluxes = sol.fluxes.copy()
    return fluxes


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    model = cobra.io.read_sbml_model(MODEL_PATH)

    biomass = model.reactions.get_by_id(BIOMASS_ID)
    phb = model.reactions.get_by_id(PHB_REACTION_ID)

    describe_reaction(biomass, "Biomass reaction")
    describe_reaction(phb, "PHB reaction")

    print("\nSolving growth-optimal state (pFBA)...")
    growth_fluxes = get_pfba_solution(model, biomass.id, growth_coupled=False)

    print("Solving PHB-optimal state (growth-coupled, pFBA)...")
    phb_fluxes = get_pfba_solution(model, phb.id, growth_coupled=True, biomass_id=biomass.id)

    exchange_ids = {r.id for r in model.exchanges}
    df = pd.DataFrame({"growth_optimal_flux": growth_fluxes, "PHB_optimal_flux": phb_fluxes})
    df["reaction_id"] = df.index
    df["reaction_name"] = [model.reactions.get_by_id(rid).name for rid in df.index]
    df["is_exchange"] = df["reaction_id"].isin(exchange_ids)
    df["flux_shift"] = df["PHB_optimal_flux"] - df["growth_optimal_flux"]

    csv_path = os.path.join(OUTPUT_DIR, "step_7.14_flux_comparison_full.csv")
    df.to_csv(csv_path, index=False)
    print(f"\nFull flux comparison table written to: {csv_path}")

    def top_flux_subset(sub_df, n):
        sub_df = sub_df.copy()
        sub_df["max_abs_flux"] = sub_df[["growth_optimal_flux", "PHB_optimal_flux"]].abs().max(axis=1)
        return sub_df.sort_values("max_abs_flux", ascending=False).head(n)

    exch_top = top_flux_subset(df[df["is_exchange"]], TOP_N)
    internal_top = top_flux_subset(df[~df["is_exchange"]], TOP_N)

    def make_comparison_bar(sub_df, title, png_name):
        fig, ax = plt.subplots(figsize=(10, max(5, 0.4 * len(sub_df))))
        y = range(len(sub_df))
        height = 0.35
        ax.barh([i + height / 2 for i in y], sub_df["growth_optimal_flux"], height,
                label="Growth-optimal", color="tab:blue")
        ax.barh([i - height / 2 for i in y], sub_df["PHB_optimal_flux"], height,
                label="PHB-optimal", color="tab:green")
        ax.set_yticks(list(y))
        ax.set_yticklabels(sub_df["reaction_id"])
        ax.invert_yaxis()
        ax.axvline(0, color="black", lw=0.8)
        ax.set_xlabel("Flux (mmol/gDW/h)")
        ax.set_title(title)
        ax.legend()
        fig.tight_layout()
        path = os.path.join(OUTPUT_DIR, png_name)
        fig.savefig(path, dpi=200)
        plt.close(fig)
        return path

    exch_png = make_comparison_bar(
        exch_top, f"Top {TOP_N} Exchange Reactions — Growth-Optimal vs. PHB-Optimal",
        "step_7.14_top_exchange_flux_comparison.png",
    )
    internal_png = make_comparison_bar(
        internal_top, f"Top {TOP_N} Internal Reactions — Growth-Optimal vs. PHB-Optimal",
        "step_7.14_top_internal_flux_comparison.png",
    )

    print(f"Exchange comparison figure written to: {exch_png}")
    print(f"Internal comparison figure written to: {internal_png}")


if __name__ == "__main__":
    main()
