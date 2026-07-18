"""
step_7.16.py

Single-reaction knockout screen: systematically sets each non-boundary
reaction's flux to zero (one at a time), re-optimizes growth, and
classifies the reaction by how much growth drops. Reactions whose
knockout collapses growth below ESSENTIAL_THRESHOLD (as a fraction of
wild-type growth) are flagged "essential"; those causing a large but
non-lethal drop are flagged "critical".

Also records the growth-coupled PHB flux after each knockout, so you can
see (a) which knockouts are lethal, and (b) among the non-lethal ones,
which shift the PHB phenotype the most — useful groundwork for the
strain-design steps that follow (step_7.17, step_7.18).

Outputs:
    step_7_outputs/step_7.16_knockout_screen_full.csv
    step_7_outputs/step_7.16_essential_critical_bar_chart.png
"""

import os
import pandas as pd
import cobra
from cobra.flux_analysis import single_reaction_deletion
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

ESSENTIAL_THRESHOLD = 0.01     # growth < 1% of wild-type -> "essential"
CRITICAL_THRESHOLD = 0.50      # growth < 50% of wild-type (but not essential) -> "critical"
GROWTH_COUPLED_FRACTION = 0.90
TOP_N_PLOT = 30


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    model = cobra.io.read_sbml_model(MODEL_PATH)

    biomass = model.reactions.get_by_id(BIOMASS_ID)
    phb = model.reactions.get_by_id(PHB_REACTION_ID)

    model.objective = biomass.id
    wt_growth = model.slim_optimize()
    print(f"\nWild-type growth rate: {wt_growth:.6f} /h")

    print("Running single-reaction knockout screen (this can take a while on "
          "a ~2000-reaction model)...")
    deletion_results = single_reaction_deletion(model)
    # cobra returns 'ids' as a frozenset per row; unpack to a single reaction id
    deletion_results = deletion_results.reset_index(drop=True)
    deletion_results["reaction_id"] = deletion_results["ids"].apply(lambda s: next(iter(s)))
    deletion_results["growth_fraction_of_wt"] = deletion_results["growth"] / wt_growth if wt_growth > 0 else 0.0

    def classify(frac):
        if pd.isna(frac) or frac < ESSENTIAL_THRESHOLD:
            return "essential"
        if frac < CRITICAL_THRESHOLD:
            return "critical"
        return "non-essential"

    deletion_results["classification"] = deletion_results["growth_fraction_of_wt"].apply(classify)

    print("\nClassification counts:")
    print(deletion_results["classification"].value_counts().to_string())

    # --- PHB flux after knockout, for non-lethal knockouts ---------------------
    print("\nComputing growth-coupled PHB flux after each non-lethal knockout "
          "(this is the slower part)...")
    phb_after = []
    for _, row in deletion_results.iterrows():
        rxn_id = row["reaction_id"]
        if row["classification"] == "essential" or phb is None:
            phb_after.append(float("nan"))
            continue
        with model as m:
            m.reactions.get_by_id(rxn_id).knock_out()
            m.objective = biomass.id
            g = m.slim_optimize()
            if not g or g <= 1e-9:
                phb_after.append(float("nan"))
                continue
            fixed_growth = GROWTH_COUPLED_FRACTION * g
            biomass_m = m.reactions.get_by_id(biomass.id)
            biomass_m.lower_bound = fixed_growth
            biomass_m.upper_bound = fixed_growth
            m.objective = phb.id
            phb_flux = m.slim_optimize()
            phb_after.append(phb_flux if phb_flux else 0.0)
    deletion_results["growth_coupled_PHB_after_knockout"] = phb_after

    # --- attach reaction names -------------------------------------------------
    deletion_results["reaction_name"] = deletion_results["reaction_id"].apply(
        lambda rid: model.reactions.get_by_id(rid).name
    )

    csv_path = os.path.join(OUTPUT_DIR, "step_7.16_knockout_screen_full.csv")
    out_cols = ["reaction_id", "reaction_name", "growth", "growth_fraction_of_wt",
                "classification", "growth_coupled_PHB_after_knockout"]
    deletion_results[out_cols].to_csv(csv_path, index=False)
    print(f"\nFull knockout screen written to: {csv_path}")

    # --- Bar chart of most essential/critical reactions ------------------------
    ranked = deletion_results[deletion_results["classification"] != "non-essential"].copy()
    ranked = ranked.sort_values("growth_fraction_of_wt").head(TOP_N_PLOT)

    fig, ax = plt.subplots(figsize=(9, max(5, 0.3 * len(ranked))))
    colors = ["tab:red" if c == "essential" else "tab:orange" for c in ranked["classification"]]
    ax.barh(ranked["reaction_id"], ranked["growth_fraction_of_wt"], color=colors)
    ax.set_xlabel("Growth after knockout (fraction of wild-type)")
    ax.set_title(f"Top {len(ranked)} Essential (red) / Critical (orange) Reactions")
    ax.invert_yaxis()
    fig.tight_layout()

    png_path = os.path.join(OUTPUT_DIR, "step_7.16_essential_critical_bar_chart.png")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    print(f"Bar chart written to: {png_path}")


if __name__ == "__main__":
    main()
