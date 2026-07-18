"""
step_7.12.py

Pathway bottleneck analysis: identifies "high-flux burden points" in and
around the PHB biosynthetic pathway — reactions that are either carrying
unusually large flux relative to their neighborhood, or are operating
very close to their maximum allowed flux (bound-constrained), at the PHB
production optimum. Both are signs of a potential bottleneck: a reaction
that would need to be up-expressed / engineered / relieved of its bound
to increase PHB flux further.

Requires step_7.11 to have been run first (uses its pathway-distance
ranking to scope the analysis to the PHB neighborhood).

Method:
    1. Run pFBA at the PHB production optimum (growth-coupled at 90% of
       max growth, consistent with earlier scripts).
    2. Restrict to reactions within BOTTLENECK_MAX_DISTANCE of PHB
       synthesis (from step_7.11's output).
    3. For each such reaction compute:
         - |flux| (the raw "burden")
         - flux/bound ratio: |flux| / max(|lower_bound|, |upper_bound|)
           (how close to its capacity ceiling the reaction is running;
           values near 1.0 flag a binding capacity bottleneck)
    4. Rank by a combined bottleneck score and report/plot the top hits.

Outputs:
    step_7_outputs/step_7.12_bottleneck_analysis.csv
    step_7_outputs/step_7.12_bottleneck_bar_chart.png
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

BOTTLENECK_MAX_DISTANCE = 8   # only consider reactions within this many hops of PHB synthesis
GROWTH_COUPLED_FRACTION = 0.90
TOP_N = 20

DISTANCE_CSV = os.path.join(OUTPUT_DIR, "step_7.11_phb_pathway_distances.csv")


def load_pathway_distances():
    if not os.path.exists(DISTANCE_CSV):
        raise SystemExit(
            f"Could not find {DISTANCE_CSV}. Run step_7.11.py "
            f"first — this script reuses its pathway-distance ranking."
        )
    return pd.read_csv(DISTANCE_CSV)


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    model = cobra.io.read_sbml_model(MODEL_PATH)

    biomass = model.reactions.get_by_id(BIOMASS_ID)
    phb = model.reactions.get_by_id(PHB_REACTION_ID)

    describe_reaction(biomass, "Biomass reaction")
    describe_reaction(phb, "PHB reaction")

    dist_df = load_pathway_distances()
    scoped_ids = set(
        dist_df.loc[dist_df["distance_from_PHB"] <= BOTTLENECK_MAX_DISTANCE, "reaction_id"]
    )
    print(f"Scoping bottleneck analysis to {len(scoped_ids)} reactions within "
          f"{BOTTLENECK_MAX_DISTANCE} hops of PHB synthesis.")

    # --- Solve at the PHB production optimum (growth-coupled), pFBA ---------
    with model:
        model.objective = biomass.id
        max_growth = model.slim_optimize()
        fixed_growth = GROWTH_COUPLED_FRACTION * max_growth
        biomass_rxn = model.reactions.get_by_id(biomass.id)
        biomass_rxn.lower_bound = fixed_growth
        biomass_rxn.upper_bound = fixed_growth
        model.objective = phb.id
        solution = pfba(model)

    rows = []
    dist_lookup = dict(zip(dist_df["reaction_id"], dist_df["distance_from_PHB"]))
    for rxn_id in scoped_ids:
        rxn = model.reactions.get_by_id(rxn_id)
        flux = solution.fluxes[rxn_id]
        capacity = max(abs(rxn.lower_bound), abs(rxn.upper_bound))
        ratio = abs(flux) / capacity if capacity > 1e-9 else 0.0
        rows.append({
            "reaction_id": rxn_id,
            "reaction_name": rxn.name,
            "distance_from_PHB": dist_lookup.get(rxn_id),
            "flux": flux,
            "abs_flux": abs(flux),
            "capacity_bound": capacity,
            "flux_to_capacity_ratio": ratio,
        })

    df = pd.DataFrame(rows)
    # combined bottleneck score: normalized flux magnitude x flux/capacity ratio
    max_abs_flux = df["abs_flux"].max() if df["abs_flux"].max() > 0 else 1.0
    df["bottleneck_score"] = (df["abs_flux"] / max_abs_flux) * df["flux_to_capacity_ratio"]
    df = df.sort_values("bottleneck_score", ascending=False)

    csv_path = os.path.join(OUTPUT_DIR, "step_7.12_bottleneck_analysis.csv")
    df.to_csv(csv_path, index=False)
    print(f"\nFull bottleneck table written to: {csv_path}")

    top = df.head(TOP_N)
    print(f"\nTop {TOP_N} bottleneck candidates:")
    print(top[["reaction_id", "distance_from_PHB", "abs_flux", "flux_to_capacity_ratio", "bottleneck_score"]]
          .to_string(index=False))

    # --- Bar chart ---------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9, max(5, 0.35 * len(top))))
    y_pos = range(len(top))
    bars = ax.barh(y_pos, top["bottleneck_score"], color="tab:red")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top["reaction_id"])
    ax.invert_yaxis()
    ax.set_xlabel("Bottleneck score (normalized |flux| x flux/capacity ratio)")
    ax.set_title(f"Top {TOP_N} Pathway Bottleneck Candidates at PHB Optimum")
    for bar, (_, row) in zip(bars, top.iterrows()):
        ax.text(bar.get_width(), bar.get_y() + bar.get_height() / 2,
                f"  d={int(row['distance_from_PHB'])}, ratio={row['flux_to_capacity_ratio']:.2f}",
                va="center", fontsize=7)
    fig.tight_layout()

    png_path = os.path.join(OUTPUT_DIR, "step_7.12_bottleneck_bar_chart.png")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    print(f"\nBar chart written to: {png_path}")


if __name__ == "__main__":
    main()
