"""
step_7.13.py

PHB pathway flexibility bar chart: runs Flux Variability Analysis (FVA)
on every reaction within FLEXIBILITY_MAX_DISTANCE hops of PHB synthesis
(reusing step_7.11's pathway-distance ranking) and plots each reaction's
feasible flux range (max - min), a measure of how much "slack" or
flexibility the network has at that step. Narrow ranges indicate rigid,
tightly-constrained reactions (little room for flux to reroute); wide
ranges indicate flexible steps with alternative routing options.

Requires step_7.11 to have been run first.

Outputs:
    step_7_outputs/step_7.13_pathway_flexibility.csv
    step_7_outputs/step_7.13_pathway_flexibility_bar_chart.png
"""

import os
import pandas as pd
import cobra
from cobra.flux_analysis import flux_variability_analysis
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from step_7_config import MODEL_PATH, OUTPUT_DIR

FLEXIBILITY_MAX_DISTANCE = 6
MAX_REACTIONS_TO_PLOT = 150
FVA_FRACTION_OF_OPTIMUM = 0.9  # matches the growth-coupled fraction used elsewhere

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

    dist_df = load_pathway_distances()
    scoped = dist_df[dist_df["distance_from_PHB"] <= FLEXIBILITY_MAX_DISTANCE].copy()
    scoped = scoped.sort_values("distance_from_PHB")

    # Cap item list length inline to keep plot layout clean
    raw_ids = scoped["reaction_id"].tolist()
    if len(raw_ids) > MAX_REACTIONS_TO_PLOT:
        print(f"  (Capping plot from {len(raw_ids)} to {MAX_REACTIONS_TO_PLOT} reactions)")
        scoped_ids = raw_ids[:MAX_REACTIONS_TO_PLOT]
    else:
        scoped_ids = raw_ids

    scoped = scoped[scoped["reaction_id"].isin(scoped_ids)]
    print(f"Running FVA on {len(scoped_ids)} reactions within {FLEXIBILITY_MAX_DISTANCE} "
          f"hops of PHB synthesis (fraction_of_optimum={FVA_FRACTION_OF_OPTIMUM})...")

    fva = flux_variability_analysis(
        model, reaction_list=scoped_ids, fraction_of_optimum=FVA_FRACTION_OF_OPTIMUM
    )
    fva["reaction_id"] = fva.index
    fva["flux_range"] = fva["maximum"] - fva["minimum"]

    merged = scoped.merge(fva[["reaction_id", "minimum", "maximum", "flux_range"]], on="reaction_id")
    merged = merged.sort_values(["distance_from_PHB", "reaction_id"])

    csv_path = os.path.join(OUTPUT_DIR, "step_7.13_pathway_flexibility.csv")
    merged.to_csv(csv_path, index=False)
    print(f"Flexibility table written to: {csv_path}")

    # --- Bar chart, ordered by pathway distance, colored by distance --------
    # Compute capped dynamic figure size based on item count
    fig_width = max(10.0, min(24.0, len(merged) * 0.2))
    fig, ax = plt.subplots(figsize=(fig_width, 6))
    max_dist = merged["distance_from_PHB"].max()
    cmap = cm.get_cmap("viridis_r")
    colors = [cmap(d / max_dist) if max_dist > 0 else cmap(0.0) for d in merged["distance_from_PHB"]]

    ax.bar(range(len(merged)), merged["flux_range"], color=colors)
    ax.set_xticks(range(len(merged)))
    label_fontsize = 6 if len(merged) <= 100 else 4
    ax.set_xticklabels(merged["reaction_id"], rotation=90, fontsize=label_fontsize)
    ax.set_ylabel("Feasible flux range (max - min, mmol/gDW/h)")
    ax.set_title(
        f"PHB Pathway Flexibility (reactions within {FLEXIBILITY_MAX_DISTANCE} hops of PHB synthesis)\n"
        f"color = distance from PHB synthesis (darker = closer)"
    )

    sm = cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=max_dist))
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="Distance from PHB synthesis")
    fig.tight_layout()

    png_path = os.path.join(OUTPUT_DIR, "step_7.13_pathway_flexibility_bar_chart.png")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    print(f"Bar chart written to: {png_path}")

    rigid = merged.sort_values("flux_range").head(10)
    flexible = merged.sort_values("flux_range", ascending=False).head(10)
    print("\nMost rigid (least flexible) reactions:")
    print(rigid[["reaction_id", "distance_from_PHB", "flux_range"]].to_string(index=False))
    print("\nMost flexible reactions:")
    print(flexible[["reaction_id", "distance_from_PHB", "flux_range"]].to_string(index=False))


if __name__ == "__main__":
    main()
