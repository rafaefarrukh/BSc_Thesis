"""
step_7.19_flux_span_direction_classification.py

Flux-span distribution & direction classification: runs FVA (at 100% of
the PHB optimum by default) across every non-boundary reaction and
classifies each one by its feasible flux range:

    - blocked          : min ~ 0 and max ~ 0 (cannot carry any flux)
    - fixed             : min == max, nonzero (fully determined flux)
    - forward_only      : min >= 0, max > 0 (can only carry forward flux)
    - backward_only     : max <= 0, min < 0 (can only carry reverse flux)
    - bidirectional      : min < 0 < max (genuinely reversible in this state)

Also plots the distribution of flux "spans" (max - min) across all
reactions, which is a quick way to see how tightly vs. loosely
constrained the network is overall (a long tail of high-span reactions
suggests a lot of unused flexibility/redundancy; a sharply peaked
near-zero distribution suggests a highly determined network).

By default this runs on the WHOLE model (not just the PHB neighborhood),
since direction/span classification is a general network property. This
is the most expensive script in the Step 7 series if REACTION_SCOPE is
left at "all" on a ~2000-reaction model — reduce scope or use loopless=False
(the default) to keep runtime down.

Outputs:
    step_7_outputs/step_7.19_flux_span_classification.csv
    step_7_outputs/step_7.19_flux_span_and_direction.png
"""

import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cobra.flux_analysis import flux_variability_analysis

from step_7_utils import (
    MODEL_PATH,
    OUTPUT_DIR,
    describe_reaction,
    ensure_output_dir,
    find_phb_reaction,
    load_model,
)

PHB_REACTION_ID = "PHBS_syn_1"
FRACTION_OF_OPTIMUM = 1.0
LOOPLESS = False   # set True for loopless FVA (slower, removes futile-cycle flux)
REACTION_SCOPE = "all"   # "all" or "pathway" (reuses step_7.12's output if available)
PATHWAY_MAX_DISTANCE = 6
EPS = 1e-6

DISTANCE_CSV = os.path.join(OUTPUT_DIR, "step_7.12_phb_pathway_distances.csv")


def get_reaction_scope(model):
    if REACTION_SCOPE == "pathway" and os.path.exists(DISTANCE_CSV):
        dist_df = pd.read_csv(DISTANCE_CSV)
        scoped = dist_df[dist_df["distance_from_PHB"] <= PATHWAY_MAX_DISTANCE]
        return scoped["reaction_id"].tolist()
    return [r.id for r in model.reactions if not r.boundary]


def classify(row, eps=EPS):
    lo, hi = row["minimum"], row["maximum"]
    if abs(lo) < eps and abs(hi) < eps:
        return "blocked"
    if hi - lo < eps:
        return "fixed"
    if lo >= -eps and hi > eps:
        return "forward_only"
    if hi <= eps and lo < -eps:
        return "backward_only"
    if lo < -eps < hi:
        return "bidirectional"
    return "other"


def main():
    ensure_output_dir(OUTPUT_DIR)
    model = load_model(MODEL_PATH)
    phb = model.reactions.get_by_id(PHB_REACTION_ID) if PHB_REACTION_ID else find_phb_reaction(model)
    if phb is None:
        raise ValueError("Could not auto-detect a PHB reaction. Set PHB_REACTION_ID at the top of this script.")
    describe_reaction(phb, "PHB reaction (FVA objective)")

    model.objective = phb.id
    reaction_ids = get_reaction_scope(model)
    print(f"Running FVA (loopless={LOOPLESS}, fraction_of_optimum={FRACTION_OF_OPTIMUM}) "
          f"on {len(reaction_ids)} reactions (this is the slow step)...")

    fva = flux_variability_analysis(
        model, reaction_list=reaction_ids, fraction_of_optimum=FRACTION_OF_OPTIMUM, loopless=LOOPLESS
    )
    fva["reaction_id"] = fva.index
    fva["flux_span"] = fva["maximum"] - fva["minimum"]
    fva["direction_class"] = fva.apply(classify, axis=1)

    csv_path = os.path.join(OUTPUT_DIR, "step_7.19_flux_span_classification.csv")
    fva.to_csv(csv_path, index=False)
    print(f"\nClassification table written to: {csv_path}")

    counts = fva["direction_class"].value_counts()
    print("\nDirection classification counts:")
    print(counts.to_string())

    # --- Figure: span histogram + classification bar chart --------------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    nonzero_spans = fva.loc[fva["flux_span"] > EPS, "flux_span"]
    axes[0].hist(nonzero_spans, bins=40, color="tab:blue", alpha=0.85)
    axes[0].set_xlabel("Flux span (max - min, mmol/gDW/h)")
    axes[0].set_ylabel("Reaction count")
    axes[0].set_title(
        f"Flux-Span Distribution\n(excludes {int((fva['flux_span'] <= EPS).sum())} blocked/fixed reactions)"
    )
    axes[0].set_yscale("log")

    order = ["blocked", "fixed", "forward_only", "backward_only", "bidirectional", "other"]
    order = [c for c in order if c in counts.index]
    colors = {
        "blocked": "tab:gray", "fixed": "tab:orange", "forward_only": "tab:blue",
        "backward_only": "tab:red", "bidirectional": "tab:green", "other": "black",
    }
    axes[1].bar(order, [counts[c] for c in order], color=[colors[c] for c in order])
    axes[1].set_ylabel("Reaction count")
    axes[1].set_title("Reaction Direction Classification")
    for i, c in enumerate(order):
        axes[1].text(i, counts[c], str(counts[c]), ha="center", va="bottom", fontsize=9)

    fig.suptitle(
        f"Flux-Span Distribution & Direction Classification "
        f"(FVA @ {int(FRACTION_OF_OPTIMUM*100)}% of PHB optimum, loopless={LOOPLESS})",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.93])

    png_path = os.path.join(OUTPUT_DIR, "step_7.19_flux_span_and_direction.png")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    print(f"\nFigure written to: {png_path}")


if __name__ == "__main__":
    main()
