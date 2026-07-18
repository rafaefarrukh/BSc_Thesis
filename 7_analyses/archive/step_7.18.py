"""
step_7.18_fva_multifraction.py

Flux Variability Analysis (FVA) at three fractions of the optimal PHB
flux (100%, 90%, 50%), computed both the regular way and "loopless"
(excluding thermodynamically infeasible internal flux loops — relevant
here since the MEMOTE report flagged 134 reactions in stoichiometrically
balanced cycles). For each of the 6 (fraction x loopless) combinations
this reports summary statistics, and separately plots how the flux range
of reactions in the PHB pathway neighborhood (from step_7.12) changes as
the fraction is relaxed from 100% down to 50%.

By default the FVA scope is restricted to the PHB pathway neighborhood
(reusing step_7.12's distance ranking) to keep runtime reasonable — set
REACTION_SCOPE = "all" below to run on every reaction in the model
instead (much slower: FVA solves roughly 2 LPs per reaction per
fraction/loopless combination).

Requires step_7.12 to have been run first if REACTION_SCOPE == "pathway".

Outputs:
    step_7_outputs/step_7.18_fva_frac<pct>_<loopless|standard>.csv   (x6)
    step_7_outputs/step_7.18_fva_summary_stats.csv
    step_7_outputs/step_7.18_pathway_ranges_across_fractions.png
"""

import os

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cobra.flux_analysis import flux_variability_analysis

from step_7_utils import (
    MODEL_PATH,
    OUTPUT_DIR,
    cap_item_list,
    capped_figsize,
    describe_reaction,
    ensure_output_dir,
    find_phb_reaction,
    load_model,
)

PHB_REACTION_ID = "PHBS_syn_1"
FRACTIONS = [1.0, 0.9, 0.5]

# "pathway" = restrict to reactions within PATHWAY_MAX_DISTANCE hops of PHB
# synthesis (reuses step_7.12's output); "all" = every reaction in the model
REACTION_SCOPE = "pathway"
PATHWAY_MAX_DISTANCE = 6
MAX_REACTIONS_IN_SCOPE = 150   # hard cap — see step7_utils.capped_figsize for why this matters

DISTANCE_CSV = os.path.join(OUTPUT_DIR, "step_7.12_phb_pathway_distances.csv")


def get_reaction_scope(model):
    if REACTION_SCOPE == "all":
        ids = [r.id for r in model.reactions if not r.boundary]
        return cap_item_list(ids, MAX_REACTIONS_IN_SCOPE, label="reactions"), None

    if not os.path.exists(DISTANCE_CSV):
        print(f"WARNING: {DISTANCE_CSV} not found (run step_7.12 first). "
              f"Falling back to REACTION_SCOPE='all' for this run.")
        ids = [r.id for r in model.reactions if not r.boundary]
        return cap_item_list(ids, MAX_REACTIONS_IN_SCOPE, label="reactions"), None

    dist_df = pd.read_csv(DISTANCE_CSV)
    scoped = dist_df[dist_df["distance_from_PHB"] <= PATHWAY_MAX_DISTANCE].sort_values("distance_from_PHB")
    ids = cap_item_list(scoped["reaction_id"].tolist(), MAX_REACTIONS_IN_SCOPE, label="pathway reactions")
    scoped = scoped[scoped["reaction_id"].isin(ids)]
    return ids, scoped


def summarize_fva(fva_df, eps=1e-9):
    ranges = fva_df["maximum"] - fva_df["minimum"]
    blocked = ((fva_df["minimum"].abs() < eps) & (fva_df["maximum"].abs() < eps)).sum()
    return {
        "n_reactions": len(fva_df),
        "n_blocked": int(blocked),
        "pct_blocked": 100 * blocked / len(fva_df) if len(fva_df) else 0.0,
        "mean_range": ranges.mean(),
        "median_range": ranges.median(),
        "total_range": ranges.sum(),
        "max_range": ranges.max(),
    }


def main():
    ensure_output_dir(OUTPUT_DIR)
    model = load_model(MODEL_PATH)
    phb = model.reactions.get_by_id(PHB_REACTION_ID) if PHB_REACTION_ID else find_phb_reaction(model)
    if phb is None:
        raise ValueError("Could not auto-detect a PHB reaction. Set PHB_REACTION_ID at the top of this script.")
    describe_reaction(phb, "PHB reaction (FVA objective)")

    reaction_ids, dist_df = get_reaction_scope(model)
    print(f"FVA scope: {len(reaction_ids)} reactions "
          f"({'pathway neighborhood' if REACTION_SCOPE == 'pathway' and dist_df is not None else 'all reactions'})")

    model.objective = phb.id

    summary_rows = []
    fva_results = {}   # (fraction, loopless) -> DataFrame
    for fraction in FRACTIONS:
        for loopless in (False, True):
            label = f"frac{int(fraction*100)}_{'loopless' if loopless else 'standard'}"
            print(f"\nRunning FVA: fraction_of_optimum={fraction}, loopless={loopless} ...")
            fva = flux_variability_analysis(
                model, reaction_list=reaction_ids, fraction_of_optimum=fraction, loopless=loopless
            )
            fva_results[(fraction, loopless)] = fva

            csv_path = os.path.join(OUTPUT_DIR, f"step_7.18_fva_{label}.csv")
            fva.to_csv(csv_path)
            print(f"  written to {csv_path}")

            stats = summarize_fva(fva)
            stats.update({"fraction_of_optimum": fraction, "loopless": loopless})
            summary_rows.append(stats)

    summary_df = pd.DataFrame(summary_rows)[
        ["fraction_of_optimum", "loopless", "n_reactions", "n_blocked", "pct_blocked",
         "mean_range", "median_range", "total_range", "max_range"]
    ]
    summary_csv = os.path.join(OUTPUT_DIR, "step_7.18_fva_summary_stats.csv")
    summary_df.to_csv(summary_csv, index=False)
    print(f"\nFVA summary statistics (all fraction x loopless combinations):")
    print(summary_df.to_string(index=False))
    print(f"\nSummary stats written to: {summary_csv}")

    # --- PHB pathway flux ranges across fractions (loopless variant) ---------
    if dist_df is not None:
        fig_width = capped_figsize(len(reaction_ids), per_item=0.2, min_size=10.0, max_size=24.0)
        fig, ax = plt.subplots(figsize=(fig_width, 6))
        ordered_ids = dist_df.sort_values(["distance_from_PHB", "reaction_id"])["reaction_id"].tolist()
        ordered_ids = [rid for rid in ordered_ids if rid in reaction_ids]

        width = 0.25
        x = np.arange(len(ordered_ids))
        colors = {1.0: "tab:blue", 0.9: "tab:orange", 0.5: "tab:red"}
        for i, fraction in enumerate(FRACTIONS):
            fva = fva_results[(fraction, True)]  # loopless variant
            ranges = (fva.loc[ordered_ids, "maximum"] - fva.loc[ordered_ids, "minimum"]).values
            ax.bar(x + (i - 1) * width, ranges, width, label=f"{int(fraction*100)}% of PHB optimum",
                   color=colors.get(fraction, None))

        ax.set_xticks(x)
        label_fontsize = 6 if len(ordered_ids) <= 100 else 4
        ax.set_xticklabels(ordered_ids, rotation=90, fontsize=label_fontsize)
        ax.set_ylabel("Flux range (max - min, loopless FVA)")
        ax.set_title("PHB Pathway Flux Ranges Across PHB-Optimality Fractions")
        ax.legend()
        fig.tight_layout()

        png_path = os.path.join(OUTPUT_DIR, "step_7.18_pathway_ranges_across_fractions.png")
        fig.savefig(png_path, dpi=200)
        plt.close(fig)
        print(f"Pathway-ranges-across-fractions figure written to: {png_path}")
    else:
        print("\n(Skipping pathway-ranges figure: no pathway-distance scoping available. "
              "Run step_7.12 first and set REACTION_SCOPE='pathway' to get this plot.)")


if __name__ == "__main__":
    main()
