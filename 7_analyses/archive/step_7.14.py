"""
step_7.14_subsystem_flux_coupling.py

Subsystem flux coupling: identifies reactions whose fluxes are
stoichiometrically linked (they tend to move together, in a fixed or
near-fixed ratio, across the feasible flux space) — a strong signal of a
shared linear pathway segment, an obligate cofactor-regeneration loop, or
redundant/duplicate reactions.

Method (sampling-based proxy for exact flux coupling analysis):
    1. Restrict to reactions within COUPLING_MAX_DISTANCE of PHB synthesis
       (reusing step_7.12's pathway-distance ranking) plus, optionally,
       any other reactions sharing a subsystem annotation with those.
    2. Draw many random feasible flux samples from the model (via
       cobra's ACHR sampler) under a fixed, reasonable set of medium
       bounds.
    3. Compute the pairwise Pearson correlation of flux across samples
       for every reaction pair in scope.
    4. Report/plot pairs with |r| above COUPLING_THRESHOLD as
       "flux-coupled".

Note: this is a statistical proxy, not exact LP-based flux coupling
analysis (FCA) — it will reliably catch strong linear/obligate coupling
but can miss coupling that only manifests in extreme corners of the flux
cone that random sampling under-visits, and can occasionally flag
incidental correlation. Treat results as candidates to confirm by
inspecting the stoichiometry directly.

MEMORY SAFETY NOTE: an earlier version of this script sized the
correlation heatmap figure directly from the number of reactions in
scope (`0.3 * n` inches per axis). If the pathway-distance neighborhood
turns out larger than expected (e.g. a hub metabolite wasn't fully
excluded by the currency-metabolite filter), that scales to a
multi-hundred-inch figure at 200 dpi and can exhaust available RAM
before it ever gets to plotting. This version hard-caps both the number
of reactions considered (MAX_REACTIONS_IN_SCOPE) and the rendered figure
size (via step7_utils.capped_figsize), and frees the full sample matrix
from memory as soon as the relevant columns are extracted.

Requires step_7.12 to have been run first.

Outputs:
    step_7_outputs/step_7.14_flux_coupling_pairs.csv
    step_7_outputs/step_7.14_flux_coupling_heatmap.png
"""

import os
import gc

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cobra.sampling import sample

from step_7_utils import MODEL_PATH, OUTPUT_DIR, cap_item_list, capped_figsize, ensure_output_dir, load_model

COUPLING_MAX_DISTANCE = 4     # scope: reactions within this many hops of PHB synthesis
MAX_REACTIONS_IN_SCOPE = 120   # HARD CAP — see the memory-safety note in the module docstring
N_SAMPLES = 150                 # random flux samples to draw (more = slower/more RAM, more reliable)
SAMPLING_METHOD = "achr"       # 'achr' (single-threaded, low memory) or 'optgp' (parallel, more memory)
COUPLING_THRESHOLD = 0.90      # |Pearson r| above this = "coupled"

DISTANCE_CSV = os.path.join(OUTPUT_DIR, "step_7.12_phb_pathway_distances.csv")


def load_pathway_distances():
    if not os.path.exists(DISTANCE_CSV):
        raise SystemExit(
            f"Could not find {DISTANCE_CSV}. Run step_7.12_phb_pathway_distance_and_map.py "
            f"first — this script reuses its pathway-distance ranking."
        )
    return pd.read_csv(DISTANCE_CSV)


def main():
    ensure_output_dir(OUTPUT_DIR)
    model = load_model(MODEL_PATH)

    dist_df = load_pathway_distances()
    scoped_ids = list(
        dist_df.loc[dist_df["distance_from_PHB"] <= COUPLING_MAX_DISTANCE, "reaction_id"]
    )
    print(f"Found {len(scoped_ids)} reactions within {COUPLING_MAX_DISTANCE} hops of PHB synthesis.")
    scoped_ids = cap_item_list(scoped_ids, MAX_REACTIONS_IN_SCOPE, label="reactions in scope")
    print(f"Scoping flux coupling analysis to {len(scoped_ids)} reactions.")

    print(f"\nDrawing {N_SAMPLES} random flux samples via '{SAMPLING_METHOD}' "
          f"(this is the slow step)...")
    samples = sample(model, n=N_SAMPLES, method=SAMPLING_METHOD)

    scoped_ids = [rid for rid in scoped_ids if rid in samples.columns]
    # extract only the columns we need, in float32, then drop the full
    # (all-reactions x n_samples) matrix immediately to free memory
    sub_samples = samples[scoped_ids].astype("float32")
    del samples
    gc.collect()

    # Drop reactions with ~zero variance (they're not informative for
    # correlation and can produce NaNs)
    variances = sub_samples.var()
    informative_ids = variances[variances > 1e-8].index.tolist()
    dropped = set(scoped_ids) - set(informative_ids)
    if dropped:
        print(f"  (dropping {len(dropped)} reactions with ~zero flux variance in the samples)")
    sub_samples = sub_samples[informative_ids]

    print("Computing pairwise correlation matrix...")
    corr = sub_samples.corr(method="pearson")

    # --- Extract coupled pairs above threshold --------------------------------
    pairs = []
    ids = corr.columns.tolist()
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            r = corr.iloc[i, j]
            if pd.notna(r) and abs(r) >= COUPLING_THRESHOLD:
                pairs.append({
                    "reaction_1": ids[i], "reaction_2": ids[j],
                    "pearson_r": r,
                    "coupling_type": "positive (co-directional)" if r > 0 else "negative (anti-directional)",
                })
    pairs_df = pd.DataFrame(pairs).sort_values("pearson_r", key=abs, ascending=False)

    csv_path = os.path.join(OUTPUT_DIR, "step_7.14_flux_coupling_pairs.csv")
    pairs_df.to_csv(csv_path, index=False)
    print(f"\nFound {len(pairs_df)} reaction pairs with |r| >= {COUPLING_THRESHOLD}.")
    print(f"Coupled-pairs table written to: {csv_path}")
    if len(pairs_df):
        print(pairs_df.head(20).to_string(index=False))

    # --- Heatmap of the full correlation matrix in scope ----------------------
    fig_size = capped_figsize(len(ids), per_item=0.25, min_size=8.0, max_size=16.0)
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    im = ax.imshow(corr.values, cmap="coolwarm", vmin=-1, vmax=1)
    label_fontsize = 6 if len(ids) <= 80 else 4
    ax.set_xticks(range(len(ids)))
    ax.set_xticklabels(ids, rotation=90, fontsize=label_fontsize)
    ax.set_yticks(range(len(ids)))
    ax.set_yticklabels(ids, fontsize=label_fontsize)
    fig.colorbar(im, ax=ax, label="Pearson r (sampled flux correlation)")
    ax.set_title(
        f"Flux Coupling — Reactions within {COUPLING_MAX_DISTANCE} hops of PHB synthesis\n"
        f"({N_SAMPLES} samples, {SAMPLING_METHOD} sampler)"
    )
    fig.tight_layout()

    png_path = os.path.join(OUTPUT_DIR, "step_7.14_flux_coupling_heatmap.png")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    print(f"Heatmap written to: {png_path}")


if __name__ == "__main__":
    main()
