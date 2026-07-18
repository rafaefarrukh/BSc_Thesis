"""
step_7.17_random_sampling_achr.py

Random flux sampling: draws many random points from the model's feasible
flux polytope using cobra's ACHR (Artificial Centering Hit-and-Run)
sampler, under the model's current (default, as-loaded) bounds. This
"maps" the feasible flux space in a way that a single FBA/pFBA solution
cannot — FBA gives one optimal vertex, sampling shows the shape and
spread of everything that's actually feasible.

From the samples this script reports:
    - per-reaction summary statistics (mean, std, min, max flux)
    - the growth-rate and PHB-flux marginal distributions
    - a 2D PCA projection of the sampled flux vectors (a rough visual
      map of the feasible space's shape/structure)

PCA is implemented directly via numpy SVD (no scikit-learn dependency).

Outputs:
    step_7_outputs/step_7.17_samples_raw.csv            (full sample matrix)
    step_7_outputs/step_7.17_samples_summary_stats.csv   (per-reaction stats)
    step_7_outputs/step_7.17_sampling_overview.png
"""

import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cobra.sampling import sample

from step_7_utils import (
    MODEL_PATH,
    OUTPUT_DIR,
    describe_reaction,
    ensure_output_dir,
    find_biomass_reaction,
    find_phb_reaction,
    load_model,
)

PHB_REACTION_ID = "PHBS_syn_1"
N_SAMPLES = 300
SAMPLING_METHOD = "achr"   # 'achr' (no extra deps) or 'optgp' (parallel, more memory)
SAVE_FULL_SAMPLE_MATRIX = True


def pca_2d(matrix):
    """Simple PCA via SVD. matrix: (n_samples, n_features), already
    column-centered internally. Returns (n_samples, 2) projection and the
    fraction of variance explained by each of the first two components."""
    X = matrix - matrix.mean(axis=0, keepdims=True)
    # drop zero-variance columns to avoid numerical issues
    keep = X.std(axis=0) > 1e-10
    X = X[:, keep]
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    projection = U[:, :2] * S[:2]
    explained = (S**2) / np.sum(S**2)
    return projection, explained[:2]


def main():
    ensure_output_dir(OUTPUT_DIR)
    model = load_model(MODEL_PATH)
    biomass = find_biomass_reaction(model)
    phb = model.reactions.get_by_id(PHB_REACTION_ID) if PHB_REACTION_ID else find_phb_reaction(model)

    describe_reaction(biomass, "Biomass reaction")
    describe_reaction(phb, "PHB reaction")

    print(f"\nDrawing {N_SAMPLES} random flux samples via '{SAMPLING_METHOD}' "
          f"under the model's current default bounds (this can take a while)...")
    samples = sample(model, n=N_SAMPLES, method=SAMPLING_METHOD)

    if SAVE_FULL_SAMPLE_MATRIX:
        raw_csv = os.path.join(OUTPUT_DIR, "step_7.17_samples_raw.csv")
        samples.to_csv(raw_csv, index=False)
        print(f"Full sample matrix written to: {raw_csv}")

    # --- Per-reaction summary stats -------------------------------------------
    summary = pd.DataFrame({
        "reaction_id": samples.columns,
        "mean_flux": samples.mean(axis=0).values,
        "std_flux": samples.std(axis=0).values,
        "min_flux": samples.min(axis=0).values,
        "max_flux": samples.max(axis=0).values,
    })
    summary_csv = os.path.join(OUTPUT_DIR, "step_7.17_samples_summary_stats.csv")
    summary.to_csv(summary_csv, index=False)
    print(f"Per-reaction summary stats written to: {summary_csv}")

    # --- PCA projection ----------------------------------------------------------
    print("Computing PCA projection of sampled flux vectors...")
    projection, explained = pca_2d(samples.values)

    # --- Figure: PCA scatter + growth/PHB marginal histograms ----------------
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    color_vals = samples[phb.id].values if phb is not None and phb.id in samples.columns else None
    sc = axes[0].scatter(projection[:, 0], projection[:, 1], c=color_vals, cmap="viridis", s=15, alpha=0.7)
    if color_vals is not None:
        fig.colorbar(sc, ax=axes[0], label=f"{phb.id} flux")
    axes[0].set_xlabel(f"PC1 ({explained[0]*100:.1f}% var)")
    axes[0].set_ylabel(f"PC2 ({explained[1]*100:.1f}% var)")
    axes[0].set_title("PCA of Sampled Feasible Flux Space")

    if biomass.id in samples.columns:
        axes[1].hist(samples[biomass.id], bins=30, color="tab:blue", alpha=0.8)
        axes[1].set_xlabel(f"{biomass.id} flux (growth rate, 1/h)")
        axes[1].set_ylabel("Sample count")
        axes[1].set_title("Growth Rate Distribution")

    if phb is not None and phb.id in samples.columns:
        axes[2].hist(samples[phb.id], bins=30, color="tab:green", alpha=0.8)
        axes[2].set_xlabel(f"{phb.id} flux")
        axes[2].set_ylabel("Sample count")
        axes[2].set_title("PHB Flux Distribution")

    fig.suptitle(f"Random Flux Sampling ({SAMPLING_METHOD.upper()}, n={N_SAMPLES})", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    png_path = os.path.join(OUTPUT_DIR, "step_7.17_sampling_overview.png")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    print(f"\nSampling overview figure written to: {png_path}")


if __name__ == "__main__":
    main()
