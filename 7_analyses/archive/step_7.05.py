"""
step_7.5_phenotypic_phase_plane.py

Phenotypic Phase Plane (PhPP) analysis (Edwards & Palsson, 2001-style).

Varies two uptake fluxes over a 2D grid (by default: carbon source uptake
and oxygen uptake) and computes the maximum achievable growth rate at each
grid point. The resulting surface typically breaks into a small number of
regions ("phases") with distinct linear shadow-price behavior (e.g.
carbon-limited vs. oxygen-limited vs. futile/infeasible regions).
Phase transition points/lines are approximated here by detecting where the
local slope of growth rate vs. the varied fluxes changes sharply.

Outputs:
    step_7_outputs/step_7.5_phpp_grid.csv           (full growth-rate grid)
    step_7_outputs/step_7.5_phpp_heatmap.png         (2D phenotype map)
    step_7_outputs/step_7.5_phpp_transition_points.csv
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cobra.flux_analysis import production_envelope

from step_7_utils import (
    MODEL_PATH,
    OUTPUT_DIR,
    ensure_output_dir,
    find_biomass_reaction,
    find_first_match,
    load_model,
    GLUCOSE_EXCHANGE_CANDIDATES,
    OXYGEN_EXCHANGE_CANDIDATES,
)

# --------------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------------
# Set these explicitly if auto-detection picks the wrong exchange reactions,
# e.g. CARBON_EXCHANGE_ID = "EX_glc__D_e"
CARBON_EXCHANGE_ID = "EX_glc__D_e"
OXYGEN_EXCHANGE_ID = "EX_o2_e"

POINTS = 20            # grid resolution along each axis
CARBON_UPTAKE_RANGE = (0, 20)   # mmol/gDW/h, magnitude of uptake (bound made negative)
OXYGEN_UPTAKE_RANGE = (0, 20)   # mmol/gDW/h

# Slope-change threshold (relative) used to flag phase-transition grid points
TRANSITION_THRESHOLD = 0.15


def setup_exchange_bounds(model, carbon_rxn, oxygen_rxn):
    """Open the two exchange reactions up to the ranges we want to scan;
    production_envelope will then vary within [lower_bound, upper_bound]."""
    carbon_rxn.lower_bound = -CARBON_UPTAKE_RANGE[1]
    carbon_rxn.upper_bound = -CARBON_UPTAKE_RANGE[0]
    oxygen_rxn.lower_bound = -OXYGEN_UPTAKE_RANGE[1]
    oxygen_rxn.upper_bound = -OXYGEN_UPTAKE_RANGE[0]


def find_transition_points(grid_df, x_col, y_col, z_col):
    """
    Very simple phase-transition detector: reshape the (x, y, z) samples
    onto a grid, then flag points where the local gradient of z changes
    sharply relative to its neighbors (a proxy for a shadow-price / phase
    boundary crossing). This is a heuristic, not an exact LP-based
    phase-boundary computation — use it to visually locate transitions,
    then inspect nearby points more closely if needed.
    """
    x_vals = np.sort(grid_df[x_col].unique())
    y_vals = np.sort(grid_df[y_col].unique())
    Z = grid_df.pivot(index=y_col, columns=x_col, values=z_col).reindex(
        index=y_vals, columns=x_vals
    ).values

    dzdx = np.gradient(Z, x_vals, axis=1)
    dzdy = np.gradient(Z, y_vals, axis=0)

    # second derivative magnitude as a transition indicator
    d2zdx2 = np.gradient(dzdx, x_vals, axis=1)
    d2zdy2 = np.gradient(dzdy, y_vals, axis=0)
    curvature = np.abs(d2zdx2) + np.abs(d2zdy2)

    threshold = TRANSITION_THRESHOLD * np.nanmax(curvature) if np.nanmax(curvature) > 0 else np.inf
    transition_mask = curvature > threshold

    transitions = []
    for i, y in enumerate(y_vals):
        for j, x in enumerate(x_vals):
            if transition_mask[i, j]:
                transitions.append({x_col: x, y_col: y, z_col: Z[i, j], "curvature": curvature[i, j]})

    return x_vals, y_vals, Z, transitions


def main():
    ensure_output_dir(OUTPUT_DIR)
    model = load_model(MODEL_PATH)
    biomass = find_biomass_reaction(model)

    carbon_rxn = (
        model.reactions.get_by_id(CARBON_EXCHANGE_ID)
        if CARBON_EXCHANGE_ID
        else find_first_match(model, GLUCOSE_EXCHANGE_CANDIDATES)
    )
    oxygen_rxn = (
        model.reactions.get_by_id(OXYGEN_EXCHANGE_ID)
        if OXYGEN_EXCHANGE_ID
        else find_first_match(model, OXYGEN_EXCHANGE_CANDIDATES)
    )

    if carbon_rxn is None or oxygen_rxn is None:
        raise ValueError(
            "Could not auto-detect carbon/oxygen exchange reactions. "
            "Set CARBON_EXCHANGE_ID and OXYGEN_EXCHANGE_ID at the top of this "
            "script to the correct exchange reaction IDs (e.g. 'EX_glc__D_e', 'EX_o2_e')."
        )

    print(f"Carbon exchange: {carbon_rxn.id} ({carbon_rxn.name})")
    print(f"Oxygen exchange: {oxygen_rxn.id} ({oxygen_rxn.name})")

    setup_exchange_bounds(model, carbon_rxn, oxygen_rxn)

    print("Computing 2D phenotypic phase plane (this runs FVA across a grid)...")
    grid = production_envelope(
        model,
        reactions=[carbon_rxn.id, oxygen_rxn.id],
        objective=biomass.id,
        points=POINTS,
    )

    csv_path = os.path.join(OUTPUT_DIR, "step_7.5_phpp_grid.csv")
    grid.to_csv(csv_path, index=False)

    x_vals, y_vals, Z, transitions = find_transition_points(
        grid, carbon_rxn.id, oxygen_rxn.id, "flux_maximum"
    )

    # --- Plot heatmap -----------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 6))
    mesh = ax.pcolormesh(x_vals, y_vals, Z, shading="auto", cmap="viridis")
    fig.colorbar(mesh, ax=ax, label="Max growth rate (1/h)")

    if transitions:
        tx = [t[carbon_rxn.id] for t in transitions]
        ty = [t[oxygen_rxn.id] for t in transitions]
        ax.scatter(tx, ty, color="red", s=12, marker="x", label="Candidate phase transitions")
        ax.legend(loc="upper right")

    ax.set_xlabel(f"{carbon_rxn.id} flux (mmol/gDW/h, uptake shown as negative bound)")
    ax.set_ylabel(f"{oxygen_rxn.id} flux (mmol/gDW/h, uptake shown as negative bound)")
    ax.set_title("Phenotypic Phase Plane — Growth Rate")
    fig.tight_layout()

    png_path = os.path.join(OUTPUT_DIR, "step_7.5_phpp_heatmap.png")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    # --- Transition points CSV ---------------------------------------------------
    trans_csv = os.path.join(OUTPUT_DIR, "step_7.5_phpp_transition_points.csv")
    import csv as _csv
    with open(trans_csv, "w", newline="") as fh:
        writer = _csv.writer(fh)
        writer.writerow([carbon_rxn.id, oxygen_rxn.id, "growth_rate", "curvature"])
        for t in transitions:
            writer.writerow([t[carbon_rxn.id], t[oxygen_rxn.id], t["flux_maximum"], t["curvature"]])

    print(f"\nGrid data written to:        {csv_path}")
    print(f"Heatmap written to:          {png_path}")
    print(f"Candidate transition points: {trans_csv} ({len(transitions)} flagged)")
    print(
        "\nNote: transition points are a heuristic (local curvature of the growth "
        "surface). For exact phase boundaries, examine shadow prices / reduced "
        "costs of the two exchange reactions at neighboring grid points."
    )


if __name__ == "__main__":
    main()
