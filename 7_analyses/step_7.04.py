import os
import csv as _csv
import numpy as np
import pandas as pd
import cobra
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from step_7_config import (
    MODEL_PATH,
    OUTPUT_DIR,
    BIOMASS_ID,
    GLUCOSE_EXCHANGE_ID,
    MALTOSE_EXCHANGE_ID,
    OXYGEN_EXCHANGE_ID,
)

POINTS = 20
CARBON_UPTAKE_RANGE = (0, 20)
OXYGEN_UPTAKE_RANGE = (0, 20)
TRANSITION_THRESHOLD = 0.15

def setup_exchange_bounds(carbon_rxn, oxygen_rxn, x, y):
    """Sets bounds for a specific grid point."""
    carbon_rxn.lower_bound = -x
    carbon_rxn.upper_bound = -x
    oxygen_rxn.lower_bound = -y
    oxygen_rxn.upper_bound = -y

def find_transition_points(grid_df, x_col, y_col, z_col):
    x_vals = np.sort(grid_df[x_col].unique())
    y_vals = np.sort(grid_df[y_col].unique())
    Z = grid_df.pivot(index=y_col, columns=x_col, values=z_col).reindex(
        index=y_vals, columns=x_vals
    ).values

    dzdx = np.gradient(Z, x_vals, axis=1)
    dzdy = np.gradient(Z, y_vals, axis=0)
    d2zdx2 = np.gradient(dzdx, x_vals, axis=1)
    d2zdy2 = np.gradient(dzdy, y_vals, axis=0)
    curvature = np.abs(d2zdx2) + np.abs(d2zdy2)

    nan_mask = np.isnan(Z)
    contaminated = nan_mask.copy()
    contaminated[:, 1:] |= nan_mask[:, :-1]
    contaminated[:, :-1] |= nan_mask[:, 1:]
    contaminated[1:, :] |= nan_mask[:-1, :]
    contaminated[:-1, :] |= nan_mask[1:, :]

    clean_curvature = np.where(contaminated, np.nan, curvature)
    max_clean = np.nanmax(clean_curvature) if np.any(~np.isnan(clean_curvature)) else 0
    threshold = TRANSITION_THRESHOLD * max_clean if max_clean > 0 else np.inf
    transition_mask = (clean_curvature > threshold) & ~contaminated

    transitions = []
    for i, y in enumerate(y_vals):
        for j, x in enumerate(x_vals):
            if transition_mask[i, j]:
                transitions.append({x_col: x, y_col: y, z_col: Z[i, j], "curvature": curvature[i, j]})
    return x_vals, y_vals, Z, transitions

def generate_phpp_data(model, biomass, carbon_rxn, oxygen_rxn, label):
    print(f"\n--- Processing PhPP data for {label.upper()} ---")

    x_vals = np.linspace(CARBON_UPTAKE_RANGE[0], CARBON_UPTAKE_RANGE[1], POINTS)
    y_vals = np.linspace(OXYGEN_UPTAKE_RANGE[0], OXYGEN_UPTAKE_RANGE[1], POINTS)

    results = []
    print(f"Computing 2D phenotypic phase plane (manual grid, {POINTS*POINTS} points)...")

    for x in x_vals:
        for y in y_vals:
            with model:
                setup_exchange_bounds(carbon_rxn, oxygen_rxn, x, y)
                growth = model.slim_optimize()
                results.append({carbon_rxn.id: x, oxygen_rxn.id: y, "flux_maximum": growth})

    grid = pd.DataFrame(results)
    csv_path = os.path.join(OUTPUT_DIR, f"step_7.04_{label}_phpp_grid.csv")
    grid.to_csv(csv_path, index=False)

    x_v, y_v, Z, transitions = find_transition_points(grid, carbon_rxn.id, oxygen_rxn.id, "flux_maximum")

    trans_csv = os.path.join(OUTPUT_DIR, f"step_7.04_{label}_phpp_transition_points.csv")
    with open(trans_csv, "w", newline="") as fh:
        writer = _csv.writer(fh)
        writer.writerow([carbon_rxn.id, oxygen_rxn.id, "growth_rate", "curvature"])
        for t in transitions:
            writer.writerow([t[carbon_rxn.id], t[oxygen_rxn.id], t["flux_maximum"], t["curvature"]])

    return x_v, y_v, Z, transitions, carbon_rxn.id, oxygen_rxn.id

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    model = cobra.io.read_sbml_model(MODEL_PATH)
    biomass = model.reactions.get_by_id(BIOMASS_ID)
    oxygen_rxn = model.reactions.get_by_id(OXYGEN_EXCHANGE_ID)
    glucose_rxn = model.reactions.get_by_id(GLUCOSE_EXCHANGE_ID)
    maltose_rxn = model.reactions.get_by_id(MALTOSE_EXCHANGE_ID)

    glc_data = generate_phpp_data(model, biomass, glucose_rxn, oxygen_rxn, "glucose")
    mal_data = generate_phpp_data(model, biomass, maltose_rxn, oxygen_rxn, "maltose")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    for item in [{"data": glc_data, "ax": ax1, "title": "Glucose Phenotypic Phase Plane"},
                 {"data": mal_data, "ax": ax2, "title": "Maltose Phenotypic Phase Plane"}]:
        x_v, y_v, Z, transitions, carbon_id, oxygen_id = item["data"]
        mesh = item["ax"].pcolormesh(x_v, y_v, Z, shading="auto", cmap="viridis")
        fig.colorbar(mesh, ax=item["ax"], label="Max growth rate (1/h)")
        if transitions:
            item["ax"].scatter([t[carbon_id] for t in transitions], [t[oxygen_id] for t in transitions],
                               color="red", s=12, marker="x", label="Candidate phase transitions")
            item["ax"].legend(loc="upper right")
        item["ax"].set_xlabel(f"{carbon_id} flux")
        item["ax"].set_ylabel(f"{oxygen_id} flux")
        item["ax"].set_title(item["title"])

    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "step_7.04_phpp_combined_heatmap.png"), dpi=200)
    plt.close(fig)

if __name__ == "__main__":
    main()
