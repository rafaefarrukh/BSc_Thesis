"""
step_7.9_2D_robustness_heatmaps.py

2D robustness analysis: for two nutrient pairs — (O2 x glucose) and
(nitrogen x glucose) — sweeps both uptake bounds over a grid and records,
at each grid point:
    - maximum growth rate
    - growth-coupled PHB flux (growth fixed at 90% of that point's max
      growth, then PHB maximized)

This shows how the two nutrients jointly shape the feasible growth/PHB
phenotype space (complementing the 1D sweeps in step_7.8).

Outputs:
    step_7_outputs/step_7.9_2D_o2_glucose.csv
    step_7_outputs/step_7.9_2D_nitrogen_glucose.csv
    step_7_outputs/step_7.9_2D_robustness_heatmaps.png
"""

import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from step_7_utils import (
    MODEL_PATH,
    OUTPUT_DIR,
    describe_reaction,
    ensure_output_dir,
    find_biomass_reaction,
    find_first_match,
    find_phb_reaction,
    load_model,
    GLUCOSE_EXCHANGE_CANDIDATES,
    NITROGEN_EXCHANGE_CANDIDATES,
    OXYGEN_EXCHANGE_CANDIDATES,
    MALTOSE_EXCHANGE_ID,
)

# --------------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------------
PHB_REACTION_ID = "PHBS_syn_1"
CARBON_EXCHANGE_ID = "EX_glc__D_e"
NITROGEN_EXCHANGE_ID = "EX_nh4_e"
OXYGEN_EXCHANGE_ID = "EX_o2_e"
MALTOSE_ID = MALTOSE_EXCHANGE_ID  # "EX_malt_e" — included as a 3rd robustness grid

GROWTH_COUPLED_FRACTION = 0.90
GRID_POINTS = 15  # per axis -> GRID_POINTS^2 FBA solves per heatmap pair

GLC_RANGE = (0.0, 20.0)
O2_RANGE = (0.0, 25.0)
NH4_RANGE = (0.0, 10.0)
MALTOSE_RANGE = (0.0, 20.0)
NH4_NONLIMITING = 20.0  # background N level used in the O2 x glucose grid
O2_NONLIMITING = 20.0   # background O2 level used in the nitrogen x maltose grid


def get_reactions(model):
    biomass = find_biomass_reaction(model)
    phb = model.reactions.get_by_id(PHB_REACTION_ID) if PHB_REACTION_ID else find_phb_reaction(model)
    carbon = (
        model.reactions.get_by_id(CARBON_EXCHANGE_ID)
        if CARBON_EXCHANGE_ID
        else find_first_match(model, GLUCOSE_EXCHANGE_CANDIDATES)
    )
    nitrogen = (
        model.reactions.get_by_id(NITROGEN_EXCHANGE_ID)
        if NITROGEN_EXCHANGE_ID
        else find_first_match(model, NITROGEN_EXCHANGE_CANDIDATES)
    )
    oxygen = (
        model.reactions.get_by_id(OXYGEN_EXCHANGE_ID)
        if OXYGEN_EXCHANGE_ID
        else find_first_match(model, OXYGEN_EXCHANGE_CANDIDATES)
    )
    maltose = model.reactions.get_by_id(MALTOSE_ID) if MALTOSE_ID in model.reactions else None
    missing = [
        name for name, rxn in
        [("PHB", phb), ("carbon", carbon), ("nitrogen", nitrogen), ("oxygen", oxygen)]
        if rxn is None
    ]
    if missing:
        raise ValueError(
            f"Could not auto-detect reaction(s): {missing}. Set the corresponding "
            f"*_REACTION_ID / *_EXCHANGE_ID constants at the top of this script."
        )
    if maltose is None:
        print(f"  NOTE: maltose exchange '{MALTOSE_ID}' not found in model — "
              f"nitrogen x maltose grid will be skipped.")
    return biomass, phb, carbon, nitrogen, oxygen, maltose


def growth_and_phb_at_point(model, biomass, phb, carbon, nitrogen, oxygen, maltose,
                             glc_bound, o2_bound, nh4_bound, malt_bound=0.0):
    with model as m:
        m.reactions.get_by_id(carbon.id).lower_bound = -glc_bound
        m.reactions.get_by_id(oxygen.id).lower_bound = -o2_bound
        m.reactions.get_by_id(nitrogen.id).lower_bound = -nh4_bound
        if maltose is not None:
            m.reactions.get_by_id(maltose.id).lower_bound = -malt_bound

        m.objective = biomass.id
        max_growth = m.slim_optimize()
        max_growth = max_growth if max_growth and max_growth > 1e-9 else 0.0

        phb_flux = 0.0
        if max_growth > 0:
            fixed_growth = GROWTH_COUPLED_FRACTION * max_growth
            biomass_m = m.reactions.get_by_id(biomass.id)
            biomass_m.lower_bound = fixed_growth
            biomass_m.upper_bound = fixed_growth
            m.objective = phb.id
            phb_flux = m.slim_optimize()
            phb_flux = phb_flux if phb_flux and phb_flux > 1e-9 else 0.0

    return max_growth, phb_flux


def sweep_2d(model, biomass, phb, carbon, nitrogen, oxygen, maltose,
             x_label, x_range, y_label, y_range, fixed_third):
    """Generic 2D sweep. x_label/y_label in {'glucose','oxygen','nitrogen','maltose'};
    fixed_third: dict with the bound to use for whichever nutrient isn't
    part of this particular grid."""
    x_vals = np.linspace(x_range[0], x_range[1], GRID_POINTS)
    y_vals = np.linspace(y_range[0], y_range[1], GRID_POINTS)

    rows = []
    for yv in y_vals:
        for xv in x_vals:
            bounds = dict(fixed_third)
            bounds[x_label] = xv
            bounds[y_label] = yv
            growth, phb_flux = growth_and_phb_at_point(
                model, biomass, phb, carbon, nitrogen, oxygen, maltose,
                bounds.get("glucose", 0.0), bounds.get("oxygen", 0.0),
                bounds.get("nitrogen", 0.0), malt_bound=bounds.get("maltose", 0.0),
            )
            rows.append({x_label: xv, y_label: yv, "growth_rate_1_h": growth, "PHB_flux": phb_flux})
    return pd.DataFrame(rows), x_vals, y_vals


def to_grid(df, x_label, y_label, z_col, x_vals, y_vals):
    return df.pivot(index=y_label, columns=x_label, values=z_col).reindex(
        index=y_vals, columns=x_vals
    ).values


def main():
    ensure_output_dir(OUTPUT_DIR)
    model = load_model(MODEL_PATH)
    biomass, phb, carbon, nitrogen, oxygen, maltose = get_reactions(model)

    describe_reaction(biomass, "Biomass reaction")
    describe_reaction(phb, "PHB reaction")
    if maltose is not None:
        describe_reaction(maltose, "Maltose exchange")

    print(f"\nRunning O2 x glucose grid ({GRID_POINTS}x{GRID_POINTS} = {GRID_POINTS**2} FBA points)...")
    o2_glc_df, o2_x, o2_y = sweep_2d(
        model, biomass, phb, carbon, nitrogen, oxygen, maltose,
        "glucose", GLC_RANGE, "oxygen", O2_RANGE,
        fixed_third={"glucose": 0, "oxygen": 0, "nitrogen": NH4_NONLIMITING, "maltose": 0},
    )
    o2_glc_csv = os.path.join(OUTPUT_DIR, "step_7.9_2D_o2_glucose.csv")
    o2_glc_df.to_csv(o2_glc_csv, index=False)
    print(f"  written to {o2_glc_csv}")

    print(f"\nRunning nitrogen x glucose grid ({GRID_POINTS}x{GRID_POINTS} = {GRID_POINTS**2} FBA points)...")
    n_glc_df, n_x, n_y = sweep_2d(
        model, biomass, phb, carbon, nitrogen, oxygen, maltose,
        "glucose", GLC_RANGE, "nitrogen", NH4_RANGE,
        fixed_third={"glucose": 0, "oxygen": O2_RANGE[1], "nitrogen": 0, "maltose": 0},
    )
    n_glc_csv = os.path.join(OUTPUT_DIR, "step_7.9_2D_nitrogen_glucose.csv")
    n_glc_df.to_csv(n_glc_csv, index=False)
    print(f"  written to {n_glc_csv}")

    # --- Nitrogen x maltose grid (glucose closed to isolate maltose as the
    # sole carbon source) --------------------------------------------------
    n_malt_df = n_malt_x = n_malt_y = None
    if maltose is not None:
        print(f"\nRunning nitrogen x maltose grid ({GRID_POINTS}x{GRID_POINTS} = {GRID_POINTS**2} FBA points)...")
        n_malt_df, n_malt_x, n_malt_y = sweep_2d(
            model, biomass, phb, carbon, nitrogen, oxygen, maltose,
            "maltose", MALTOSE_RANGE, "nitrogen", NH4_RANGE,
            fixed_third={"glucose": 0, "oxygen": O2_NONLIMITING, "nitrogen": 0, "maltose": 0},
        )
        n_malt_csv = os.path.join(OUTPUT_DIR, "step_7.9_2D_nitrogen_maltose.csv")
        n_malt_df.to_csv(n_malt_csv, index=False)
        print(f"  written to {n_malt_csv}")
    else:
        print("\nSkipping nitrogen x maltose grid (maltose exchange not found in model).")

    # --- Combined heatmap figure (2x2, or 3x2 if maltose grid available) -----
    n_rows = 3 if maltose is not None else 2
    fig, axes = plt.subplots(n_rows, 2, figsize=(12, 5 * n_rows))

    growth_o2 = to_grid(o2_glc_df, "glucose", "oxygen", "growth_rate_1_h", o2_x, o2_y)
    phb_o2 = to_grid(o2_glc_df, "glucose", "oxygen", "PHB_flux", o2_x, o2_y)
    growth_n = to_grid(n_glc_df, "glucose", "nitrogen", "growth_rate_1_h", n_x, n_y)
    phb_n = to_grid(n_glc_df, "glucose", "nitrogen", "PHB_flux", n_x, n_y)

    panels = [
        (axes[0, 0], growth_o2, o2_x, o2_y, "Growth rate (1/h): O2 x Glucose", "Oxygen uptake", "Glucose uptake"),
        (axes[0, 1], phb_o2, o2_x, o2_y, "PHB flux: O2 x Glucose", "Oxygen uptake", "Glucose uptake"),
        (axes[1, 0], growth_n, n_x, n_y, "Growth rate (1/h): Nitrogen x Glucose", "NH4 uptake", "Glucose uptake"),
        (axes[1, 1], phb_n, n_x, n_y, "PHB flux: Nitrogen x Glucose", "NH4 uptake", "Glucose uptake"),
    ]
    if maltose is not None:
        growth_malt = to_grid(n_malt_df, "maltose", "nitrogen", "growth_rate_1_h", n_malt_x, n_malt_y)
        phb_malt = to_grid(n_malt_df, "maltose", "nitrogen", "PHB_flux", n_malt_x, n_malt_y)
        panels.extend([
            (axes[2, 0], growth_malt, n_malt_x, n_malt_y,
             "Growth rate (1/h): Nitrogen x Maltose", "NH4 uptake", "Maltose uptake"),
            (axes[2, 1], phb_malt, n_malt_x, n_malt_y,
             "PHB flux: Nitrogen x Maltose", "NH4 uptake", "Maltose uptake"),
        ])

    for ax, Z, xv, yv, title, ylabel, xlabel in panels:
        mesh = ax.pcolormesh(xv, yv, Z, shading="auto", cmap="viridis")
        fig.colorbar(mesh, ax=ax)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(f"{xlabel} (mmol/gDW/h)")
        ax.set_ylabel(f"{ylabel} (mmol/gDW/h)")

    fig.suptitle("2D Robustness Heatmaps", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    png_path = os.path.join(OUTPUT_DIR, "step_7.9_2D_robustness_heatmaps.png")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    print(f"\nCombined heatmap figure written to: {png_path}")


if __name__ == "__main__":
    main()
