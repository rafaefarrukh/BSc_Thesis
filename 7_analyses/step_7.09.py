"""
step_7.09.py

Aeration x carbon-source screen: combines a set of aeration levels
(anaerobic / microaerobic / aerobic / hyperaerobic O2 uptake bounds) with
a set of alternate carbon sources (whichever of the candidates in
step_7_config.ALTERNATE_CARBON_SOURCES actually exist as exchange reactions
in this model) and computes, for every combination:
    - maximum growth rate
    - growth-coupled PHB flux

This identifies which carbon source / aeration combinations support the
best growth and/or PHB production.

Outputs:
    step_7_outputs/step_7.09_aeration_carbon_screen.csv
    step_7_outputs/step_7.09_aeration_carbon_source_summary.png
"""

import os
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
    PHB_REACTION_ID,
    NITROGEN_EXCHANGE_ID,
    OXYGEN_EXCHANGE_ID,
    ALTERNATE_CARBON_SOURCES,
    describe_reaction,
    growth_coupled_phb,
)

# --------------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------------
GROWTH_COUPLED_FRACTION = 0.90
CARBON_UPTAKE_BOUND = 20.0   # mmol/gDW/h, applied to whichever carbon source is active
NH4_NONLIMITING = 20.0

AERATION_LEVELS = {
    "Anaerobic": 0.0,
    "Microaerobic": 2.0,
    "Aerobic": 20.0,
    "Hyperaerobic": 40.0,
}


def available_carbon_sources(model):
    available = {}
    for name, rxn_id in ALTERNATE_CARBON_SOURCES.items():
        if rxn_id in model.reactions:
            available[name] = rxn_id
        else:
            print(f"  (skipping '{name}': exchange reaction '{rxn_id}' not found in model)")
    if not available:
        raise ValueError(
            "None of the candidate carbon sources in ALTERNATE_CARBON_SOURCES were "
            "found in this model. Edit that dict in step_7_config.py to match this "
            "model's actual exchange reaction IDs."
        )
    return available


def evaluate_combination(model, biomass, phb, nitrogen, oxygen, carbon_rxn_id, o2_bound):
    with model as m:
        # close all other carbon-source exchanges to isolate this one,
        # then open only the target carbon source
        for name, rxn_id in ALTERNATE_CARBON_SOURCES.items():
            if rxn_id in m.reactions:
                m.reactions.get_by_id(rxn_id).lower_bound = 0.0

        m.reactions.get_by_id(carbon_rxn_id).lower_bound = -CARBON_UPTAKE_BOUND
        m.reactions.get_by_id(oxygen.id).lower_bound = -o2_bound
        m.reactions.get_by_id(nitrogen.id).lower_bound = -NH4_NONLIMITING

        growth, phb_flux = growth_coupled_phb(
            m, biomass=biomass, phb=phb, fraction=GROWTH_COUPLED_FRACTION
        )

    return growth, phb_flux


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    model = cobra.io.read_sbml_model(MODEL_PATH)

    biomass = model.reactions.get_by_id(BIOMASS_ID)
    phb = model.reactions.get_by_id(PHB_REACTION_ID)
    nitrogen = model.reactions.get_by_id(NITROGEN_EXCHANGE_ID)
    oxygen = model.reactions.get_by_id(OXYGEN_EXCHANGE_ID)

    describe_reaction(biomass, "Biomass reaction")
    describe_reaction(phb, "PHB reaction")

    print("\nChecking which candidate carbon sources exist in this model...")
    carbon_sources = available_carbon_sources(model)
    print(f"  found: {list(carbon_sources.keys())}")

    rows = []
    for c_name, c_rxn_id in carbon_sources.items():
        for aer_name, o2_bound in AERATION_LEVELS.items():
            growth, phb_flux = evaluate_combination(model, biomass, phb, nitrogen, oxygen, c_rxn_id, o2_bound)
            rows.append({
                "carbon_source": c_name,
                "aeration": aer_name,
                "o2_uptake_bound": o2_bound,
                "max_growth_rate_1_h": growth,
                "growth_coupled_PHB_flux": phb_flux,
            })
            print(f"  {c_name:14s} x {aer_name:14s} -> growth={growth:.4f} /h, PHB={phb_flux:.4f}")

    df = pd.DataFrame(rows)
    csv_path = os.path.join(OUTPUT_DIR, "step_7.09_aeration_carbon_screen.csv")
    df.to_csv(csv_path, index=False)
    print(f"\nFull results written to: {csv_path}")

    # --- Summary heatmaps: growth and PHB, carbon source x aeration --------
    growth_pivot = df.pivot(index="carbon_source", columns="aeration", values="max_growth_rate_1_h")
    phb_pivot = df.pivot(index="carbon_source", columns="aeration", values="growth_coupled_PHB_flux")
    aeration_order = list(AERATION_LEVELS.keys())
    growth_pivot = growth_pivot[aeration_order]
    phb_pivot = phb_pivot[aeration_order]

    fig, axes = plt.subplots(1, 2, figsize=(13, 0.6 * len(carbon_sources) + 3))
    for ax, pivot, title, cmap in [
        (axes[0], growth_pivot, "Max Growth Rate (1/h)", "Blues"),
        (axes[1], phb_pivot, "Growth-Coupled PHB Flux", "Greens"),
    ]:
        im = ax.imshow(pivot.values, cmap=cmap, aspect="auto")
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels(pivot.columns, rotation=30, ha="right")
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels(pivot.index)
        ax.set_title(title)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        for i in range(pivot.shape[0]):
            for j in range(pivot.shape[1]):
                val = pivot.values[i, j]
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        color="white" if val > np.nanmax(pivot.values) * 0.5 else "black", fontsize=8)

    fig.suptitle("Aeration x Carbon-Source Screen", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    png_path = os.path.join(OUTPUT_DIR, "step_7.09_aeration_carbon_source_summary.png")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    print(f"Summary figure written to: {png_path}")


if __name__ == "__main__":
    main()
