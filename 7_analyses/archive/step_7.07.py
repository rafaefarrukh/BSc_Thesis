"""
step_7.7_medium_optimization.py

Medium optimization / synthesis analysis:

  1. Screens a small set of physiologically relevant medium conditions
     (aerobic / microaerobic / anaerobic x nitrogen-replete / nitrogen-
     limited, plus a carbon-excess "feeding" condition) and computes, for
     each, the growth-coupled PHB production (growth fixed at 90% of that
     condition's max growth, then PHB maximized — consistent with the
     scenario used in step_7.2).
  2. Computes a minimal medium (cobra's minimal_medium) able to sustain a
     target growth rate, to identify which exchange reactions are actually
     required.

Outputs:
    step_7_outputs/step_7.7_medium_screen_results.csv
    step_7_outputs/step_7.7_minimal_medium.csv
    step_7_outputs/step_7.7_medium_optimization_summary.png
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cobra.medium import minimal_medium

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
)

# --------------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------------
PHB_REACTION_ID = "PHBS_syn_1"
CARBON_EXCHANGE_ID = "EX_glc__D_e"
NITROGEN_EXCHANGE_ID = "EX_nh4_e"
OXYGEN_EXCHANGE_ID = "EX_o2_e"

GROWTH_COUPLED_FRACTION = 0.90     # matches step_7.2's growth-coupled scenario
MINIMAL_MEDIUM_GROWTH_TARGET_FRACTION = 0.50  # fraction of max growth to sustain

# Uptake bound magnitudes (mmol/gDW/h) used to define each medium condition
GLC_EXCESS = 20.0
GLC_LIMITED = 5.0
NH4_REPLETE = 20.0
NH4_LIMITED = 1.0
O2_AEROBIC = 20.0
O2_MICROAEROBIC = 2.0
O2_ANAEROBIC = 0.0

MEDIA_CONDITIONS = {
    "Aerobic, N-replete":        {"glc": GLC_EXCESS,  "nh4": NH4_REPLETE, "o2": O2_AEROBIC},
    "Aerobic, N-limited":        {"glc": GLC_EXCESS,  "nh4": NH4_LIMITED, "o2": O2_AEROBIC},
    "Microaerobic, N-replete":   {"glc": GLC_EXCESS,  "nh4": NH4_REPLETE, "o2": O2_MICROAEROBIC},
    "Anaerobic, N-replete":      {"glc": GLC_EXCESS,  "nh4": NH4_REPLETE, "o2": O2_ANAEROBIC},
    "Carbon-excess feed, N-limited": {"glc": GLC_EXCESS * 2, "nh4": NH4_LIMITED, "o2": O2_AEROBIC},
    "Carbon-limited, N-replete": {"glc": GLC_LIMITED, "nh4": NH4_REPLETE, "o2": O2_AEROBIC},
}


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
    return biomass, phb, carbon, nitrogen, oxygen


def screen_media(model, biomass, phb, carbon, nitrogen, oxygen):
    rows = []
    for name, bounds in MEDIA_CONDITIONS.items():
        with model as m:
            m.reactions.get_by_id(carbon.id).lower_bound = -bounds["glc"]
            m.reactions.get_by_id(nitrogen.id).lower_bound = -bounds["nh4"]
            m.reactions.get_by_id(oxygen.id).lower_bound = -bounds["o2"]

            m.objective = biomass.id
            max_growth = m.slim_optimize()
            max_growth = max_growth if max_growth and max_growth > 1e-9 else 0.0

            phb_flux = 0.0
            fixed_growth = 0.0
            if max_growth > 0:
                fixed_growth = GROWTH_COUPLED_FRACTION * max_growth
                biomass_m = m.reactions.get_by_id(biomass.id)
                biomass_m.lower_bound = fixed_growth
                biomass_m.upper_bound = fixed_growth
                m.objective = phb.id
                phb_flux = m.slim_optimize()
                phb_flux = phb_flux if phb_flux and phb_flux > 1e-9 else 0.0

        rows.append({
            "medium": name,
            "glc_uptake_bound": bounds["glc"],
            "nh4_uptake_bound": bounds["nh4"],
            "o2_uptake_bound": bounds["o2"],
            "max_growth_rate_1_h": max_growth,
            "growth_coupled_growth_1_h": fixed_growth,
            "growth_coupled_PHB_flux": phb_flux,
        })
    return pd.DataFrame(rows)


def compute_minimal_medium(model, biomass, carbon, nitrogen, oxygen):
    with model as m:
        # start from an open, generous medium so minimal_medium has room to search
        m.reactions.get_by_id(carbon.id).lower_bound = -GLC_EXCESS
        m.reactions.get_by_id(nitrogen.id).lower_bound = -NH4_REPLETE
        m.reactions.get_by_id(oxygen.id).lower_bound = -O2_AEROBIC
        m.objective = biomass.id
        max_growth = m.slim_optimize()
        target = MINIMAL_MEDIUM_GROWTH_TARGET_FRACTION * max_growth
        medium_solution = minimal_medium(m, target, minimize_components=True)
    return medium_solution, target, max_growth


def main():
    ensure_output_dir(OUTPUT_DIR)
    model = load_model(MODEL_PATH)
    biomass, phb, carbon, nitrogen, oxygen = get_reactions(model)

    describe_reaction(biomass, "Biomass reaction")
    describe_reaction(phb, "PHB reaction")
    describe_reaction(carbon, "Carbon exchange")
    describe_reaction(nitrogen, "Nitrogen exchange")
    describe_reaction(oxygen, "Oxygen exchange")

    print("\nScreening medium conditions...")
    screen_df = screen_media(model, biomass, phb, carbon, nitrogen, oxygen)
    screen_csv = os.path.join(OUTPUT_DIR, "step_7.7_medium_screen_results.csv")
    screen_df.to_csv(screen_csv, index=False)
    print(screen_df.to_string(index=False))
    print(f"\nMedium screen written to: {screen_csv}")

    print("\nComputing minimal medium...")
    medium_solution, target_growth, max_growth = compute_minimal_medium(
        model, biomass, carbon, nitrogen, oxygen
    )
    minimal_csv = os.path.join(OUTPUT_DIR, "step_7.7_minimal_medium.csv")
    if medium_solution is not None:
        medium_solution.rename("uptake_flux").to_csv(minimal_csv, header=True)
        print(f"Minimal medium sustaining {target_growth:.4f} /h "
              f"({MINIMAL_MEDIUM_GROWTH_TARGET_FRACTION * 100:.0f}% of max growth "
              f"{max_growth:.4f} /h) requires {len(medium_solution)} exchange reactions.")
    else:
        print("No feasible minimal medium found for the requested growth target.")
    print(f"Minimal medium written to: {minimal_csv}")

    # --- Summary figure ---------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    x = np.arange(len(screen_df))
    width = 0.35
    ax1 = axes[0]
    ax1.bar(x - width / 2, screen_df["max_growth_rate_1_h"], width, label="Max growth rate", color="tab:blue")
    ax1.bar(x + width / 2, screen_df["growth_coupled_PHB_flux"], width, label="PHB flux (growth-coupled)", color="tab:green")
    ax1.set_xticks(x)
    ax1.set_xticklabels(screen_df["medium"], rotation=35, ha="right", fontsize=8)
    ax1.set_ylabel("Rate (1/h or mmol/gDW/h)")
    ax1.set_title("Growth and PHB Flux Across Candidate Media")
    ax1.legend()
    ax1.grid(alpha=0.3, axis="y")

    ax2 = axes[1]
    if medium_solution is not None and len(medium_solution) > 0:
        mm_sorted = medium_solution.sort_values(ascending=False)
        ax2.barh(mm_sorted.index, mm_sorted.values, color="tab:purple")
        ax2.set_xlabel("Required uptake flux (mmol/gDW/h)")
        ax2.set_title(
            f"Minimal Medium\n(sustaining {target_growth:.3f} /h, "
            f"{MINIMAL_MEDIUM_GROWTH_TARGET_FRACTION * 100:.0f}% of max growth)"
        )
        ax2.grid(alpha=0.3, axis="x")
    else:
        ax2.text(0.5, 0.5, "No feasible minimal medium found", ha="center", va="center")
        ax2.axis("off")

    fig.tight_layout()
    png_path = os.path.join(OUTPUT_DIR, "step_7.7_medium_optimization_summary.png")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    print(f"\nSummary figure written to: {png_path}")


if __name__ == "__main__":
    main()
