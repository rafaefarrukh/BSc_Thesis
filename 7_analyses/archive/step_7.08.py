"""
step_7.8_1D_robustness.py

1D robustness (single-parameter sensitivity) analysis. Independently
sweeps each of the three main uptake exchanges — glucose (carbon),
O2 (aeration), NH4 (nitrogen) — across a range of uptake bounds while
holding the other two at generous, non-limiting levels. At each point
records:
    - maximum growth rate
    - growth-coupled PHB flux (growth fixed at 90% of that point's max
      growth, then PHB maximized — consistent with step_7.2/step_7.7)

This single sweep serves both the "1D robustness" and the "oxygen
sensitivity sweep / nitrogen-limitation sweep" analyses: the O2 and NH4
columns of the same sweep directly are the O2 and N sensitivity sweeps,
with critical thresholds (onset of limitation, saturation point)
identified automatically below.

Outputs:
    step_7_outputs/step_7.8_1D_robustness_glucose.csv
    step_7_outputs/step_7.8_1D_robustness_oxygen.csv
    step_7_outputs/step_7.8_1D_robustness_nitrogen.csv
    step_7_outputs/step_7.8_1D_robustness.png
    step_7_outputs/step_7.8_thresholds_summary.txt
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
MALTOSE_ID = MALTOSE_EXCHANGE_ID  # "EX_malt_e" — included as a 4th robustness axis

GROWTH_COUPLED_FRACTION = 0.90
N_POINTS = 30

# "Non-limiting" background level used for the nutrients NOT being swept
GLC_NONLIMITING = 20.0
O2_NONLIMITING = 20.0
NH4_NONLIMITING = 20.0
MALTOSE_NONLIMITING = 0.0  # maltose closed by default when it's not the swept nutrient

# Sweep ranges (mmol/gDW/h, uptake magnitude)
GLC_RANGE = (0.0, 20.0)
O2_RANGE = (0.0, 25.0)
NH4_RANGE = (0.0, 10.0)
MALTOSE_RANGE = (0.0, 20.0)

# Growth rate below this fraction of the sweep's own maximum is considered
# "no growth" for threshold-detection purposes
GROWTH_EPS_FRACTION = 0.01


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
              f"maltose robustness sweep will be skipped.")
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


def sweep_nutrient(model, biomass, phb, carbon, nitrogen, oxygen, maltose,
                    varied, value_range):
    """varied: one of 'glucose', 'oxygen', 'nitrogen', 'maltose'."""
    values = np.linspace(value_range[0], value_range[1], N_POINTS)
    rows = []
    for v in values:
        glc = v if varied == "glucose" else GLC_NONLIMITING
        o2 = v if varied == "oxygen" else O2_NONLIMITING
        nh4 = v if varied == "nitrogen" else NH4_NONLIMITING
        malt = v if varied == "maltose" else MALTOSE_NONLIMITING
        # When sweeping maltose specifically, isolate it as the sole carbon
        # source (close glucose) so the sweep reflects maltose's own
        # growth-supporting capacity rather than glucose masking it.
        glc_for_point = 0.0 if varied == "maltose" else glc
        growth, phb_flux = growth_and_phb_at_point(
            model, biomass, phb, carbon, nitrogen, oxygen, maltose,
            glc_for_point, o2, nh4, malt_bound=malt
        )
        rows.append({"uptake_bound": v, "growth_rate_1_h": growth, "PHB_flux": phb_flux})
    return pd.DataFrame(rows)


def detect_thresholds(df, label):
    """Identify (a) the uptake level at which growth first becomes
    essentially nonzero (limitation threshold) and (b) the uptake level
    beyond which growth stops increasing appreciably (saturation point)."""
    max_growth = df["growth_rate_1_h"].max()
    if max_growth <= 0:
        return f"{label}: no growth achieved anywhere in the swept range."

    eps = GROWTH_EPS_FRACTION * max_growth
    nonzero = df[df["growth_rate_1_h"] > eps]
    limitation_threshold = nonzero["uptake_bound"].iloc[0] if len(nonzero) else None

    # saturation: first point within 99% of the max growth rate
    saturated = df[df["growth_rate_1_h"] >= 0.99 * max_growth]
    saturation_point = saturated["uptake_bound"].iloc[0] if len(saturated) else None

    return (
        f"{label}: growth becomes nonzero at uptake >= "
        f"{limitation_threshold:.3f} mmol/gDW/h; reaches ~99% of max growth "
        f"({max_growth:.4f} /h) at uptake >= {saturation_point:.3f} mmol/gDW/h."
    )


def main():
    ensure_output_dir(OUTPUT_DIR)
    model = load_model(MODEL_PATH)
    biomass, phb, carbon, nitrogen, oxygen, maltose = get_reactions(model)

    describe_reaction(biomass, "Biomass reaction")
    describe_reaction(phb, "PHB reaction")
    describe_reaction(carbon, "Carbon exchange")
    describe_reaction(oxygen, "Oxygen exchange")
    describe_reaction(nitrogen, "Nitrogen exchange")
    if maltose is not None:
        describe_reaction(maltose, "Maltose exchange")

    sweep_specs = [
        ("Glucose", "glucose", GLC_RANGE),
        ("Oxygen", "oxygen", O2_RANGE),
        ("Nitrogen (NH4)", "nitrogen", NH4_RANGE),
    ]
    if maltose is not None:
        sweep_specs.append(("Maltose", "maltose", MALTOSE_RANGE))

    sweeps = {}
    threshold_lines = []
    for label, varied, value_range in sweep_specs:
        print(f"\nSweeping {label} uptake ({N_POINTS} points)...")
        df = sweep_nutrient(model, biomass, phb, carbon, nitrogen, oxygen, maltose, varied, value_range)
        sweeps[label] = df
        csv_path = os.path.join(OUTPUT_DIR, f"step_7.8_1D_robustness_{varied}.csv")
        df.to_csv(csv_path, index=False)
        print(f"  written to {csv_path}")
        threshold_lines.append(detect_thresholds(df, label))

    # --- Combined figure: 2 rows (growth, PHB) x N columns (nutrients) -------
    n_cols = len(sweep_specs)
    fig, axes = plt.subplots(2, n_cols, figsize=(14 * n_cols / 3, 7), sharex="col")
    if n_cols == 1:
        axes = axes.reshape(2, 1)
    for col, (label, df) in enumerate(sweeps.items()):
        axes[0, col].plot(df["uptake_bound"], df["growth_rate_1_h"], color="tab:blue", lw=2)
        axes[0, col].set_title(f"{label} uptake")
        axes[0, col].set_ylabel("Max growth rate (1/h)" if col == 0 else "")
        axes[0, col].grid(alpha=0.3)

        axes[1, col].plot(df["uptake_bound"], df["PHB_flux"], color="tab:green", lw=2)
        axes[1, col].set_xlabel("Uptake bound (mmol/gDW/h)")
        axes[1, col].set_ylabel("PHB flux, growth-coupled" if col == 0 else "")
        axes[1, col].grid(alpha=0.3)

    fig.suptitle("1D Robustness — Growth and PHB Flux vs. Single-Nutrient Uptake", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    png_path = os.path.join(OUTPUT_DIR, "step_7.8_1D_robustness.png")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    # --- Threshold summary --------------------------------------------------
    txt_path = os.path.join(OUTPUT_DIR, "step_7.8_thresholds_summary.txt")
    with open(txt_path, "w") as fh:
        fh.write("STEP 7.8 — 1D Robustness: limitation / saturation thresholds\n")
        fh.write("=" * 70 + "\n")
        fh.write("\n".join(threshold_lines))

    print("\n" + "\n".join(threshold_lines))
    print(f"\nFigure written to: {png_path}")
    print(f"Threshold summary written to: {txt_path}")


if __name__ == "__main__":
    main()
