"""
step_7.24_global_sensitivity_sobol.py

Global Sensitivity Analysis (Sobol method): unlike the earlier one-factor-
at-a-time sweeps (step_7.8) and pairwise grids (step_7.9), Sobol analysis
varies ALL uncertain inputs simultaneously across their full ranges and
decomposes the variance of the model output into contributions from each
input individually (first-order index S1), plus each input's total
contribution including interactions with other inputs (total-order index
ST). A large gap between ST and S1 for an input means it matters mostly
through interactions with other inputs, not on its own.

Inputs varied here (uptake bound magnitudes, mmol/gDW/h):
    - glucose uptake
    - O2 uptake
    - NH4 uptake
Output analyzed: growth-coupled PHB flux (growth fixed at
GROWTH_COUPLED_FRACTION of its own max at each sampled point, then PHB
maximized) — same definition used throughout step_7.2 onward.

Requires: SALib (pip install SALib)

Outputs:
    step_7_outputs/step_7.24_sobol_samples_and_output.csv
    step_7_outputs/step_7.24_sobol_indices.csv
    step_7_outputs/step_7.24_sobol_indices_bar_chart.png
"""

import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from SALib.sample import saltelli
    from SALib.analyze import sobol
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "This script requires SALib. Install it with: pip install SALib"
    ) from exc

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

PHB_REACTION_ID = "PHBS_syn_1"
CARBON_EXCHANGE_ID = "EX_glc__D_e"
NITROGEN_EXCHANGE_ID = "EX_nh4_e"
OXYGEN_EXCHANGE_ID = "EX_o2_e"

GROWTH_COUPLED_FRACTION = 0.90

# Sobol base sample size N -> total model evaluations = N * (2D + 2) for D
# parameters (here D=3 -> 8N evaluations). Keep modest; this is the slow
# part of the script since each evaluation is 1-2 LP solves.
N_BASE_SAMPLES = 128

PROBLEM = {
    "num_vars": 3,
    "names": ["glucose_uptake", "o2_uptake", "nh4_uptake"],
    "bounds": [[0.0, 20.0], [0.0, 25.0], [0.0, 10.0]],
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


def evaluate(model, biomass, phb, carbon, nitrogen, oxygen, glc, o2, nh4):
    with model as m:
        m.reactions.get_by_id(carbon.id).lower_bound = -glc
        m.reactions.get_by_id(oxygen.id).lower_bound = -o2
        m.reactions.get_by_id(nitrogen.id).lower_bound = -nh4

        m.objective = biomass.id
        max_growth = m.slim_optimize()
        max_growth = max_growth if max_growth and max_growth > 1e-9 else 0.0
        if max_growth <= 0:
            return 0.0

        fixed_growth = GROWTH_COUPLED_FRACTION * max_growth
        biomass_m = m.reactions.get_by_id(biomass.id)
        biomass_m.lower_bound = fixed_growth
        biomass_m.upper_bound = fixed_growth
        m.objective = phb.id
        phb_flux = m.slim_optimize()
    return phb_flux if phb_flux and phb_flux > 1e-9 else 0.0


def main():
    ensure_output_dir(OUTPUT_DIR)
    model = load_model(MODEL_PATH)
    biomass, phb, carbon, nitrogen, oxygen = get_reactions(model)

    describe_reaction(biomass, "Biomass reaction")
    describe_reaction(phb, "PHB reaction")
    describe_reaction(carbon, "Carbon exchange")
    describe_reaction(oxygen, "Oxygen exchange")
    describe_reaction(nitrogen, "Nitrogen exchange")

    param_values = saltelli.sample(PROBLEM, N_BASE_SAMPLES)
    n_evals = param_values.shape[0]
    print(f"\nGenerated {n_evals} Sobol sample points (N_BASE_SAMPLES={N_BASE_SAMPLES}, "
          f"{PROBLEM['num_vars']} parameters). Evaluating model (slow step)...")

    outputs = np.zeros(n_evals)
    for i, (glc, o2, nh4) in enumerate(param_values):
        outputs[i] = evaluate(model, biomass, phb, carbon, nitrogen, oxygen, glc, o2, nh4)
        if (i + 1) % 200 == 0:
            print(f"  ...{i + 1}/{n_evals} evaluations complete")

    samples_df = pd.DataFrame(param_values, columns=PROBLEM["names"])
    samples_df["PHB_flux_output"] = outputs
    samples_csv = os.path.join(OUTPUT_DIR, "step_7.24_sobol_samples_and_output.csv")
    samples_df.to_csv(samples_csv, index=False)
    print(f"\nSample points and outputs written to: {samples_csv}")

    print("Computing Sobol sensitivity indices...")
    sobol_indices = sobol.analyze(PROBLEM, outputs, print_to_console=False)

    indices_df = pd.DataFrame({
        "parameter": PROBLEM["names"],
        "S1": sobol_indices["S1"],
        "S1_conf": sobol_indices["S1_conf"],
        "ST": sobol_indices["ST"],
        "ST_conf": sobol_indices["ST_conf"],
    })
    indices_csv = os.path.join(OUTPUT_DIR, "step_7.24_sobol_indices.csv")
    indices_df.to_csv(indices_csv, index=False)
    print(f"\nSobol indices:\n{indices_df.to_string(index=False)}")
    print(f"\nSobol indices written to: {indices_csv}")

    # Second-order interaction indices, if available
    if "S2" in sobol_indices:
        s2 = sobol_indices["S2"]
        interactions = []
        names = PROBLEM["names"]
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                interactions.append({"pair": f"{names[i]} x {names[j]}", "S2": s2[i, j]})
        interactions_df = pd.DataFrame(interactions)
        print(f"\nSecond-order interaction indices:\n{interactions_df.to_string(index=False)}")

    # --- Bar chart -----------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(PROBLEM["names"]))
    width = 0.35
    ax.bar(x - width / 2, indices_df["S1"], width, yerr=indices_df["S1_conf"],
           label="First-order (S1)", color="tab:blue", capsize=4)
    ax.bar(x + width / 2, indices_df["ST"], width, yerr=indices_df["ST_conf"],
           label="Total-order (ST)", color="tab:orange", capsize=4)
    ax.set_xticks(x)
    ax.set_xticklabels(indices_df["parameter"])
    ax.set_ylabel("Sobol sensitivity index")
    ax.set_title("Global Sensitivity Analysis (Sobol) — PHB Flux Output")
    ax.legend()
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()

    png_path = os.path.join(OUTPUT_DIR, "step_7.24_sobol_indices_bar_chart.png")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    print(f"\nBar chart written to: {png_path}")


if __name__ == "__main__":
    main()
