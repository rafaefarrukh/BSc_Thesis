"""
step_7.23_metabolic_control_analysis.py

Metabolic Control Analysis (MCA) — flux control coefficients.

Classical MCA defines the flux control coefficient of enzyme/reaction i
on flux J as:
    C_i^J = (dJ / J) / (dv_i / v_i) = (v_i / J) * (dJ / dv_i)
i.e. the fractional change in a target flux J per fractional change in
reaction i's capacity, evaluated at steady state. It normally requires a
kinetic model (enzyme elasticities). Since this is a purely stoichiometric
(FBA) genome-scale model, this script computes an FBA-based numerical
analogue: it perturbs each reaction's flux bound by a small percentage,
re-optimizes the PHB objective (growth-coupled), and takes the resulting
finite-difference sensitivity as an approximation to the control
coefficient. This is a standard, widely-used proxy in the metabolic
engineering literature when kinetic parameters aren't available.

Scope: reactions within CONTROL_MAX_DISTANCE hops of PHB synthesis
(reuses step_7.12's pathway-distance ranking) — flux control coefficients
are most interpretable, and cheapest to compute reliably, close to the
pathway of interest.

Requires step_7.12 to have been run first.

Outputs:
    step_7_outputs/step_7.23_flux_control_coefficients.csv
    step_7_outputs/step_7.23_control_coefficients_bar_chart.png
"""

import os

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
    find_phb_reaction,
    load_model,
)

PHB_REACTION_ID = "PHBS_syn_1"
CONTROL_MAX_DISTANCE = 3
GROWTH_COUPLED_FRACTION = 0.90
PERTURBATION_FRACTION = 0.01   # +/-1% perturbation of each reaction's bound magnitude
TOP_N_PLOT = 25

DISTANCE_CSV = os.path.join(OUTPUT_DIR, "step_7.12_phb_pathway_distances.csv")


def phb_at_growth_coupled_optimum(model, biomass, phb):
    with model as m:
        m.objective = biomass.id
        max_growth = m.slim_optimize()
        if not max_growth or max_growth <= 1e-9:
            return 0.0
        fixed_growth = GROWTH_COUPLED_FRACTION * max_growth
        biomass_m = m.reactions.get_by_id(biomass.id)
        biomass_m.lower_bound = fixed_growth
        biomass_m.upper_bound = fixed_growth
        m.objective = phb.id
        phb_flux = m.slim_optimize()
    return phb_flux if phb_flux else 0.0


def perturbed_bound_capacity(rxn, direction, fraction):
    """Return (new_lb, new_ub) after scaling the reaction's *capacity*
    (the bound with larger magnitude) up (+) or down (-) by `fraction`,
    preserving the sign/directionality of the original bounds."""
    lb, ub = rxn.lower_bound, rxn.upper_bound
    scale = 1 + direction * fraction
    new_lb = lb * scale if lb < 0 else lb
    new_ub = ub * scale if ub > 0 else ub
    return new_lb, new_ub


def main():
    ensure_output_dir(OUTPUT_DIR)
    model = load_model(MODEL_PATH)
    biomass = find_biomass_reaction(model)
    phb = model.reactions.get_by_id(PHB_REACTION_ID) if PHB_REACTION_ID else find_phb_reaction(model)
    if phb is None:
        raise ValueError("Could not auto-detect a PHB reaction. Set PHB_REACTION_ID at the top of this script.")

    describe_reaction(biomass, "Biomass reaction")
    describe_reaction(phb, "PHB reaction")

    if not os.path.exists(DISTANCE_CSV):
        raise SystemExit(
            f"Could not find {DISTANCE_CSV}. Run step_7.12_phb_pathway_distance_and_map.py "
            f"first — this script reuses its pathway-distance ranking."
        )
    dist_df = pd.read_csv(DISTANCE_CSV)
    scoped_ids = dist_df.loc[dist_df["distance_from_PHB"] <= CONTROL_MAX_DISTANCE, "reaction_id"].tolist()
    print(f"Computing flux control coefficients for {len(scoped_ids)} reactions within "
          f"{CONTROL_MAX_DISTANCE} hops of PHB synthesis...")

    J0 = phb_at_growth_coupled_optimum(model, biomass, phb)
    print(f"Baseline (growth-coupled) PHB flux J0 = {J0:.6f}")
    if J0 <= 1e-9:
        raise SystemExit("Baseline PHB flux is ~0; control coefficients are not meaningful here. "
                          "Check PHB_REACTION_ID / medium bounds.")

    rows = []
    for rid in scoped_ids:
        rxn = model.reactions.get_by_id(rid)
        v_i = rxn.upper_bound if rxn.upper_bound > 0 else rxn.lower_bound
        if abs(v_i) < 1e-9:
            continue  # reaction has no meaningful capacity to perturb

        with model as m:
            m_rxn = m.reactions.get_by_id(rid)
            new_lb, new_ub = perturbed_bound_capacity(m_rxn, +1, PERTURBATION_FRACTION)
            m_rxn.lower_bound, m_rxn.upper_bound = new_lb, new_ub
            J_plus = phb_at_growth_coupled_optimum(m, biomass, phb)

        with model as m:
            m_rxn = m.reactions.get_by_id(rid)
            new_lb, new_ub = perturbed_bound_capacity(m_rxn, -1, PERTURBATION_FRACTION)
            m_rxn.lower_bound, m_rxn.upper_bound = new_lb, new_ub
            J_minus = phb_at_growth_coupled_optimum(m, biomass, phb)

        dJ = (J_plus - J_minus) / 2.0
        d_capacity = PERTURBATION_FRACTION * abs(v_i)
        # scaled (dimensionless) control coefficient: (dJ/J0) / (d_capacity/v_i)
        control_coeff = (dJ / J0) / (d_capacity / abs(v_i)) if d_capacity > 0 else 0.0

        rows.append({
            "reaction_id": rid,
            "reaction_name": rxn.name,
            "distance_from_PHB": dist_df.set_index("reaction_id").loc[rid, "distance_from_PHB"]
            if rid in dist_df["reaction_id"].values else None,
            "J_plus": J_plus, "J_minus": J_minus,
            "flux_control_coefficient": control_coeff,
        })

    df = pd.DataFrame(rows).sort_values("flux_control_coefficient", key=abs, ascending=False)
    csv_path = os.path.join(OUTPUT_DIR, "step_7.23_flux_control_coefficients.csv")
    df.to_csv(csv_path, index=False)
    print(f"\nFlux control coefficient table written to: {csv_path}")

    top = df.head(TOP_N_PLOT)
    print(f"\nTop {len(top)} reactions by |flux control coefficient| on PHB flux:")
    print(top[["reaction_id", "distance_from_PHB", "flux_control_coefficient"]].to_string(index=False))

    fig, ax = plt.subplots(figsize=(9, max(5, 0.3 * len(top))))
    colors = ["tab:green" if v > 0 else "tab:red" for v in top["flux_control_coefficient"]]
    ax.barh(top["reaction_id"], top["flux_control_coefficient"], color=colors)
    ax.axvline(0, color="black", lw=0.8)
    ax.invert_yaxis()
    ax.set_xlabel("Flux control coefficient on PHB flux (FBA finite-difference approximation)")
    ax.set_title(f"Top {len(top)} Flux Control Coefficients (green = positive control, red = negative)")
    fig.tight_layout()

    png_path = os.path.join(OUTPUT_DIR, "step_7.23_control_coefficients_bar_chart.png")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    print(f"\nBar chart written to: {png_path}")


if __name__ == "__main__":
    main()
