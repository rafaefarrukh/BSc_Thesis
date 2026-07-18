"""
step_7.6_dynamic_fba.py

Dynamic FBA (dFBA) using the static optimization approach (SOA;
Mahadevan et al., 2002): at each small time step, extracellular
substrate concentrations set kinetic (Michaelis-Menten) uptake bounds
on the model, FBA/pFBA is solved to get instantaneous rates, and the
system state (biomass, substrates, product) is advanced by explicit
Euler integration.

Three scenarios are simulated (PHB accumulation in many bacteria is
triggered by nitrogen starvation while carbon remains available, so
these are chosen to probe exactly that behavior):

  1. "standard"   — batch growth, nitrogen never limiting. Biomass
                     objective throughout. Expect growth until carbon
                     (glucose) is exhausted; PHB stays low.
  2. "n_limited"   — nitrogen supplied at a low initial concentration,
                     so it runs out mid-simulation while glucose is
                     still available. Biomass objective throughout —
                     growth naturally stalls once nitrogen is gone, and
                     any PHB pathway flux the model can still support
                     is exposed because leftover carbon has nowhere
                     else to go once growth is nitrogen-constrained.
  3. "two_stage"   — mimics a fed-batch PHB production process: Stage 1
                     is a nitrogen-replete growth phase (biomass
                     objective); at a fixed switch time the nitrogen
                     feed is cut and the objective is switched to
                     maximize PHB flux (subject to a minimum
                     maintenance-level growth constraint), modeling the
                     industrial "growth phase then production phase"
                     strategy.

Outputs:
    step_7_outputs/step_7.6_dfba_<scenario>.csv   (one file per scenario)
    step_7_outputs/step_7.6_dfba_timecourse.png   (combined comparison figure)
"""

import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cobra.flux_analysis import pfba

from step_7_utils import (
    MODEL_PATH,
    OUTPUT_DIR,
    ensure_output_dir,
    find_biomass_reaction,
    find_phb_reaction,
    find_first_match,
    load_model,
    michaelis_menten,
    GLUCOSE_EXCHANGE_CANDIDATES,
    OXYGEN_EXCHANGE_CANDIDATES,
    NITROGEN_EXCHANGE_CANDIDATES,
)

# --------------------------------------------------------------------------
# CONFIG — override any of these if auto-detection picks the wrong reaction
# --------------------------------------------------------------------------
PHB_REACTION_ID = "PHBS_syn_1"
CARBON_EXCHANGE_ID = "EX_glc__D_e"
NITROGEN_EXCHANGE_ID = "EX_nh4_e"
OXYGEN_EXCHANGE_ID = "EX_o2_e"

# Simulation settings
SIM_TIME_H = 24.0     # total simulated time, hours
DT_H = 0.1            # time step, hours

# Initial state (batch flask, arbitrary but reasonable starting point —
# adjust to match real experimental conditions if available)
X0 = 0.05      # initial biomass, gDW/L
GLC0 = 20.0    # initial glucose, mM
NH4_STANDARD0 = 20.0   # initial ammonium, mM (nitrogen-replete)
NH4_LIMITED0 = 1.0     # initial ammonium, mM (nitrogen-limited scenarios)

# Uptake kinetics (Michaelis-Menten); tune Vmax/Km to organism-specific
# values if known, these are generic placeholder values
GLC_VMAX = 10.0   # mmol/gDW/h
GLC_KM = 0.5      # mM
NH4_VMAX = 5.0    # mmol/gDW/h
NH4_KM = 0.1      # mM
O2_UPTAKE_BOUND = 20.0  # mmol/gDW/h, constant aerobic supply (not depleted)

# Two-stage scenario: time at which nitrogen feed is cut and objective
# switches from growth to PHB maximization
STAGE_SWITCH_TIME_H = 8.0
# Minimum growth (maintenance) enforced during the PHB-maximizing stage,
# as a fraction of the growth rate achieved right before the switch
STAGE2_MIN_GROWTH_FRACTION = 0.05


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


def run_dfba(model, biomass, phb, carbon, nitrogen, oxygen, scenario):
    """Run one dFBA scenario. Returns a pandas DataFrame time course."""
    n_steps = int(SIM_TIME_H / DT_H)

    X = X0
    glc = GLC0
    nh4 = NH4_LIMITED0 if scenario in ("n_limited", "two_stage") else NH4_STANDARD0
    phb_conc = 0.0

    records = []
    stage2_active = False
    stage2_min_growth = 0.0

    for step in range(n_steps + 1):
        t = step * DT_H
        records.append({
            "time_h": t, "biomass_gL": X, "glucose_mM": glc,
            "ammonium_mM": nh4, "PHB_mM": phb_conc,
        })
        if step == n_steps:
            break

        with model as m:
            # --- kinetic uptake bounds from current concentrations ---
            v_glc_kinetic = michaelis_menten(glc, GLC_VMAX, GLC_KM)
            v_nh4_kinetic = michaelis_menten(nh4, NH4_VMAX, NH4_KM)

            # cap so a substrate can't be driven below zero within this step
            # (classic dFBA numerical-stability correction)
            v_glc_cap = glc / (X * DT_H) if X > 0 else v_glc_kinetic
            v_nh4_cap = nh4 / (X * DT_H) if X > 0 else v_nh4_kinetic
            v_glc = min(v_glc_kinetic, v_glc_cap)
            v_nh4 = min(v_nh4_kinetic, v_nh4_cap)

            carbon_rxn_m = m.reactions.get_by_id(carbon.id)
            nitrogen_rxn_m = m.reactions.get_by_id(nitrogen.id)
            oxygen_rxn_m = m.reactions.get_by_id(oxygen.id)
            biomass_rxn_m = m.reactions.get_by_id(biomass.id)
            phb_rxn_m = m.reactions.get_by_id(phb.id)

            carbon_rxn_m.lower_bound = -v_glc
            nitrogen_rxn_m.lower_bound = -v_nh4
            oxygen_rxn_m.lower_bound = -O2_UPTAKE_BOUND

            # --- scenario-specific objective handling ---
            if scenario == "two_stage" and t >= STAGE_SWITCH_TIME_H:
                if not stage2_active:
                    stage2_active = True
                # enforce a maintenance-level minimum growth, then maximize PHB
                biomass_rxn_m.lower_bound = stage2_min_growth
                m.objective = phb_rxn_m.id
            else:
                m.objective = biomass_rxn_m.id

            try:
                sol = pfba(m)
                status_ok = True
            except Exception:
                status_ok = False

        if not status_ok or sol.status != "optimal":
            # infeasible step (e.g. all substrates exhausted) -> no more change
            growth_rate = 0.0
            v_glc_actual = 0.0
            v_phb = 0.0
        else:
            growth_rate = sol.fluxes[biomass.id]
            v_glc_actual = -sol.fluxes[carbon.id]  # positive = uptake
            v_nh4_actual = -sol.fluxes[nitrogen.id]
            v_phb = sol.fluxes[phb.id]
            if scenario == "two_stage" and t < STAGE_SWITCH_TIME_H:
                stage2_min_growth = growth_rate * STAGE2_MIN_GROWTH_FRACTION

        if not status_ok or sol.status != "optimal":
            v_nh4_actual = 0.0

        # --- Euler update ---
        dX = growth_rate * X * DT_H
        dGlc = -v_glc_actual * X * DT_H if status_ok else 0.0
        dNH4 = -v_nh4_actual * X * DT_H if status_ok else 0.0
        dPHB = v_phb * X * DT_H if status_ok else 0.0

        X = max(X + dX, 0.0)
        glc = max(glc + dGlc, 0.0)
        nh4 = max(nh4 + dNH4, 0.0)
        phb_conc = max(phb_conc + dPHB, 0.0)

    return pd.DataFrame.from_records(records)


def main():
    ensure_output_dir(OUTPUT_DIR)
    model = load_model(MODEL_PATH)
    biomass, phb, carbon, nitrogen, oxygen = get_reactions(model)

    print(f"Biomass:  {biomass.id}")
    print(f"PHB:      {phb.id}")
    print(f"Carbon:   {carbon.id}")
    print(f"Nitrogen: {nitrogen.id}")
    print(f"Oxygen:   {oxygen.id}")

    scenarios = ["standard", "n_limited", "two_stage"]
    results = {}
    for scenario in scenarios:
        print(f"\nRunning dFBA scenario: {scenario} ...")
        df = run_dfba(model, biomass, phb, carbon, nitrogen, oxygen, scenario)
        results[scenario] = df
        csv_path = os.path.join(OUTPUT_DIR, f"step_7.6_dfba_{scenario}.csv")
        df.to_csv(csv_path, index=False)
        print(f"  final biomass={df['biomass_gL'].iloc[-1]:.3f} g/L, "
              f"final PHB={df['PHB_mM'].iloc[-1]:.3f} mM  -> {csv_path}")

    # --- Combined time-course figure -----------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    labels = {"standard": "Standard", "n_limited": "N-limited", "two_stage": "Two-stage"}
    colors = {"standard": "tab:blue", "n_limited": "tab:orange", "two_stage": "tab:green"}

    panels = [
        ("biomass_gL", "Biomass (g/L)", axes[0, 0]),
        ("glucose_mM", "Glucose (mM)", axes[0, 1]),
        ("ammonium_mM", "Ammonium (mM)", axes[1, 0]),
        ("PHB_mM", "PHB (mM)", axes[1, 1]),
    ]
    for col, ylabel, ax in panels:
        for scenario in scenarios:
            df = results[scenario]
            ax.plot(df["time_h"], df[col], label=labels[scenario], color=colors[scenario], lw=2)
        ax.set_ylabel(ylabel)
        ax.set_xlabel("Time (h)")
        ax.grid(alpha=0.3)

    axes[1, 0].axvline(STAGE_SWITCH_TIME_H, color="gray", linestyle="--", lw=1,
                        label="Two-stage switch")
    axes[0, 0].legend(loc="best", fontsize=9)
    fig.suptitle("Dynamic FBA — Standard vs. N-limited vs. Two-stage PHB Production", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    png_path = os.path.join(OUTPUT_DIR, "step_7.6_dfba_timecourse.png")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    print(f"\nCombined time-course figure written to: {png_path}")


if __name__ == "__main__":
    main()
