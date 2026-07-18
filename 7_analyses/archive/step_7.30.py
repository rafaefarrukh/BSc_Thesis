"""
STEP 7.30 — CORRECTED DYNAMIC FBA (MASS-BALANCE FIX)

The original step_7.6 "standard" and "N-limited" dFBA runs showed biomass
climbing to ~1740 g/L and glucose staying pinned at exactly 20 mM for the
entire simulation. Both are signs that the exchange-flux -> concentration
update wasn't actually coupled: uptake fluxes (mmol/gDW/h) must be
multiplied by biomass (gDW/L) AND the timestep (h) to get a concentration
change (mM), and growth-rate integration needs a small enough timestep
(with substep re-solving) to stay numerically stable and to correctly
deplete substrates.

Fixes applied here:
  1. Explicit unit-consistent Euler update:
       dX/dt = mu * X                        (gDW/L)
       dS/dt = -q_S * X   for each substrate  (mM, since q_S in mmol/gDW/h)
  2. Uptake bounds are capped each step so a substrate cannot be consumed
     below zero within that timestep (bound = min(default_bound,
     S / (X * dt)) ), which is what actually produces depletion kinetics
     instead of a flat concentration line.
  3. A small, fixed timestep (dt = 0.05 h) with re-solved FBA every step.
  4. A sanity cap: if biomass exceeds a configurable ceiling (default
     100 gDW/L, well above typical shake-flask/bioreactor titers), the
     simulation stops early and prints a warning rather than silently
     producing an unphysical curve.

Input:  ../models/model_7.xml
Output: ./step_7.30_dfba_standard_corrected.csv
        ./step_7.30_dfba_n_limited_corrected.csv
        ./step_7.30_dfba_two_stage_corrected.csv
        ./step_7.30_dfba_timecourse_corrected.png
        ./step_7.30_summary.txt
"""

import cobra
import pandas as pd
import matplotlib.pyplot as plt

MODEL_PATH = "../models/model_7.xml"
BIOMASS_ID = "Growth"
PHB_ID = "PHBS_syn_1"
GLC_EX = "EX_glc__D_e"
NH4_EX = "EX_nh4_e"

DT = 0.05          # h, integration timestep
T_MAX = 24.0        # h, total simulated time
X0 = 0.05           # gDW/L, initial biomass
S0_GLC = 20.0        # mM, initial glucose
S0_NH4 = 20.0        # mM, initial ammonium (N-replete start)
S0_NH4_LIMITED = 1.0  # mM, initial ammonium (N-limited start)
BIOMASS_SANITY_CAP = 100.0  # gDW/L; stop and warn if exceeded

print(f"Loading model from: {MODEL_PATH}")
model = cobra.io.read_sbml_model(MODEL_PATH)
print(f"Loaded model '{model.id}': {len(model.reactions)} reactions, "
      f"{len(model.metabolites)} metabolites, {len(model.genes)} genes")

glc_default_lb = model.reactions.get_by_id(GLC_EX).lower_bound
nh4_default_lb = model.reactions.get_by_id(NH4_EX).lower_bound
print(f"Default uptake caps: glucose lb={glc_default_lb}, NH4 lb={nh4_default_lb}")


def run_dfba(scenario_name, s0_glc, s0_nh4, two_stage_switch_h=None):
    """
    two_stage_switch_h: if set, at that simulated time the biomass reaction
    is bounded to near-zero (mimicking an N-starvation-triggered production
    phase) so carbon is redirected to PHB, matching the intent of the
    original 'two_stage' scenario.
    """
    X, glc, nh4, phb = X0, s0_glc, s0_nh4, 0.0
    t = 0.0
    rows = [dict(time_h=t, biomass_gL=X, glucose_mM=glc, ammonium_mM=nh4, PHB_mM=phb)]
    capped = False

    while t < T_MAX:
        with model:
            # Cap uptake so we can't draw a substrate below zero this step
            glc_cap = min(abs(glc_default_lb), glc / (X * DT)) if (X > 0 and glc > 0) else 0.0
            nh4_cap = min(abs(nh4_default_lb), nh4 / (X * DT)) if (X > 0 and nh4 > 0) else 0.0
            model.reactions.get_by_id(GLC_EX).lower_bound = -glc_cap
            model.reactions.get_by_id(NH4_EX).lower_bound = -nh4_cap

            if two_stage_switch_h is not None and t >= two_stage_switch_h:
                # Production phase: cap growth hard, let PHB objective dominate
                model.reactions.get_by_id(BIOMASS_ID).upper_bound = 0.01
                model.objective = PHB_ID
            else:
                model.objective = BIOMASS_ID

            sol = model.optimize()
            if sol.status != "optimal":
                print(f"  [{scenario_name}] infeasible at t={t:.2f}h — stopping early.")
                break

            mu = sol.fluxes[BIOMASS_ID]
            q_glc = sol.fluxes[GLC_EX]     # negative = uptake
            q_nh4 = sol.fluxes[NH4_EX]
            q_phb = sol.fluxes[PHB_ID]

        # Euler update, unit-consistent (mmol/gDW/h * gDW/L * h = mmol/L = mM)
        X_new = max(X + mu * X * DT, 0.0)
        glc_new = max(glc + q_glc * X * DT, 0.0)
        nh4_new = max(nh4 + q_nh4 * X * DT, 0.0)
        phb_new = max(phb + q_phb * X * DT, 0.0)

        X, glc, nh4, phb = X_new, glc_new, nh4_new, phb_new
        t += DT
        rows.append(dict(time_h=round(t, 3), biomass_gL=X, glucose_mM=glc,
                          ammonium_mM=nh4, PHB_mM=phb))

        if X > BIOMASS_SANITY_CAP:
            print(f"  [{scenario_name}] WARNING: biomass exceeded sanity cap "
                  f"({BIOMASS_SANITY_CAP} gDW/L) at t={t:.2f}h — stopping early. "
                  "This usually means the growth rate is unrealistically high "
                  "for the timestep; consider revisiting biomass/maintenance terms.")
            capped = True
            break

    df = pd.DataFrame(rows)
    print(f"[{scenario_name}] final biomass={X:.3f} g/L, final glucose={glc:.3f} mM, "
          f"final PHB={phb:.3f} mM" + ("  (STOPPED EARLY — see warning)" if capped else ""))
    return df


print("=" * 70)
print("STEP 7.30 — CORRECTED DYNAMIC FBA")
print("=" * 70)

print("\nRunning corrected dFBA scenario: standard ...")
df_standard = run_dfba("standard", S0_GLC, S0_NH4)
df_standard.to_csv("./step_7.30_dfba_standard_corrected.csv", index=False)

print("\nRunning corrected dFBA scenario: n_limited ...")
df_nlim = run_dfba("n_limited", S0_GLC, S0_NH4_LIMITED)
df_nlim.to_csv("./step_7.30_dfba_n_limited_corrected.csv", index=False)

print("\nRunning corrected dFBA scenario: two_stage ...")
df_two = run_dfba("two_stage", S0_GLC, S0_NH4, two_stage_switch_h=8.0)
df_two.to_csv("./step_7.30_dfba_two_stage_corrected.csv", index=False)

# ---- Combined figure ----
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
scenarios = {"Standard": df_standard, "N-limited": df_nlim, "Two-stage": df_two}
for name, df in scenarios.items():
    axes[0, 0].plot(df["time_h"], df["biomass_gL"], label=name)
    axes[0, 1].plot(df["time_h"], df["glucose_mM"], label=name)
    axes[1, 0].plot(df["time_h"], df["ammonium_mM"], label=name)
    axes[1, 1].plot(df["time_h"], df["PHB_mM"], label=name)
axes[0, 0].set_title("Biomass (g/L)"); axes[0, 0].set_xlabel("Time (h)")
axes[0, 1].set_title("Glucose (mM)"); axes[0, 1].set_xlabel("Time (h)")
axes[1, 0].set_title("Ammonium (mM)"); axes[1, 0].set_xlabel("Time (h)")
axes[1, 1].set_title("PHB (mM)"); axes[1, 1].set_xlabel("Time (h)")
for ax in axes.flat:
    ax.legend()
fig.suptitle("Corrected Dynamic FBA — Standard vs. N-limited vs. Two-stage")
fig.tight_layout()
fig.savefig("./step_7.30_dfba_timecourse_corrected.png", dpi=150)
print("\nCombined figure written to: ./step_7.30_dfba_timecourse_corrected.png")

with open("./step_7.30_summary.txt", "w") as f:
    f.write("STEP 7.30 — Corrected dFBA summary\n")
    f.write("=" * 60 + "\n")
    f.write("Fixes applied vs. step_7.6:\n")
    f.write("  1. Uptake bounds capped per-step by available substrate/(X*dt),\n")
    f.write("     which now correctly depletes glucose/ammonium over time.\n")
    f.write("  2. Explicit unit-consistent Euler update for biomass and each\n")
    f.write("     substrate/product pool.\n")
    f.write(f"  3. Small fixed timestep (dt={DT} h) with re-solved FBA each step.\n")
    f.write(f"  4. Sanity cap at {BIOMASS_SANITY_CAP} gDW/L to catch runaway growth\n")
    f.write("     early rather than silently reporting unphysical titers.\n\n")
    for name, df in scenarios.items():
        f.write(f"{name}: final biomass={df['biomass_gL'].iloc[-1]:.3f} g/L, "
                f"final glucose={df['glucose_mM'].iloc[-1]:.3f} mM, "
                f"final PHB={df['PHB_mM'].iloc[-1]:.3f} mM\n")

print("\nSummary written to: ./step_7.30_summary.txt")
