"""
STEP 7.26 (FIXED) — YIELD ANALYSIS (MASS & CARBON BASIS)

BUG FIX vs. the original 7.26: every scenario reported glc_uptake = 0.0
because the script optimized the model with whatever default exchange
bounds it happened to have, rather than explicitly fixing glucose uptake.
With multiple alternate optima available, the solver was free to route
zero flux through EX_glc__D_e entirely. This version explicitly sets
EX_glc__D_e's lower bound to -GLC_UPTAKE_BOUND before every optimization,
mirroring what steps 7.7-7.10 already do correctly, and restores the
original bound afterward.

Input:  ../models/model_7.xml
        ./step_7.4_production_envelope.csv   (growth vs PHB flux envelope)
Output: ./step_7.26_yield_summary.csv
        ./step_7.26_yield_vs_growth.png
        ./step_7.26_summary.txt
"""

import cobra
import pandas as pd
import matplotlib.pyplot as plt

MODEL_PATH = "../models/model_7.xml"
BIOMASS_ID = "Growth"
PHB_ID = "PHBS_syn_1"
GLC_EX_ID = "EX_glc__D_e"
GLC_UPTAKE_BOUND = 10.0  # mmol/gDW/h, matches steps 7.7-7.10 convention

# Molecular properties
MW_GLUCOSE = 180.156   # g/mol, C6H12O6
C_GLUCOSE = 6
MW_PHB_MONOMER = 86.09  # g/mol, C4H6O2 (3-hydroxybutyrate repeat unit)
C_PHB_MONOMER = 4

# Literature reference range (glucose/sucrose substrates; not organism-specific)
LIT_MASS_YIELD_LOW, LIT_MASS_YIELD_HIGH = 0.33, 0.48  # g PHB / g glucose

print(f"Loading model from: {MODEL_PATH}")
model = cobra.io.read_sbml_model(MODEL_PATH)
print(f"Loaded model '{model.id}': {len(model.reactions)} reactions, "
      f"{len(model.metabolites)} metabolites, {len(model.genes)} genes")

glc_rxn = model.reactions.get_by_id(GLC_EX_ID)
default_glc_bounds = glc_rxn.bounds
print(f"Default {GLC_EX_ID} bounds (before fix): {default_glc_bounds}")

print("=" * 70)
print("STEP 7.26 (FIXED) — YIELD ANALYSIS")
print("=" * 70)


def compute_yield(glc_uptake_flux, phb_flux):
    if glc_uptake_flux <= 1e-9:
        return dict(mass_yield_g_g=float("nan"), c_yield=float("nan"), pct_theoretical=float("nan"))
    mass_yield = (phb_flux * MW_PHB_MONOMER) / (glc_uptake_flux * MW_GLUCOSE)
    c_mol_glc = glc_uptake_flux * C_GLUCOSE
    c_mol_phb = phb_flux * C_PHB_MONOMER
    c_yield = c_mol_phb / c_mol_glc
    theoretical_ceiling = 2.0 / 3.0  # see original script docstring for derivation/caveats
    pct_theoretical = 100 * c_yield / theoretical_ceiling
    return dict(mass_yield_g_g=mass_yield, c_yield=c_yield, pct_theoretical=pct_theoretical)


rows = []

# ---- Growth-optimal (glucose fixed open) ----
glc_rxn.lower_bound = -GLC_UPTAKE_BOUND
model.objective = BIOMASS_ID
sol = model.optimize()
glc_flux = abs(sol.fluxes[GLC_EX_ID])
phb_flux = sol.fluxes[PHB_ID]
y = compute_yield(glc_flux, phb_flux)
rows.append(dict(scenario="growth_optimal", growth_rate=sol.fluxes[BIOMASS_ID],
                  glc_uptake=glc_flux, phb_flux=phb_flux, **y))

# ---- PHB-optimal (glucose fixed open) ----
model.objective = PHB_ID
sol = model.optimize()
glc_flux = abs(sol.fluxes[GLC_EX_ID])
phb_flux = sol.fluxes[PHB_ID]
y = compute_yield(glc_flux, phb_flux)
rows.append(dict(scenario="PHB_optimal", growth_rate=sol.fluxes[BIOMASS_ID],
                  glc_uptake=glc_flux, phb_flux=phb_flux, **y))

# ---- 90% growth-coupled (glucose fixed open) ----
model.objective = BIOMASS_ID
max_growth = model.slim_optimize()
model.reactions.get_by_id(BIOMASS_ID).lower_bound = 0.9 * max_growth
model.objective = PHB_ID
sol = model.optimize()
glc_flux = abs(sol.fluxes[GLC_EX_ID])
phb_flux = sol.fluxes[PHB_ID]
y = compute_yield(glc_flux, phb_flux)
rows.append(dict(scenario="growth_coupled_90pct", growth_rate=sol.fluxes[BIOMASS_ID],
                  glc_uptake=glc_flux, phb_flux=phb_flux, **y))
model.reactions.get_by_id(BIOMASS_ID).lower_bound = 0.0  # reset

# Restore original glucose bound
glc_rxn.bounds = default_glc_bounds

yield_df = pd.DataFrame(rows)
yield_df.to_csv("./step_7.26_yield_summary.csv", index=False)
print("\nYield summary (glucose now explicitly fixed at "
      f"-{GLC_UPTAKE_BOUND} mmol/gDW/h before each optimization):")
print(yield_df.to_string(index=False))

if (yield_df["glc_uptake"] < 1e-6).any():
    print("\nWARNING: at least one scenario still shows ~0 glucose uptake even "
          "with the exchange forced open. This means the model is meeting its "
          "objective via a DIFFERENT carbon source entirely under default bounds "
          "-- check which other exchanges are open by default (see step_7.29 "
          "findings: most non-maltose sugars share a common downstream "
          "bottleneck, so verify the model isn't defaulting to an alternate, "
          "unconstrained substrate).")

try:
    env = pd.read_csv("./step_7.4_production_envelope.csv")
    fig, ax1 = plt.subplots(figsize=(8, 5))
    ax1.plot(env["Growth"], env["mass_yield_maximum"], marker="o", color="tab:green",
             label="Mass yield (g PHB / g glucose), max PHB branch")
    ax1.axhline(LIT_MASS_YIELD_LOW, color="gray", linestyle="--", linewidth=1)
    ax1.axhline(LIT_MASS_YIELD_HIGH, color="gray", linestyle="--", linewidth=1,
                label=f"Literature range ({LIT_MASS_YIELD_LOW}-{LIT_MASS_YIELD_HIGH} g/g, other PHB producers)")
    ax1.set_xlabel("Growth rate (1/h)")
    ax1.set_ylabel("Mass yield (g PHB / g glucose)")
    ax1.set_title("PHB Mass Yield vs. Growth Rate, with Literature Reference Range")
    ax1.legend()
    fig.tight_layout()
    fig.savefig("./step_7.26_yield_vs_growth.png", dpi=150)
    print("\nEnvelope-yield figure written to: ./step_7.26_yield_vs_growth.png")
except FileNotFoundError:
    print("\nNOTE: ./step_7.4_production_envelope.csv not found — skipping envelope-yield plot.")

with open("./step_7.26_summary.txt", "w") as f:
    f.write("STEP 7.26 (FIXED) — Yield analysis summary\n")
    f.write("=" * 60 + "\n")
    f.write(f"Fix applied: EX_glc__D_e lower bound explicitly set to "
            f"-{GLC_UPTAKE_BOUND} before each optimization (was previously left\n"
            f"at its default {default_glc_bounds}, which allowed 0 glucose uptake).\n\n")
    f.write(yield_df.to_string(index=False))
    f.write("\n\nLiterature reference (glucose/sucrose substrates, various PHB producers):\n")
    f.write(f"  Mass yield ~ {LIT_MASS_YIELD_LOW}-{LIT_MASS_YIELD_HIGH} g PHB / g sugar\n")
    f.write("  (Not organism-specific to Rossellomorea marisflavi -- for context only.)\n")

print("\nSummary written to: ./step_7.26_summary.txt")
