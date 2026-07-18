import cobra
from cobra.flux_analysis import flux_variability_analysis
import pandas as pd

def assess_unbounded_fraction(model, threshold=999):
    """Runs FVA and returns the fraction of unbounded non-blocked reactions."""
    # Get active (non-blocked) reactions
    fva_res = flux_variability_analysis(model, fraction_of_optimum=0)

    # Identify non-blocked reactions (flux > 1e-6 in either direction)
    active_mask = (fva_res['maximum'].abs() > 1e-6) | (fva_res['minimum'].abs() > 1e-6)
    active_reactions = fva_res[active_mask]

    if len(active_reactions) == 0:
        return 0.0, []

    # Find unbounded reactions (hitting the +/- 1000 bounds)
    unbounded_mask = (active_reactions['maximum'] >= threshold) | (active_reactions['minimum'] <= -threshold)
    unbounded_rxns = active_reactions[unbounded_mask].index.tolist()

    fraction = (len(unbounded_rxns) / len(active_reactions)) * 100
    return fraction, unbounded_rxns

# ==========================================
# 1. Load Your Model
# ==========================================
# Replace 'your_draft_model.xml' with your actual file path
model = cobra.io.read_sbml_model("../models/model_3.xml")

print(f"Original Model: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites.")

# ==========================================
# 2. Assess Initial Unbounded Flux (Baseline)
# ==========================================
print("\nRunning initial FVA (this may take a moment)...")
initial_fraction, initial_unbounded = assess_unbounded_fraction(model)
print(f"Initial unbounded fraction: {initial_fraction:.2f}% ({len(initial_unbounded)} reactions)")

# ==========================================
# 3. Apply Targeted Fixes
# ==========================================


# Fix A: Block Lumped / Duplicate Reactions
# Based on your Memote output, these are likely redundant duplicates causing loops.
duplicates_to_block = [
    "ACONT",     # Lumped aconitase (keep ACONTa and ACONTb)
    "ACTD_1",    # Duplicate
    "ECOAH1_1",  # Duplicate
    "GLUDy_1",   # Duplicate
    "HACD1_1",   # Duplicate
    "HACD1i",    # Duplicate
    "HACD8i",    # Duplicate
    "NADDPp_1",  # Duplicate
    "PYNP1_1",   # Duplicate
    "DURIPP_1",  # Duplicate
    "PGMT_2",    # Duplicate
    "PGMT_B",    # Duplicate
    "MTAM_1",    # Duplicate
    "OCBT_1",    # Duplicate
    "PTHPS_1",   # Duplicate
    "UAGPT2_1",  # Duplicate
    "r2465_1"    # Typical tool-generated artifact
]

blocked_count = 0
for rxn_id in duplicates_to_block:
    if rxn_id in model.reactions:
        model.reactions.get_by_id(rxn_id).bounds = (0, 0)
        blocked_count += 1
print(f"\nBlocked {blocked_count} suspected duplicate/lumped reactions.")

# Fix B: Constrain Free-Energy Transporters
# Protons and water shuttling freely across compartments cause infinite ATP/NADH loops.
# We make these strictly irreversible (outward or inward only) rather than freely reversible.
directional_transporters = [
    "Htex", "H2Otpp", "CO2t", "CO2tex", "NH4t", "NH4tex"
]

fixed_transports = 0
for rxn_id in directional_transporters:
    if rxn_id in model.reactions:
        rxn = model.reactions.get_by_id(rxn_id)
        # If it's reversible (-1000, 1000), force it to be forward only (0, 1000)
        # You may need to tweak this to (-1000, 0) depending on the actual physiological direction.
        if rxn.lower_bound < 0 and rxn.upper_bound > 0:
            rxn.lower_bound = 0
            fixed_transports += 1
print(f"Fixed directionality for {fixed_transports} transport reactions.")

# ==========================================
# 4. Assess Improvement (Post-Fix)
# ==========================================
print("\nRunning post-fix FVA (this may take a moment)...")
new_fraction, new_unbounded = assess_unbounded_fraction(model)
print(f"New unbounded fraction: {new_fraction:.2f}% ({len(new_unbounded)} reactions)")

improvement = initial_fraction - new_fraction
print(f"Total improvement: {improvement:.2f}% drop in unbounded reactions.")

# Optional: Save the fixed model
cobra.io.write_sbml_model(model, "../models/model_4.xml")
# print("Saved fixed model to 'fixed_draft_model.xml'")
