"""
STEP 7.31 — MASS/CHARGE BALANCE CHECK, PHB-PATHWAY NEIGHBORHOOD ONLY

step_7.1 found 318/1711 non-boundary reactions (18.6%) are mass/charge
imbalanced model-wide. That number alone doesn't tell you whether the flux
values you actually care about (PHB synthesis and its immediate neighbors,
as used throughout steps 7.3, 7.12, 7.13, 7.16, 7.23) are trustworthy. This
script cross-references the global imbalance list against the PHB-distance
table from step_7.12 and reports imbalances specifically within 0-3 hops of
PHB synthesis, since that's the neighborhood most of the mechanistic
(non-global) analyses in this pipeline actually rely on.

Input:  ../models/model_7.xml
        ./step_7.1_mass_charge_imbalances.csv
        ./step_7.12_phb_pathway_distances.csv
Output: ./step_7.31_phb_neighborhood_imbalances.csv
        ./step_7.31_summary.txt
"""

import cobra
import pandas as pd

MODEL_PATH = "../models/model_7.xml"
MAX_DISTANCE = 3  # matches the "core neighborhood" used by steps 7.13-7.15/7.18

print(f"Loading model from: {MODEL_PATH}")
model = cobra.io.read_sbml_model(MODEL_PATH)
print(f"Loaded model '{model.id}': {len(model.reactions)} reactions, "
      f"{len(model.metabolites)} metabolites, {len(model.genes)} genes")

print("=" * 70)
print("STEP 7.31 — PHB-NEIGHBORHOOD MASS/CHARGE BALANCE CHECK")
print("=" * 70)

try:
    imbalances = pd.read_csv("./step_7.1_mass_charge_imbalances.csv")
except FileNotFoundError:
    print("./step_7.1_mass_charge_imbalances.csv not found — recomputing "
          "mass/charge balance directly from the model instead.")
    rows = []
    for rxn in model.reactions:
        if rxn.boundary:
            continue
        try:
            balance = rxn.check_mass_balance()
        except Exception:
            balance = {"note": "could not evaluate (missing formula/charge)"}
        if balance:
            rows.append(dict(reaction_id=rxn.id, reaction_name=rxn.name, imbalance=str(balance)))
    imbalances = pd.DataFrame(rows)

try:
    distances = pd.read_csv("./step_7.12_phb_pathway_distances.csv")
except FileNotFoundError:
    raise SystemExit(
        "ERROR: ./step_7.12_phb_pathway_distances.csv is required for this script "
        "(run step_7.12 first)."
    )

merged = imbalances.merge(distances[["reaction_id", "distance_from_PHB"]],
                           on="reaction_id", how="inner")
neighborhood = merged[merged["distance_from_PHB"] <= MAX_DISTANCE].sort_values("distance_from_PHB")
neighborhood.to_csv("./step_7.31_phb_neighborhood_imbalances.csv", index=False)

n_total_imbalanced = len(imbalances)
n_in_neighborhood = len(distances[distances["distance_from_PHB"] <= MAX_DISTANCE])
n_imbalanced_in_neighborhood = len(neighborhood)

print(f"\nTotal mass/charge-imbalanced reactions (model-wide): {n_total_imbalanced}")
print(f"Reactions within {MAX_DISTANCE} hops of PHB synthesis: {n_in_neighborhood}")
print(f"Of those, imbalanced: {n_imbalanced_in_neighborhood} "
      f"({100 * n_imbalanced_in_neighborhood / max(n_in_neighborhood, 1):.1f}%)")

if n_imbalanced_in_neighborhood > 0:
    print("\nImbalanced reactions closest to PHB synthesis:")
    print(neighborhood.head(20).to_string(index=False))
else:
    print("\nNo imbalances found within the PHB pathway neighborhood — "
          "flux values from steps 7.3/7.13/7.16/7.23 are not undermined by "
          "mass/charge balance issues at their core.")

# Explicitly check the PHB reaction and its two closest neighbors
core_ids = ["PHBS_syn_1", "AACOAR_syn"]
for rid in core_ids:
    if rid in model.reactions:
        rxn = model.reactions.get_by_id(rid)
        try:
            bal = rxn.check_mass_balance()
        except Exception as e:
            bal = f"could not evaluate: {e}"
        status = "BALANCED" if not bal else f"IMBALANCED: {bal}"
        print(f"\nCore reaction check — {rid} ({rxn.name}): {status}")

with open("./step_7.31_summary.txt", "w") as f:
    f.write("STEP 7.31 — PHB-neighborhood mass/charge balance summary\n")
    f.write("=" * 60 + "\n")
    f.write(f"Total mass/charge-imbalanced reactions (model-wide): {n_total_imbalanced}\n")
    f.write(f"Reactions within {MAX_DISTANCE} hops of PHB synthesis: {n_in_neighborhood}\n")
    f.write(f"Of those, imbalanced: {n_imbalanced_in_neighborhood} "
            f"({100 * n_imbalanced_in_neighborhood / max(n_in_neighborhood, 1):.1f}%)\n\n")
    if n_imbalanced_in_neighborhood > 0:
        f.write("Interpretation: flux values for the imbalanced reactions listed in\n"
                "step_7.31_phb_neighborhood_imbalances.csv (and anything directly\n"
                "downstream of them) should be treated as less reliable until their\n"
                "metabolite formulas/charges are corrected in the model.\n")
    else:
        f.write("Interpretation: the PHB pathway core is mass/charge-balanced, so\n"
                "flux magnitudes reported near PHB synthesis in steps 7.3, 7.13,\n"
                "7.16, and 7.23 are not undermined by this particular model issue\n"
                "(though the loopless-cycle caveat from step_7.25 still applies).\n")

print("\nSummary written to: ./step_7.31_summary.txt")
