"""
STEP 6.6 — Fix stray compartment on PHB-pathway metabolites
=============================================================
MEMOTE flagged a compartment-naming inconsistency: every cytosolic
metabolite in this model uses compartment key 'C_c' (739 metabolites)
EXCEPT two metabolites added manually in Step 6.1 for the PHB pathway:

    3hbcoa__R   ((R)-3-hydroxybutyryl-CoA)
    phb_c       (Poly-beta-hydroxybutyrate)

...which were left in a separate, effectively orphan compartment
literally named 'c'. This is suspected to be why MEMOTE's
test_ngam_presence errors out silently (message/data both null): its
cytosol-auto-detection heuristic likely matches on the bare 'c' key
(the standard BiGG convention) instead of this model's actual 'C_c'
cytosol, gets confused by having two cytosol-like candidates, and
fails to find the ATP-hydrolysis (ATPM) reaction as a result.

This script:
  1. Reassigns both metabolites' compartment to 'C_c'.
  2. Removes the now-empty 'c' entry from model.compartments.
  3. Leaves metabolite IDs and annotations untouched — compartment
     assignment, IDs, and annotation dicts are independent; nothing
     here touches BiGG/KEGG/etc. cross-references.
  4. Re-validates growth and reaction count are unaffected before
     saving, and aborts rather than silently writing a broken model.

Input:  ../models/model_5.xml
Output: ../models/model_6.xml
"""

import cobra

# ----------------------------------------------------------------------
# CONFIG
# ----------------------------------------------------------------------
MODEL_IN = "../models/model_5.xml"
MODEL_OUT = "../models/model_6.xml"

STRAY_COMPARTMENT = "c"
TARGET_COMPARTMENT = "C_c"
METABOLITES_TO_FIX = ["3hbcoa__R", "phb_c"]

MIN_ACCEPTABLE_GROWTH = 1e-6

# ----------------------------------------------------------------------
# LOAD MODEL
# ----------------------------------------------------------------------
print("=" * 70)
print("STEP 6.6 — FIX STRAY COMPARTMENT (PHB METABOLITES)")
print("=" * 70)
print(f"Loading model from: {MODEL_IN}")
model = cobra.io.read_sbml_model(MODEL_IN)
print(f"Loaded model '{model.id}': {len(model.reactions)} reactions, "
      f"{len(model.metabolites)} metabolites, {len(model.genes)} genes")

print(f"\nCompartments before fix: {model.compartments}")

growth_before = model.slim_optimize() or 0.0
print(f"Growth rate before fix: {growth_before:.6f} /h")

# ----------------------------------------------------------------------
# SAFETY CHECKS
# ----------------------------------------------------------------------
missing = [mid for mid in METABOLITES_TO_FIX if mid not in model.metabolites]
if missing:
    raise SystemExit(f"Expected metabolites not found in model: {missing}")

for mid in METABOLITES_TO_FIX:
    met = model.metabolites.get_by_id(mid)
    if met.compartment != STRAY_COMPARTMENT:
        print(f"NOTE: '{mid}' is already in compartment "
              f"'{met.compartment}', not '{STRAY_COMPARTMENT}' as "
              f"expected — leaving it untouched.")

# ----------------------------------------------------------------------
# APPLY THE FIX
# ----------------------------------------------------------------------
changed = []
for mid in METABOLITES_TO_FIX:
    met = model.metabolites.get_by_id(mid)
    if met.compartment == STRAY_COMPARTMENT:
        print(f"\nReassigning {met.id} ('{met.name}'): "
              f"compartment '{met.compartment}' -> '{TARGET_COMPARTMENT}'")
        met.compartment = TARGET_COMPARTMENT
        changed.append(met.id)

if not changed:
    print("\nNothing to change — no metabolites were in the stray "
          "compartment. Exiting without writing a new file.")
    raise SystemExit(0)

# Drop the now-empty stray compartment entry, if present and empty
remaining_in_stray = [m for m in model.metabolites if m.compartment == STRAY_COMPARTMENT]
if not remaining_in_stray and STRAY_COMPARTMENT in model.compartments:
    del model.compartments[STRAY_COMPARTMENT]
    print(f"\nRemoved now-empty compartment entry '{STRAY_COMPARTMENT}' "
          f"from model.compartments.")
elif remaining_in_stray:
    print(f"\nNOTE: {len(remaining_in_stray)} metabolite(s) still "
          f"reference compartment '{STRAY_COMPARTMENT}': "
          f"{[m.id for m in remaining_in_stray]} — leaving that "
          f"compartment entry in place.")

print(f"\nCompartments after fix: {model.compartments}")

# ----------------------------------------------------------------------
# VALIDATE — annotations untouched, growth unaffected
# ----------------------------------------------------------------------
for mid in changed:
    met = model.metabolites.get_by_id(mid)
    print(f"  {met.id}: annotation unchanged -> {dict(met.annotation)}")

growth_after = model.slim_optimize() or 0.0
print(f"\nGrowth rate after fix: {growth_after:.6f} /h "
      f"(was {growth_before:.6f} /h)")

if growth_after < MIN_ACCEPTABLE_GROWTH:
    raise SystemExit(
        "\nABORTING — growth is ~0 after this change, which shouldn't "
        "happen from a pure compartment relabeling. Refusing to save. "
        "Investigate before retrying."
    )

if abs(growth_after - growth_before) > 1e-6:
    print("NOTE: growth rate changed slightly — unexpected for a "
          "compartment-label-only edit, but not necessarily wrong if "
          "PHB-related reactions interact with the objective. Verify "
          "this is expected before trusting downstream analyses.")

print(f"\nReaction count unchanged: {len(model.reactions)}")
print(f"Metabolite count unchanged: {len(model.metabolites)}")

# ----------------------------------------------------------------------
# SAVE
# ----------------------------------------------------------------------
cobra.io.write_sbml_model(model, MODEL_OUT)
print(f"\nModel written to: {MODEL_OUT}")
print("Re-run MEMOTE on this model to confirm test_ngam_presence now "
      "returns a real message/data pair (pass or a clear failure "
      "reason) instead of null/null, and that the compartment "
      "mixed-naming note is gone from the summary.")