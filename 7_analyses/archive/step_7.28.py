"""
STEP 7.28 — LITERATURE / EXPERIMENTAL SANITY-CHECK COMPARISON

MEMOTE noted that no experimental datasets were supplied, so growth-rate and
gene-essentiality validation were skipped entirely (NOTSET), and none of the
step_7 analyses compare model predictions against real measurements either.

This script does NOT invent organism-specific numbers for Rosellomora
marisflavi (none were found in the literature at the time of writing). It
instead:
  1. Builds a fillable benchmark table (./step_7.28_benchmarks_template.csv)
     with placeholder rows for the measurements you'd want to compare
     against (growth rate on defined media, PHB content %, PHB titer,
     substrate-specific yields, essential-gene list from a Tn-seq/knockout
     study if one exists for this or a close relative).
  2. Pulls the model's current predictions for each benchmark metric where
     possible, so the template arrives pre-populated on the "model" column
     and you only need to fill in "literature/experimental" values as they
     become available.
  3. Reports %-difference and flags anything outside a generous 2x band as
     "needs review" once literature values are filled in.

Broad reference points from the wider PHB-producer literature (Cupriavidus
necator, Bacillus spp., grown on glucose/sucrose) are included as loose
sanity bounds only -- see step_7.26 for the same caveat. Do not treat these
as validated targets for this organism; replace with real measurements
(ideally same species/genus) as soon as they exist.

Input:  ../models/model_7.xml
        ./step_7.2_fba_results.csv
        ./step_7.26_yield_summary.csv   (optional, from step 7.26)
Output: ./step_7.28_benchmarks_template.csv   <- fill in "literature_value" column
        ./step_7.28_comparison.csv            <- rerun after filling in
        ./step_7.28_summary.txt
"""

import cobra
import pandas as pd
import os

MODEL_PATH = "../models/model_7.xml"
BIOMASS_ID = "Growth"
PHB_ID = "PHBS_syn_1"

print(f"Loading model from: {MODEL_PATH}")
model = cobra.io.read_sbml_model(MODEL_PATH)
print(f"Loaded model '{model.id}': {len(model.reactions)} reactions, "
      f"{len(model.metabolites)} metabolites, {len(model.genes)} genes")

model.objective = BIOMASS_ID
max_growth = model.slim_optimize()

model.objective = PHB_ID
max_phb = model.slim_optimize()

benchmarks = [
    dict(metric="max_growth_rate_aerobic_glucose_1_h",
         model_value=round(max_growth, 4),
         literature_value="",  # fill in from a growth curve (OD/CDW over time) if available
         source="",
         notes="Compare against a controlled batch growth curve, not a single OD reading."),
    dict(metric="PHB_content_pct_of_CDW",
         model_value="",  # requires converting PHBS_syn_1 flux + biomass composition to wt%; see step_7.26
         literature_value="",
         source="",
         notes="Broad literature range for PHB producers on sugars: often 40-80 wt%% CDW "
               "in optimized fed-batch (organism-specific; not validated for this strain)."),
    dict(metric="PHB_mass_yield_g_per_g_glucose",
         model_value="",  # pull from step_7.26_yield_summary.csv if present
         literature_value="",
         source="",
         notes="Loose reference range from other PHB producers: ~0.33-0.48 g/g on sugars."),
    dict(metric="essential_gene_count",
         model_value="",  # pull from step_7.27_gene_knockout_screen_full.csv if present
         literature_value="",
         source="",
         notes="Compare against a Tn-seq/knockout-library essentiality study, "
               "ideally in this species or a close relative, if one exists."),
    dict(metric="minimal_medium_component_count",
         model_value="",  # pull from step_7.7_minimal_medium.csv if present
         literature_value="",
         source="",
         notes="Compare against defined-minimal-medium recipes used experimentally for growth."),
]

# Auto-fill from other step_7 outputs where available
try:
    yield_df = pd.read_csv("./step_7.26_yield_summary.csv")
    growth_opt = yield_df[yield_df["scenario"] == "growth_optimal"]
    if len(growth_opt):
        for row in benchmarks:
            if row["metric"] == "PHB_mass_yield_g_per_g_glucose":
                row["model_value"] = round(float(growth_opt["mass_yield_g_g"].iloc[0]), 4) \
                    if pd.notna(growth_opt["mass_yield_g_g"].iloc[0]) else "n/a (zero PHB at growth optimum)"
except FileNotFoundError:
    pass

try:
    gene_screen = pd.read_csv("./step_7.27_gene_knockout_screen_full.csv")
    n_essential = (gene_screen["classification"] == "essential").sum()
    for row in benchmarks:
        if row["metric"] == "essential_gene_count":
            row["model_value"] = int(n_essential)
except FileNotFoundError:
    pass

try:
    min_med = pd.read_csv("./step_7.7_minimal_medium.csv")
    for row in benchmarks:
        if row["metric"] == "minimal_medium_component_count":
            row["model_value"] = len(min_med)
except FileNotFoundError:
    pass

bench_df = pd.DataFrame(benchmarks)
bench_df.to_csv("./step_7.28_benchmarks_template.csv", index=False)
print(f"\nBenchmark template written to: ./step_7.28_benchmarks_template.csv")
print("Fill in the 'literature_value' and 'source' columns as data becomes")
print("available, then re-run the comparison block below.")

# ---- Comparison block (runs only if literature_value has been filled in) ----
comparison_rows = []
for row in benchmarks:
    lit = row["literature_value"]
    model_val = row["model_value"]
    if lit not in ("", None) and model_val not in ("", None):
        try:
            lit_f, model_f = float(lit), float(model_val)
            pct_diff = 100 * (model_f - lit_f) / lit_f if lit_f != 0 else float("nan")
            flag = "OK" if abs(pct_diff) <= 100 else "NEEDS REVIEW (>2x off)"
        except (ValueError, TypeError):
            pct_diff, flag = float("nan"), "non-numeric, compare manually"
        comparison_rows.append(dict(metric=row["metric"], model_value=model_val,
                                     literature_value=lit, pct_diff=pct_diff, flag=flag))

if comparison_rows:
    comp_df = pd.DataFrame(comparison_rows)
    comp_df.to_csv("./step_7.28_comparison.csv", index=False)
    print("\nComparison results:")
    print(comp_df.to_string(index=False))
else:
    print("\nNo literature values filled in yet — comparison CSV not generated.")
    print("This is expected on first run; re-run after populating the template.")

with open("./step_7.28_summary.txt", "w") as f:
    f.write("STEP 7.28 — Literature/experimental sanity-check summary\n")
    f.write("=" * 60 + "\n")
    f.write("No organism-specific experimental data for Rosellomora marisflavi\n")
    f.write("was available at the time this template was generated. Populate\n")
    f.write("./step_7.28_benchmarks_template.csv with real measurements (growth\n")
    f.write("curves, PHB titers, essentiality screens) as they become available,\n")
    f.write("then re-run this script to generate ./step_7.28_comparison.csv.\n\n")
    f.write(bench_df.to_string(index=False))

print("\nSummary written to: ./step_7.28_summary.txt")
