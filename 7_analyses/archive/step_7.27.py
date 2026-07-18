"""
STEP 7.27 — GENE-LEVEL ESSENTIALITY SCREEN

step_7.20 screened reactions, but 208/1963 reactions (10.6%) have no GPR
rule (per step_7.1), so a reaction-level knockout for those isn't
interpretable as a genetic perturbation -- there's no gene to actually
delete. This script runs single_gene_deletion directly on the model's 1159
genes, which correctly:
  - only perturbs reactions that are actually controlled by that gene's GPR
  - accounts for isozymes/complexes (a gene knockout only kills a reaction
    if it's on every enabling path of that reaction's GPR Boolean rule)
  - reports which genes are essential/critical/non-essential in a way that
    maps directly onto real knockout-strain feasibility

Input:  ../models/model_7.xml
        ./step_7.20_knockout_screen_full.csv   (reaction-level screen, for comparison)
Output: ./step_7.27_gene_knockout_screen_full.csv
        ./step_7.27_essential_genes.csv
        ./step_7.27_reaction_vs_gene_essentiality_comparison.csv
        ./step_7.27_summary.txt
"""

import cobra
import pandas as pd
from cobra.flux_analysis import single_gene_deletion

MODEL_PATH = "../models/model_7.xml"
BIOMASS_ID = "Growth"
GROWTH_CUTOFF_ESSENTIAL = 1e-6      # fraction of WT below which a gene is "essential"
GROWTH_CUTOFF_CRITICAL = 0.10       # fraction of WT below which (but above essential) is "critical"

print(f"Loading model from: {MODEL_PATH}")
model = cobra.io.read_sbml_model(MODEL_PATH)
print(f"Loaded model '{model.id}': {len(model.reactions)} reactions, "
      f"{len(model.metabolites)} metabolites, {len(model.genes)} genes")

model.objective = BIOMASS_ID
wt_growth = model.slim_optimize()
print(f"Wild-type growth rate: {wt_growth:.6f} /h")

genes_with_gpr = {g.id for g in model.genes if len(g.reactions) > 0}
print(f"Genes associated with at least one reaction: {len(genes_with_gpr)} / {len(model.genes)}")

print("\nRunning single-gene knockout screen "
      f"(this can take a while on a ~{len(model.genes)}-gene model)...")
gene_result = single_gene_deletion(model, gene_list=model.genes, processes=1)

# cobra returns a DataFrame with 'ids' (frozenset) and 'growth', 'status'
gene_result = gene_result.reset_index(drop=True)
gene_result["gene_id"] = gene_result["ids"].apply(lambda s: next(iter(s)))
gene_result["growth_fraction_of_wt"] = gene_result["growth"] / wt_growth


def classify(frac):
    if pd.isna(frac) or frac < GROWTH_CUTOFF_ESSENTIAL:
        return "essential"
    elif frac < GROWTH_CUTOFF_CRITICAL:
        return "critical"
    else:
        return "non-essential"


gene_result["classification"] = gene_result["growth_fraction_of_wt"].apply(classify)
gene_result = gene_result[["gene_id", "growth", "growth_fraction_of_wt", "classification", "status"]]
gene_result = gene_result.sort_values("growth_fraction_of_wt")
gene_result.to_csv("./step_7.27_gene_knockout_screen_full.csv", index=False)

counts = gene_result["classification"].value_counts()
print("\nGene classification counts:")
print(counts.to_string())

essential_genes = gene_result[gene_result["classification"] == "essential"]
essential_genes.to_csv("./step_7.27_essential_genes.csv", index=False)
print(f"\n{len(essential_genes)} essential genes written to: ./step_7.27_essential_genes.csv")

# ---- Compare against reaction-level screen from step_7.20 ----
try:
    rxn_screen = pd.read_csv("./step_7.20_knockout_screen_full.csv")
    n_no_gpr = sum(len(model.reactions.get_by_id(rid).genes) == 0
                   for rid in rxn_screen["reaction_id"] if rid in model.reactions)
    essential_rxns = set(rxn_screen.loc[rxn_screen["classification"] == "essential", "reaction_id"])

    # Map essential genes -> the reactions they control, and flag whether
    # that reaction was also flagged essential at the reaction level.
    comp_rows = []
    for _, row in essential_genes.iterrows():
        gene = model.genes.get_by_id(row["gene_id"])
        for rxn in gene.reactions:
            comp_rows.append(dict(
                gene_id=row["gene_id"],
                reaction_id=rxn.id,
                reaction_essential_in_step_7_20=rxn.id in essential_rxns,
                has_gpr=len(rxn.genes) > 0,
            ))
    comp_df = pd.DataFrame(comp_rows)
    comp_df.to_csv("./step_7.27_reaction_vs_gene_essentiality_comparison.csv", index=False)

    print(f"\nReactions in step_7.20 screen with no GPR (not genetically interpretable): {n_no_gpr}")
    if len(comp_df):
        mismatch = comp_df[~comp_df["reaction_essential_in_step_7_20"]]
        print(f"Gene-essential reactions NOT flagged essential at reaction level "
              f"(likely isozyme/complex cases): {len(mismatch)}")
except FileNotFoundError:
    print("\nNOTE: ./step_7.20_knockout_screen_full.csv not found — skipping comparison.")
    comp_df = None

with open("./step_7.27_summary.txt", "w") as f:
    f.write("STEP 7.27 — Gene-level essentiality summary\n")
    f.write("=" * 60 + "\n")
    f.write(f"Wild-type growth: {wt_growth:.6f} /h\n")
    f.write(f"Total genes: {len(model.genes)} ({len(genes_with_gpr)} with >=1 reaction)\n\n")
    f.write(counts.to_string())
    f.write("\n\nThresholds: essential < {:.0e} of WT; critical < {:.2f} of WT\n"
            .format(GROWTH_CUTOFF_ESSENTIAL, GROWTH_CUTOFF_CRITICAL))
    if comp_df is not None:
        f.write(f"\nGene-essential-but-reaction-not-flagged-essential cases: "
                f"{len(comp_df[~comp_df['reaction_essential_in_step_7_20']]) if len(comp_df) else 'n/a'}\n"
                "(These typically indicate isozyme/complex GPR logic where the reaction\n"
                "survives most single-gene knockouts but not all of them.)\n")

print("\nSummary written to: ./step_7.27_summary.txt")
