###################################################################################################
# STEP 6.2.1: Annotate model with BiGG IDs
###################################################################################################

# libraries
import re
import cobra
import pandas as pd

# paths
MET_TABLE = 'BIGG/bigg_models_metabolites.txt'
REAC_TABLE = 'BIGG/bigg_models_reactions.txt'
MODEL_IN = '../models/model_2.xml'
MODEL_OUT = '../models/model_3.xml'

# build lookup tables for BiGG IDs from the BiGG tables
met_df = pd.read_csv(MET_TABLE, sep='\t')
met_lookup = {}
for _, row in met_df.dropna(subset=['old_bigg_ids']).iterrows():
    for oid in str(row['old_bigg_ids']).split(';'):
        oid = oid.strip()
        if oid:
            met_lookup[oid] = row['bigg_id']

reac_df = pd.read_csv(REAC_TABLE, sep='\t')
reac_lookup = {}
for _, row in reac_df.dropna(subset=['old_bigg_ids']).iterrows():
    for oid in str(row['old_bigg_ids']).split(';'):
        oid = oid.strip()
        if oid:
            reac_lookup[oid] = row['bigg_id']

print(f"Loaded {len(met_lookup)} metabolite ID mappings, {len(reac_lookup)} reaction ID mappings")

# import model
model = cobra.io.read_sbml_model(MODEL_IN)

# annaotate metabolites (compartment letter stripped)
met_hits = 0
for met in model.metabolites:
    clean_id = re.sub(r'_[a-z]$', '', met.id)
    if clean_id in met_lookup:
        met.annotation['bigg.metabolite'] = [met_lookup[clean_id]]
        met_hits += 1

# annotate reactions (compartment letter stripped)
reac_hits = 0
for rxn in model.reactions:
    clean_id = re.sub(r'_[a-z]$', '', rxn.id)
    if clean_id in reac_lookup:
        rxn.annotation['bigg.reaction'] = [reac_lookup[clean_id]]
        reac_hits += 1

print(f"Annotated {met_hits}/{len(model.metabolites)} metabolites, "
      f"{reac_hits}/{len(model.reactions)} reactions")

# export model
cobra.io.write_sbml_model(model, MODEL_OUT)
print(f"Annotated model saved to {MODEL_OUT}")

###################################################################################################
