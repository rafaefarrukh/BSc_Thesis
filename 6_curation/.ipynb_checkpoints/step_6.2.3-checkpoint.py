import pandas as pd
import cobra

def update_gem_annotations(model_path, csv_path, output_path):
    # Load the metabolic model (assumes SBML format)
    model = cobra.io.read_sbml_model(model_path)

    # Load the annotation data
    df = pd.read_csv(csv_path)
    df = df.fillna('')

    # Map CSV column names to memote-compatible MIRIAM prefixes
    # Modify these keys to match the exact column names generated in gene_mapping.csv
    prefix_mapping = {
        'RefSeq': 'refseq',
        'UniProt': 'uniprot',
        'KEGG': 'kegg.genes',
        'GO': 'go',
        'COG': 'cog',
        'UniRef50': 'uniref',
        'UniRef90': 'uniref'
    }

    # Iterate through each gene in the model
    for gene in model.genes:
        # Match the model gene ID against Gene_ID or Locus_Tag in the CSV
        match = df[(df['Gene_ID'] == gene.id) | (df['Locus_Tag'] == gene.id)]

        if not match.empty:
            row = match.iloc[0]

            for csv_col, miriam_prefix in prefix_mapping.items():
                if csv_col in df.columns and row[csv_col]:
                    # Split concatenated entries and clean whitespace
                    values = [v.strip() for v in str(row[csv_col]).split('|')]

                    # Cobrapy/Memote expect strings for single entries, lists for multiple
                    if len(values) == 1:
                        gene.annotation[miriam_prefix] = values[0]
                    else:
                        gene.annotation[miriam_prefix] = values

    # Save the updated model
    cobra.io.write_sbml_model(model, output_path)

if __name__ == "__main__":
    update_gem_annotations('../models/model_3.xml', 'gene_mapping.csv', '../models/model_3.xml')
