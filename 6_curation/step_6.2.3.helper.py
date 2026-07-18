import csv

def convert_gff3_to_csv(input_filepath, output_filepath):
    parsed_records = []
    database_columns = set()

    with open(input_filepath, 'r') as file:
        for line in file:
            # Skip comment lines and empty lines
            if line.startswith('#') or not line.strip():
                continue

            columns = line.strip().split('\t')
            if len(columns) != 9:
                continue

            feature_type = columns[2]
            attributes_raw = columns[8]

            # Parse the attributes column (column 9) into a dictionary
            attributes = {}
            for attr in attributes_raw.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attributes[key] = value

            # Initialize row with core identifiers
            row = {
                'Gene_ID': attributes.get('ID', ''),
                'Locus_Tag': attributes.get('locus_tag', ''),
                'Feature_Type': feature_type,
                'Name': attributes.get('Name', '')
            }

            # Parse the Dbxref string for all database IDs
            if 'Dbxref' in attributes:
                db_references = attributes['Dbxref'].split(',')
                for db_ref in db_references:
                    if ':' in db_ref:
                        db_name, db_id = db_ref.split(':', 1)
                        database_columns.add(db_name)

                        # Concatenate multiple IDs from the same database using a pipe delimiter
                        if db_name in row:
                            row[db_name] += f"|{db_id}"
                        else:
                            row[db_name] = db_id

            parsed_records.append(row)

    # Construct headers ensuring core identifiers are first, followed by alphabetized database names
    headers = ['Gene_ID', 'Locus_Tag', 'Feature_Type', 'Name'] + sorted(list(database_columns))

    # Write the parsed data to the output CSV
    with open(output_filepath, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()
        writer.writerows(parsed_records)

if __name__ == "__main__":
    # Specify the input and output file names
    convert_gff3_to_csv('../2_bakta/annotated_genome.gff3', 'gene_mapping.csv')
