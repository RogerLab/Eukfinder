"""
Translate Plast results to taxonomy using acc2tax.
Version 1.1, Dec 15, 2024
Author: D.Z.

================================================
python3 Step2.2_Parsing_Plast_results.py [-i <input_file> -d <database_path> -o <output_file>]
================================================

This script performs the following steps:
1. Reads Plast results and extracts accession numbers from the input file.
2. Uses the acc2tax tool to translate accession numbers to taxonomy (domain, phylum, genus, species).
3. Combines the original Plast results with taxonomy data.
4. Outputs a combined file containing Plast results with taxonomy information.
5. Reports eukaryotic species with more than 10 contigs.

Report Example:
   species                           Plast_count
   Blastocystis sp. subtype 4          3300
   Cyclospora cayetanensis              15
"""

import argparse
import subprocess
import pandas as pd
import os
from collections import OrderedDict
import re


def parse_acc_numbers(input_file, temp_file):
    """
    Parse accession numbers from the input file and write them to a temporary file.

    Args:
        input_file (str): Path to the Plast results file.
        temp_file (str): Path to the temporary file for storing accession numbers.

    Returns:
        df (pd.DataFrame): Parsed DataFrame with key Plast columns.
        acc_dict (dict): Dictionary mapping query IDs to subject accession IDs.
    """
    # Read Plast results into a DataFrame
    df = pd.read_csv(input_file, sep='\t', header=None)

    # Determine column structure based on the number of columns in the file
    if len(df.columns) == 12:
        df.columns = ['query ID', 'subject ID', 'ident', 'align_length', 'nb. misses',
                      'nb. gaps', 'query begin', 'query end', 'subject begin', 'subject end',
                      'e-value', 'bit score']
    elif len(df.columns) == 23:
        df.columns = ['query ID', 'subject ID', 'ident', 'align_length', 'nb. misses', 'nb. gaps',
                      'query begin', 'query end', 'subject begin', 'subject end', 'e-value', 'bit score',
                      'exact_align', 'query length', 'query frame', 'query translated', 'query coverage',
                      'nb. gaps in query', 'subject length', 'subject frame', 'subject translated',
                      'subject coverage', 'nb. gaps in subject']

    # Filter relevant columns
    df = df[['query ID', 'subject ID', 'ident', 'align_length', 'e-value', 'bit score']]

    # Retain the first occurrence of each query ID
    df = df.groupby('query ID', as_index=False).first()

    # Create a dictionary mapping query IDs to subject accession numbers
    acc_dict = {row['query ID']: row['subject ID'].split('.')[0] for _, row in df.iterrows()}

    # Write accession numbers to a temporary file
    with open(temp_file, 'w') as outfile:
        outfile.writelines(f"{acc_id}\n" for acc_id in acc_dict.values())

    return df, acc_dict


def run_acc2tax(database, temp_file, raw_output):
    """
    Execute acc2tax to retrieve taxonomy.

    Args:
        database (str): Path to the acc2tax database.
        temp_file (str): Temporary file with accession numbers.
        raw_output (str): Path to the raw taxonomy output file.
    """
    print("Running acc2tax...")
    try:
        subprocess.run(
            ["acc2tax", "-a", "--database", database, "--input", temp_file, "--output", raw_output],
            check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error running acc2tax: {e}")
        exit(1)


def parse_taxonomy_output(raw_output, acc_dict):
    """
    Parse the raw acc2tax output to extract taxonomy details.

    Args:
        raw_output (str): Path to the raw taxonomy output file.
        acc_dict (dict): Dictionary mapping query IDs to subject accession numbers.

    Returns:
        taxonomy_df (pd.DataFrame): DataFrame with taxonomy details (domain, phylum, genus, species).
    """
    taxonomy_data = {}
    with open(raw_output) as infile:
        for line in infile:
            parts = re.split(r'[\t|,]', line.strip())
            acc_id = parts[0]
            lineage = parts[2:]

            # Ensure taxonomy data has at least four fields
            if len(lineage) >= 4:
                domain, phylum, genus, species = lineage[0], lineage[1], lineage[-2], lineage[-1]
            else:
                domain, phylum, genus, species = (lineage + ['-'] * (4 - len(lineage)))

            taxonomy_data[acc_id] = [domain, phylum, genus, species]

    # Map taxonomy data back to query IDs
    result = {key: taxonomy_data.get(acc_id, ['-'] * 4) for key, acc_id in acc_dict.items()}

    # Convert to DataFrame
    return pd.DataFrame.from_dict(result, orient='index', columns=['Plast domain', 'Plast phylum', 'Plast genus', 'Plast species'])


def main(input_file, database, output_file):
    # Temporary file paths
    temp_file, raw_output = "acc2taxIN.tmp", "acc2taxOUT.raw"

    # Parse accession numbers and run acc2tax
    df, acc_dict = parse_acc_numbers(input_file, temp_file)
    run_acc2tax(database, temp_file, raw_output)

    # Parse taxonomy output and merge with the original DataFrame
    taxonomy_df = parse_taxonomy_output(raw_output, acc_dict)
    final_df = pd.merge(df, taxonomy_df, left_on='query ID', right_index=True, how='left')

    # Write the combined DataFrame to the output file
    final_df.to_csv(output_file, sep='\t', index=False)
    print(f"Taxonomy output written to: {output_file}")

    # Cleanup temporary files
    for f in [temp_file, raw_output]:
        try:
            os.remove(f)
        except OSError as e:
            print(f"Error removing {f}: {e}")

    # Count and report eukaryotic species with more than 10 contigs
    eukaryotic_counts = final_df[final_df['Plast domain'] == 'Eukaryota']['Plast species'].value_counts()
    filtered_counts = eukaryotic_counts[eukaryotic_counts > 10].reset_index()
    filtered_counts.columns = ['species', 'Plast_count']

    print("\nEukaryotic species with more than 10 contigs detected by Plast:")
    print(filtered_counts.to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Translate Plast results to taxonomy using acc2tax.",
        epilog="Usage: python3 Step2.2_Parsing_Plast_results.py [-i <input_file> -d <database_path> -o <output_file>]"
    )
    parser.add_argument("-i", "--input", required=False, help="Input Plast results file", default="Eukfinder_long.PLAST_nt.tsv")
    parser.add_argument("-d", "--database", required=False, help="Path to acc2tax database", default="/db1/extra-data-sets/Acc2tax/Acc2Tax_04_01_2024")
    parser.add_argument("-o", "--output", required=False, help="Output file path", default="parsed_Plast_Acc2tax_results.txt")

    args = parser.parse_args()
    main(args.input, args.database, args.output)
