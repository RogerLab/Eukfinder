"""
To translate Centrifuge result TaxIDs to SuperKingdom using NCBI taxonomy file.
Version 1.2, Dec 14, 2024
Author: D.Z.

================================================
python3 Step1_Parsing_centrifuge_results.py
================================================
This script performs the following tasks:
1. Copies the EUnk.fasta file from the TempEukfinder folder to the current directory, renaming it to Eukfinder_long.fasta.
2. Parses Centrifuge results to output a file "parsed_centrifuge_results_eukLong.txt".
3. Translates TaxIDs from the Centrifuge results into taxonomy details using the NCBI taxonomy database.
4. Identifies eukaryotic species with more than 10 contigs detected by Centrifuge and reports the results.

Report:
Eukaryotic species with more than 10 contigs detected by Centrifuge:

                  species          count
Blastocystis sp. subtype 4          3300
Cyclospora cayetanensis             15
"""

import re
import os
import glob
import shutil
import argparse
import pandas as pd
from Bio import SeqIO
from ete3 import NCBITaxa
from collections import OrderedDict

def copy_and_extract_headers():
    """
    Copy EUnk.fasta from TempEukfinder to the current folder, rename it to Eukfinder_long.fasta,
    and extract sequence headers.
    Returns:
        list_headers: List of sequence headers from the FASTA file.
    """
    try:
        # Locate the FASTA file in TempEukfinder
        fasta_file = glob.glob("TempEukfinder/*.EUnk.fasta")[0]
        destination_file = "Eukfinder_long.fasta"

        # Copy and rename the file
        shutil.copy(fasta_file, destination_file)
        print(f"File copied and renamed to: {destination_file}")

        # Extract headers from the FASTA file
        list_headers = [seq.id for seq in SeqIO.parse(open(destination_file), 'fasta')]
        return list_headers
    except FileNotFoundError:
        raise FileNotFoundError(f"FASTA file not found in TempEukfinder.")
    except Exception as e:
        raise Exception(f"An error occurred: {e}")


def read_centrifuge_results(file_path, list_headers):
    """
    Read Centrifuge result file into a dataframe and remove duplicates by readID.
    """
    df = pd.read_csv(file_path, sep='\t')

    df = df.loc[:,['readID', 'taxID', 'hitLength', 'queryLength']]

    df = df.drop_duplicates(subset="readID")

    # Filter the DataFrame based on the readID column
    filtered_df = df[df['readID'].isin(list_headers)]

    acc_dict = {row['readID']: row['taxID'] for _, row in filtered_df.iterrows()}

    if filtered_df.index.name != 'readID':
        filtered_df.set_index("readID", inplace=True)

    return filtered_df, acc_dict


def parse_taxonomy_output(acc_dict):
    """
    Parse taxonomy using the ETE3 NCBITaxa database.
    Args:
        acc_dict: Dictionary mapping readID to taxID.
    Returns:
        df_result: DataFrame containing taxonomy details for each readID.
    """

    ncbi = NCBITaxa()  
    taxonomy_data = {}

    for taxID in acc_dict.values():
        if taxID == 0:
            taxonomy_data[taxID] = ['-'] * 4
        else:
            lineage = ncbi.get_lineage(int(taxID))
            names = ncbi.get_taxid_translator(lineage)
            Lineage = [names[t] for t in lineage if names[t] not in ["root", "cellular organisms"]]

            # Ensure taxonomy data has four fields
            if len(Lineage) >= 4:
                taxonomy_data[taxID] = [Lineage[0], Lineage[1], Lineage[-2], Lineage[-1]]
            else:
                taxonomy_data[taxID] = Lineage + ['-'] * (4 - len(Lineage))


    # Map parsed taxonomy back to the original accession keys
    result = OrderedDict()
    for key, acc_id in acc_dict.items():
        result[key] = taxonomy_data.get(acc_id, ['-'] * 4)
 
    df_result = pd.DataFrame.from_dict(result, orient='index', columns=['centrifuge domain', 'centrifuge phylum', 'centrifuge genus', 'centrifuge species'])
    df_result = df_result.rename_axis('readID').reset_index()
    df_result.set_index('readID', drop=False)

    if df_result.index.name != 'readID':
        df_result.set_index("readID", inplace=True)
  
    return df_result
   

def translate_taxid_superkingdom():
    """
    Main function to translate TaxID to SuperKingdom.
    """
    parser = argparse.ArgumentParser(prog='script.py', description="Translate TaxID to SuperKingdom using NCBI taxonomy file", 
                                     epilog='''Usage example:
                                           python3 script.py -i centrifuge_results.txt -o output.txt''')
    parser.add_argument('-c', '--centrifuge', required=False, help="Path to Centrifuge result file")
    parser.add_argument("-d", "--database", required=False, help="Path to acc2tax database", default = "/db1/extra-data-sets/Acc2tax/Acc2Tax_04_01_2024")
    parser.add_argument('-o', '--output', required=False, help="Path to output file", default="parsed_centrifuge_results_eukLong.txt")
    args = parser.parse_args()

    # Centrifuge result file
    if args.centrifuge:
        centrifuge_file = args.centrifuge
    else:
        centrifuge_file = glob.glob("tmps_*/*Eukfinder_long_centrifuge_UP")[0]

    headers = copy_and_extract_headers()

    df_centrifuge, acc_dict = read_centrifuge_results(centrifuge_file, headers)
    #print("df_centrifuge",df_centrifuge)

    taxonomy_df = parse_taxonomy_output(acc_dict)

    # Combine original centrifuge results with taxonomy results
    final_df = pd.concat([df_centrifuge, taxonomy_df], axis=1)
    final_df.to_csv(args.output, sep='\t', index=True)

    # Count occurrences of species where domain is 'Eukaryota'
    counts_e = final_df[final_df['centrifuge domain'] == 'Eukaryota']['centrifuge species'].value_counts()

    # Sort the counts in descending order
    sorted_counts = counts_e.sort_values(ascending=False)

    # Filter species with counts > 10
    filtered_counts = sorted_counts[sorted_counts > 10]
    # Convert the Series to a DataFrame
    filtered_counts_df = filtered_counts.reset_index()

    # Rename the columns for clarity (optional)
    filtered_counts_df.columns = ['species', 'centrifuge_count']

    # Print the resulting DataFrame without index
    print("\nEukaryotic species with more than 10 contigs detected by Centrifuge:\n")
    print(filtered_counts_df)



if __name__ == '__main__':
    translate_taxid_superkingdom()
