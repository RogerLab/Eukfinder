"""
Parse results from various bioinformatics tools and generate:
1. Recovered eukaryotic nuclear genome FASTA file
2. Recovered mitochondrial genome FASTA file

Version 2.0
Author: Dandan Zhao
Date: Dec 2024

Usage:
    python3 final_parser.py -i <input_fasta> [-c <centrifuge_results>] [-d <depth_file>] \
                            [-b <binning_results>] [-p <plast_results>]

This script performs the following tasks:
1. Parses results from Centrifuge, Plast, Metaxa2, and BLAST (Mitochondrial).
2. Filters contigs based on user-defined thresholds and combines parsed results.
3. Writes two output FASTA files for recovered eukaryotic nuclear and mitochondrial genomes.
"""

import glob
import os, sys
import re
import argparse
import pandas as pd
from Bio import Entrez, SeqIO
from collections import Counter

pd.options.mode.chained_assignment = None  # Suppress the warning

def parse_fasta(fasta_file):
    """
    Parses the input FASTA file.

    Args:
        fasta_file (str): Path to the input FASTA file.

    Returns:
        iterator: FASTA sequence iterator using BioPython.
    """
    return SeqIO.parse(open(fasta_file), 'fasta')

def parse_coverage(depth_file):
    """
    Parses the depth file for contig coverage and length.

    Args:
        depth_file (str): Path to the depth file.

    Returns:
        pd.DataFrame: DataFrame with columns ['contigLen', 'totalAvgDepth'] indexed by 'contigName'.
    """
    df = pd.read_csv(depth_file, sep='\t')
    df = df.loc[:, ['contigName', 'contigLen', 'totalAvgDepth']]
    df['contigLen'] = df['contigLen'].astype(int)
    df.set_index('contigName', inplace=True)
    return df


def count_contigs_euk(df, domain_column, species_column, file):
    """
    Counts occurrences of Eukaryotic species with more than 10 contigs.

    Args:
        df (pd.DataFrame): DataFrame containing taxonomy results.
        domain_column (str): Column name indicating the domain (e.g., 'centrifuge domain').
        species_column (str): Column name indicating the species.
        file (str): Label for the count column.

    Returns:
        pd.DataFrame: DataFrame containing species and their counts.
    """
    counts_e = df[df[domain_column] == 'Eukaryota'][species_column].value_counts()
    filtered_counts = counts_e[counts_e > 10].sort_values(ascending=False).reset_index()
    filtered_counts.columns = ['species', f'{file}_count']
    print(filtered_counts)
    return filtered_counts


# Function to read parsed Centrifuge results and process it
def read_Centrifuge(infile):
    """
    Reads parsed Centrifuge results.
    Parameters:
        infile: File path to parsed Centrifuge result file.
    Returns:
        Processed DataFrame with 'readID' and 'taxID', 'SuperKingdom'columns.
    """
    df = pd.read_csv(infile, sep='\t', header=0)


    # Set 'readID' as the index
    df.index = df['readID']
    df = df.drop('readID', axis=1)
    return df

# Function to read parsed binning results
def read_binning(infile):
    """
    Reads parsed binning results.
    Parameters:
        infile: File path to parsed binning result file.
    Returns:
        Processed DataFrame with columns: 'seqID', 'MaxBin', 'metabat', 'MyCC_4mer', 'MyCC_5mer', 'MyCC_56mer' .
    """
    df = pd.read_csv(infile, sep='\t', header=0)

    # Set 'seqID' as the index
    df.index = df['seqID']
    df = df.drop('seqID', axis=1)
    return df

# Function to read parsed Plast results
def read_Plast(infile):
    df = pd.read_csv(infile, sep='\t', header=0)

    # Set 'query ID' as the index
    df.index = df['query ID']
    df = df.drop('query ID', axis=1)
    return df

def parse_Metaxa2_taxonomy(infile):
    #print("\tNow begin to parse Metaxa2_taxonomy file: %s" % infile)
    df_Metaxa2 = pd.read_csv(infile, sep='\t', header=None)
    df_Metaxa2_cols = 'qseqid classification pident length score'.strip().split(' ')
    df_Metaxa2.columns = df_Metaxa2_cols
    df_Metaxa2[['pident','length','score' ]]=df_Metaxa2[['pident','length','score' ]].apply(pd.to_numeric)
    df_Metaxa2.index = df_Metaxa2['qseqid']
    df_Metaxa2 = df_Metaxa2.drop('qseqid', axis=1)
    # filter for results with pident >= 90.0 and length >= 300
    df_filtered = df_Metaxa2[(df_Metaxa2['pident'] >= 90.0) & (df_Metaxa2['length'] >= 300)]
    #print(df_filtered)
    return df_filtered


def parse_BLAST_Mit(infile):
    #print("Now begin to parse BLAST_Mit file: %s" % infile)
    df = pd.read_csv(infile, sep='\t', header=None)
    df_cols = 'qseqid sseqid_Mit stitle_Mit evalue_Mit pident_Mit qcovhsp nident mismatch length_Mit slen qlen qstart qend sstart send staxids sscinames sskingdoms'.strip().split(' ')
    df.columns = df_cols
    df = df.loc[:,['qseqid','sseqid_Mit','stitle_Mit', 'evalue_Mit',  'pident_Mit', 'length_Mit', 'staxids', 'sscinames', 'sskingdoms']]

    #Group results by adding hit length
    df = df.groupby('qseqid', as_index=False).first()

    df[['pident_Mit','length_Mit' ]]=df[['pident_Mit','length_Mit' ]].apply(pd.to_numeric)
    df = df[(df['sskingdoms'] != 'Bacteria') & (df['sskingdoms'] != 'Archaea')]
    df['sseqid_Mit'] = df['sseqid_Mit'].astype(str)
    df_filtered = df[(df['pident_Mit'] >= 80.0) & (df['length_Mit'] >= 100)]
    #df_filtered.loc[:, 'organelle'] = 'Mitochondrion'

    df_filtered = df_filtered.copy()  # Make an explicit copy to avoid ambiguity
    df_filtered.loc[:, 'organelle'] = 'Mitochondrion'

    list_Mito_seqid = df_filtered['qseqid'].tolist()
    df_filtered.index = df_filtered['qseqid']
    df_filtered =  df_filtered.drop('qseqid', axis=1)

    return list_Mito_seqid, df_filtered

# Define a function to calculate organized counts for 'E' and 'B' in SuperKingdom
def calculate_and_organize_counts(df, column, domain_column):
    total_counts = df[column].value_counts().reset_index()
    total_counts.columns = [column, 'Total number']
    
    counts_e = df[df[domain_column] == 'Eukaryota'][column].value_counts().reset_index()
    counts_e.columns = [column, 'Count_E']
    
    counts_b = df[df[domain_column] == 'Bacteria'][column].value_counts().reset_index()
    counts_b.columns = [column, 'Count_B']

    counts_a = df[df[domain_column] == 'Archaea'][column].value_counts().reset_index()
    counts_a.columns = [column, 'Count_A']

    
    merged = pd.merge(total_counts, counts_e, on=column, how='left').fillna(0)
    merged = pd.merge(merged, counts_b, on=column, how='left').fillna(0)
    merged = pd.merge(merged, counts_a, on=column, how='left').fillna(0)

    # Convert specified columns to integer type
    merged[['Count_E', 'Count_B', 'Count_A']] = merged[['Count_E', 'Count_B', 'Count_A']].astype(int)
    
    merged['%E'] = (merged['Count_E'] / merged['Total number']).round(2)
    merged['%B'] = (merged['Count_B'] / merged['Total number']).round(2)
    merged['%A'] = (merged['Count_A'] / merged['Total number']).round(2)
    
    return merged

def detect_species(dict):
    # Flatten the dictionary values into a single list of strings
    all_strings = [item for sublist in dict.values() for item in sublist]

    # Split each string into words and collect all words
    all_words = [word for string in all_strings for word in string.split()]

    # Count the occurrences of each word
    word_counts = Counter(all_words)

    # Find the most common word (or specify filtering criteria if needed)
    species_output = max(word_counts, key=word_counts.get)

    #print(f"Species output: {species_output}")
    return species_output

def write_Mit_fasta(df, list, file, fasta_seq):
    basename = file[:file.index(".fasta")]  
    # Iterate over each species and find corresponding contigNames
    for item in list:
        contigs = df[df['stitle_Mit'].str.startswith(item)]['contigName'].tolist()

        fasta_Mit_output = "%s_Mito_%s.fas" %(basename, item)
        print(f"Write {len(contigs)} contings to {item} Mito genome file: {fasta_Mit_output}. ")

        # Re-parse the FASTA file to avoid iterator exhaustion
        fasta_seq = parse_fasta(file)

        with open(fasta_Mit_output, "w") as out:
            for seq in fasta_seq:
                if seq.id in contigs:
                    SeqIO.write([seq], out, "fasta")


def Main():
    """
    Main function to process input files and generate FASTA outputs.
    """

    parser = argparse.ArgumentParser(description="Process genome results and generate FASTA files.")
    parser.add_argument("-i", "--input", required=False, help="Input FASTA file", default="Eukfinder_long.fasta")
    parser.add_argument("-c", "--centrifuge", required=False, default="parsed_centrifuge_results_eukLong.txt")
    parser.add_argument("-d", "--depth", required=False, default="Eukfinder_long_EUnk.depth.txt")
    parser.add_argument("-b", "--binning", required=False, default="binning_results.tsv")
    parser.add_argument("-p", "--plast", required=False, default="parsed_Plast_Acc2tax_results.txt")
    parser.add_argument("-m", "--mito", required=False, default="Eukfinder_long_BLAST4Mit.out")


    args = parser.parse_args()
    file = args.input
    species_dict = {}
    basename = file[:file.index(".fasta")]  

    # Use Biopython SeqIO to read in contigs
    fasta_seq = parse_fasta(file)

    # Read Parseed Centrifuge result file
    df_centrifuge = read_Centrifuge(args.centrifuge)
    # Print the dimensions of the DataFrame
    # print(f"df_centrifuge DataFrame dimensions: {df_centrifuge.shape}")

    count_centrifuge = count_contigs_euk(df_centrifuge, 'centrifuge domain', 'centrifuge species', 'centrifuge')
    species_dict['centrifuge'] = count_centrifuge['species'].tolist()

    # parsed Plast result file
    df_Plast = read_Plast(args.plast)
    # Print the dimensions of the DataFrame
    # print(f"df_Plast DataFrame dimensions: {df_Plast.shape}")

    count_Plast = count_contigs_euk(df_Plast, 'Plast domain', 'Plast species', 'Plast')
    species_dict['Plast'] = count_Plast['species'].tolist()
    #print(species_dict)
    

    # Parse depth coverage file
    df_depth = parse_coverage(args.depth)

    # Print the dimensions of the DataFrame
    print(f"df_depth DataFrame dimensions: {df_depth.shape}")


    # parsed binning result file
    df_binning = read_binning(args.binning)
    # Print the dimensions of the DataFrame
    # print(f"df_binning DataFrame dimensions: {df_binning.shape}")



    # Parse Metaxa2 results:
    LSU_file = glob.glob("Metaxa2_results/*_metaxa2_LSU.taxonomy.txt")
    SSU_file = glob.glob("Metaxa2_results/*_metaxa2_SSU.taxonomy.txt" )
    df_LSU_filtered = parse_Metaxa2_taxonomy(LSU_file[0])
    df_LSU_filtered.columns=['tax_LSU','ident_LSU','length_LSU','score_LSU']
    df_LSU_filtered['Metaxa2 Mito'] = df_LSU_filtered['tax_LSU'].apply(lambda x: x.split(';')[0] if pd.notna(x) and x.startswith('Mitochondria') else None)

    df_LSU_filtered['Metaxa2 species LSU'] = df_LSU_filtered["tax_LSU"].apply(lambda x: x.split(';')[-2] if pd.notna(x) else None)
    #print(df_LSU_filtered)

    species_dict['LSU'] = df_LSU_filtered['Metaxa2 species LSU'].tolist()

    df_SSU_filtered = parse_Metaxa2_taxonomy(SSU_file[0])
    df_SSU_filtered.columns=['tax_SSU','ident_SSU','length_SSU','score_SSU']
    df_SSU_filtered["Metaxa2 species SSU"] = df_SSU_filtered["tax_SSU"].apply(lambda x: x.split(";")[-2] if pd.notna(x) else None)

    #print(df_SSU_filtered)
    species_dict['SSU'] = df_SSU_filtered['Metaxa2 species SSU'].tolist()
    print("\nspecies_dict after parsing Metaxa2 results:")
    print(species_dict)


    metaxa2_out = pd.concat([df_LSU_filtered,df_SSU_filtered], axis =1, sort=False)

    # Parse BLAST_Mitochondrial results:
    Mito_contigs, df_Mito = parse_BLAST_Mit(args.mito)

    df_combined = pd.concat([df_depth, df_centrifuge, df_Plast, df_binning, df_Mito, metaxa2_out], axis=1)
    #print(df_combined.head())

    df_combined = df_combined.rename_axis('contigName').reset_index()
    df_combined.set_index('contigName', drop=False)


    species_detected = detect_species(species_dict)
    print("\nspecies_detected:\t", species_detected)

    # Output results
    output_file = "Combined_binning_results_CPB_Mito.tsv"
    print(f"Write results to {output_file}.")
    df_combined.to_csv(output_file, sep='\t', index_label="contig")
    #print("Finished.")


    # Ensure 'Metaxa2 species SSU' exists in the DataFrame
    if 'Metaxa2 species SSU' in df_combined.columns:
        # Ensure the column is of type string
        df_combined['Metaxa2 species SSU'] = df_combined['Metaxa2 species SSU'].astype(str)
    
        # Filter rows where 'Metaxa2 species SSU' contains the species_detected
        filtered_rows = df_combined[df_combined['Metaxa2 species SSU'].str.contains(species_detected, na=False)]
       
    
        # Extract corresponding values in 'totalAvgDepth'
        total_avg_depth_SSU = filtered_rows['totalAvgDepth'].tolist()[0]
    
        # Print the results
        print("Depth for SSU contig:\t",total_avg_depth_SSU)
    
        # Remove contigs with larger depth than SSU contig
        df_combined = df_combined[df_combined['totalAvgDepth'] <= total_avg_depth_SSU ]

    else:
        print("Column 'Metaxa2 species SSU' does not exist in the DataFrame.")



    columns = ['Bacteria', 'Archaea']

    # Remove rows based on conditions
    df_combined = df_combined[~(    (df_combined['ident'] > 0.9) & 
                                    (df_combined['align_length'] >= 1000) & 
                                    (df_combined['Plast domain'].isin(columns)))]



    # Filter DataFrame for SSU contigs

    filtered_df_SSU = df_combined[(df_combined['ident_SSU'] > 97.0) & (df_combined['length_SSU'] > 500)]
    if not filtered_df_SSU.empty:
        contigs_SSU = filtered_df_SSU['contigName'].tolist()[0]

    # Filter DataFrame for Mito contigs
    filtered_df_Mit = df_combined[(df_combined['pident_Mit'] > 97.0) & (df_combined['length_Mit'] > 500) ]

    # Extract unique species from the first word of 'stitle_Mit'
    list_species = filtered_df_Mit['stitle_Mit'].apply(lambda x: x.split()[0]).unique().tolist()



    # Iterate over each species and find corresponding contigNames
    write_Mit_fasta(filtered_df_Mit, list_species, file, fasta_seq)


    # Remove Mitochondrial rows based on Blast vs Mitochondrial database 
    df_combined = df_combined[~(    (df_combined['pident_Mit'] > 99.0) & 
                                    (df_combined['length_Mit'] >= 1000) & 
                                    (df_combined['sskingdoms'] == 'Eukaryota'))]

    print(list_species)

    # *** If there are more than one species in list_species, put the following lines into a for loop to iterate the items in list_species

    # Add a new column list_species[0] based on conditions
    df_combined[list_species[0]] = df_combined.apply(
        lambda row: list_species[0] if (
            list_species[0] in str(row['Plast species']) or 
            list_species[0] in str(row['centrifuge species'])
        ) else None, axis=1
    )

    df_filter_species = df_combined[df_combined[list_species[0]] == list_species[0]]

    # Add the 'Domain' column based on the specified conditions
    df_combined['Domain'] = df_combined.apply(
        lambda row: "Eukaryota" if row[list_species[0]] == list_species[0] else 
                    (row['centrifuge domain'] if pd.notna(row['centrifuge domain']) else row['Plast domain']),
        axis=1
    )

    # Calculate results for each column
    columns_to_analyze = ['MyCC_5mer', 'MyCC_56mer', 'MyCC_4mer']

    results = {col: calculate_and_organize_counts(df_combined, col, 'Domain') for col in columns_to_analyze}

    # Initialize an empty list
    MyCC_euk_bins = []
    # Display results
    for col, result in results.items():
        #print(f"Results for {col}:")
        #print(result)
        #print("-")
        filtered_df = result[result['%E'] > 0.50]
        MyCC_euk_bins.extend(filtered_df[col].tolist())

    print("MyCC_euk_bins:\t", MyCC_euk_bins)


    # Count the appearance of values in the specified columns that are in MyCC_euk_bins
    df_combined['MyCC_count'] = df_combined[['MyCC_5mer', 'MyCC_56mer', 'MyCC_4mer']].apply(
        lambda row: sum(val in MyCC_euk_bins for val in row), axis=1
    )

    # Print the value_counts() for the MyCC_count column
    print(df_combined["MyCC_count"].value_counts())


    # Filter the DataFrame for rows where 'MyCC_count' >= 2
    df_filter_MyCC = df_combined[df_combined["MyCC_count"] >= 2]

    # Filter the DataFrame for rows where 'centrifuge species' is targeted species
    df_filter_centrifuge = df_combined[df_combined['centrifuge species'] == list_species[0]]

    # Filter the DataFrame for rows where 'Plast species' is targeted species
    df_filter_Plast = df_combined[df_combined['Plast species'] == list_species[0]]

    df_filtered = pd.concat([df_filter_MyCC,df_filter_centrifuge, df_filter_Plast])
    df_filtered = df_filtered.drop_duplicates(subset='contigName')

    contigs = df_filtered['contigName'].tolist()

    fasta_output = "%s_%s.fas" %(basename, list_species[0])

    print(f"Write {len(contigs)} contings to {list_species[0]} genome file: {fasta_output}. ")

    # Re-parse the FASTA file to avoid iterator exhaustion
    fasta_seq = parse_fasta(file)

    with open(fasta_output, "w") as output:
        for seq in fasta_seq:
            if seq.id in contigs:
                SeqIO.write([seq], output, "fasta")


if __name__ == '__main__':

    Main()





