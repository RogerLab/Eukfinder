'''
Generate sequenceID to taxid for building centrifuge database
================================================
Usage:
python3 Build_Centrifuge_map_from_assembly_report.py
================================================

Examples of assembly_report file downloaded from NCBI:
GCF_000151665.1_ASM15166v1_assembly_report.txt

Examples of output file:
GCF_000151665.1_ASM15166v1_genome2taxid.txt
'''




import os
import glob

#   Info  #
__author__ = 'Dandan Zhao'
__email__ = 'd.zhao@dal.ca'
__version__ = '2.0.0'
#   End Info   #


# Function to extract information from the assembly report file
def parse_assembly_report(file_path):
    """
    Extract the organism name, taxid, and a list of sequence IDs from an assembly report file.

    Parameters:
    file_path (str): The path to the assembly report file.

    Returns:
    tuple: Contains the organism name (str), taxid (str), and a list of sequence IDs (list of str).
    """
    seq_id_list = []

    with open(file_path, 'r') as file:
        print(f"Processing file: {file_path}")
        for line in file:
            # Extract Organism name from line with "# Organism name:"
            if line.startswith("# Organism name:"):
                organism_name = line.split(":")[1].strip()

            # Extract Taxid from line with "# Taxid:"
            elif line.startswith("# Taxid:"):
                taxid = line.split(":")[1].strip()
            
            # Extract sequence IDs from data lines that do not start with '#'
            elif not line.startswith("#"):
                sequence_id = line.split("\t")[6]
                seq_id_list.append(sequence_id)

    return organism_name, taxid, seq_id_list
      

if __name__ == '__main__':

    # Step 1: Initialize output file for genome list and write the header
    genome_list_output_file = "Genome_list.tsv"
    with open(genome_list_output_file, 'w') as genome_list_out:
        genome_list_out.write("Genome\tOrganism Name\tTaxid\n")

        # Step 2: Read all assembly report files and extract data
        assembly_report_files = glob.glob('*_assembly_report.txt')
        
        for report_file in assembly_report_files:
            # Extract organism name, taxid, and sequence ID list from each report file
            organism_name, taxid, seq_id_list = parse_assembly_report(report_file)

            # Generate the output file name for genome to taxid mapping
            basename = report_file.replace("_assembly_report.txt", "")
            output_file = f"{basename}_genome2taxid.txt"

            # Step 3: Write the sequence ID to taxid mapping for each report file
            with open(output_file, 'w') as seq_id_out:    
                for seq_id in seq_id_list:
                    print(f"Writing Sequence ID: {seq_id} to file: {output_file}")
                    line = f"{seq_id}\t{taxid}\n"
                    seq_id_out.write(line)

            # Step 4: Write the genome list containing all files, organism names, and taxids
            genome_list_line = f"{basename}\t{organism_name}\t{taxid}\n"
            genome_list_out.write(genome_list_line)

    print("All tasks completed successfully.")
