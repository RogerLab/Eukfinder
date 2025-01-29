import sys, os
from subprocess import PIPE, Popen

# Step 1: Read the assembly summary file once to extract both FTP link and TaxID
def readingCatalog_combined(infile):
    """
    Read the assembly summary file only once to extract both FTP link and TaxID.
    
    :param infile: Path to the NCBI assembly summary file
    :return: A dictionary with genome accession as the key and a tuple (FTP link, TaxID) as the value
    """
    with open(infile) as In:
        combined_data = {}
        for line in In:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                acc = line[0].split('.')[0]  # genome accession
                ftp_link = line[19]  # column 20 for FTP link
                tax_id = line[6]  # column 7 for TaxID
                combined_data[acc] = (ftp_link, tax_id)
        return combined_data


# Step 2: Read the genome IDs that you want to recover
def readgenomeIDs(infile):
    with open(infile) as I:
        # Create a set to store the extracted genome accessions
        genome_ids = set()
        
        for line in I:
            # Extract genome accession by splitting at the dot
            acc = line.split('.')[0]  # Extract accession part before '.'
            genome_ids.add(acc.strip())  # Add the accession to the set, removing any trailing spaces

        return genome_ids


# Step 4: Perform the genome downloading based on genome accessions
def perform_combined(combined_data, IDs, list_downloaded_GenomeID):
    """
    Downloads genomes based on genome accessions and retrieves both FTP and TaxID in one step.
    
    :param combined_data: A dictionary with genome accession as the key and a tuple (FTP link, TaxID) as the value
    :param IDs: A set of genome accessions to be downloaded and processed
    :param list_downloaded_GenomeID: A list of already downloaded GenomeIDs
    """
    for id in IDs:
        if id in list_downloaded_GenomeID:
            print(f'{id} already downloaded. Skipping...')
            continue  # Skip if already downloaded
        
        print('About to retrieve %s' % id)
        if id in combined_data.keys():
            ftp_link, tax_id = combined_data[id]
            basename = ftp_link.split('/')[-1]
            #print('Downloading from: %s' % ftp_link)
            ftpcmd = f'wget {ftp_link}/{basename}_genomic.fna.gz'
            try:
                retrieve = Popen(ftpcmd, stderr=PIPE, stdout=PIPE, shell=True)
                o, e = retrieve.communicate()
                os.system(f'gunzip {basename}_genomic.fna.gz')
            except:
                print('Failed to execute: %s' % ftpcmd)
        else:
            print('Failed to find: %s' % id)

    return


# Step 5: Create genome list by reading filenames
def Reading_Genome_filenames(infilelist):
    with open(infilelist) as In:
        L = [f.strip() for f in In if f.strip() != '']
    dict = {}
    for line in L:
        file_name = line.split('_')
        GenomeID = '%s_%s' % (file_name[0], file_name[1])
        GenomeID =GenomeID.split('.')[0]
        dict[GenomeID] = line
    print("Finish Reading_Genome_filenames\n")
    return dict


# Step 6: Parse the headers in the sequence files
def reading_headers(infile):
    list_headers = []
    with open(infile) as I:
        I = I.read().split('>')[1:]
        for read in I:
            entry = read.split('\n')
            key = entry[0].split(' ')[0]
            if not key in list_headers:
                list_headers.append(key)
        return list_headers


# Step 8: Write the Genome to TaxID map output
def writing_output_map(map, infile):
    file_out = open(infile.replace(".txt", "_genome2taxid.txt"), 'w')
    print("Writing genome2taxid dict to file.\n")
    for key, value in map.items():
        line = '%s\t%s\n' % (key, value)
        file_out.write(line)
    file_out.close()


if __name__ == '__main__':
    '''
    Usage: python combined_script.py <assembly_summary_file> <genomeID_file>
    '''
    # Step 1: Read the assembly summary file once
    combined_data = readingCatalog_combined(sys.argv[1])

    # Step 2: Read the genome IDs
    Genomes2Recover = readgenomeIDs(sys.argv[2])
    print("Step 2, %s Genomes to Recover:\n", len(Genomes2Recover))

    # === Step 3.1: Get list of downloaded .fna files ===
    list_downloaded_files = [f for f in os.listdir() if f.endswith('.fna')]
    list_downloaded_GenomeID = []

    # === Step 3.2: Parse filenames to extract GenomeID ===
    for file in list_downloaded_files:
        file_name = file.split('_')
        GenomeID = '%s_%s' % (file_name[0], file_name[1])
        GenomeID = GenomeID.split('.')[0]
        list_downloaded_GenomeID.append(GenomeID)

    # Step 3.3: Print out the downloaded GenomeIDs
    print("Downloaded Genomes: %s\n", len(list_downloaded_GenomeID))

    # Step 4: Perform the download and mapping
    print("In total %s genomes to be downloaded." % len(Genomes2Recover))
    perform_combined(combined_data, Genomes2Recover, list_downloaded_GenomeID)

    # Step 5: Create genome list (simulate the shell command `ls -lthr *.fna | awk '{print $NF}' > genome_list.txt`)
    os.system("ls -lthr *.fna | awk '{print $NF}' > genome_list.txt")

    # Step 6: Map Genome Accession to TaxID (use same combined data dictionary)
    genome_dic = Reading_Genome_filenames("genome_list.txt")
    map = {}
    for genomeID, genomefile in genome_dic.items():
        tax_id = combined_data[genomeID][-1]
        map[genomefile] = tax_id

        # Step 7: Process each genome file, map sequence headers to TaxID
        list_headers = reading_headers(genomefile)
        with open(genomefile.replace(".fna", "_seq2taxid_map.txt"), 'w') as file_out:
            for item in list_headers:
                file_out.write(f'{item}\t{tax_id}\n')

    # Step 8: Write the output map to a file
    writing_output_map(map, sys.argv[2])

    print("Finished.")
