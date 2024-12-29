# Building a Mitochondrial-Related Genome Database for Metagenomic Analysis

This database is designed to detect mitochondrial sequences in metagenomic assemblies by including mitochondrial genomes and closely related bacterial and archaeal genomes. The focus is on high-quality reference genomes at the complete assembly level. Below are the detailed steps for constructing this database:

## 1. Inclusion Criteria

The database includes:

Mitochondrial Genomes:
Specifically, the genome file mitochondrion.1.1.genomic.fna.

Bacterial Genomes:

Classes included:
- Alphaproteobacteria (closest relatives to mitochondria).
- Deltaproteobacteria.
- Gammaproteobacteria.
 
Only genomes classified as Reference Genomes at the Complete Assembly level are included.

Archaeal Genomes:

Phyla and classes included:

- Lokiarchaeota.
- Thorarchaeota.
- Odinarchaeota.
- Methanobacteria.
- Halobacteria.
- Sulfolobus.
- Thermoproteus.
 
Only Reference Genomes at the Complete Assembly level are included.


## 2. Building the Database

###Step 1: Collect Mitochondrial Genome

Download the mitochondrial genome (mitochondrion.1.1.genomic.fna) from NCBI:

Source: NCBI FTP or Nucleotide Database.

Command for direct download:

```sh
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.1.1.genomic.fna.gz
gunzip mitochondrion.1.1.genomic.fna.gz
```

###Step 2: Identify Closely Related Bacteria and Archaea Genomes

To include bacterial and archaeal genomes closely related to mitochondria, perform the following:

Access the NCBI Nucleotide Database:

Visit NCBI Nucleotide.
Set the Search Criteria to “Genome” to focus on complete genome sequences.
Refine Search for Target Organisms:

Search for each bacterial or archaeal group (e.g., “Alphaproteobacteria,” “Lokiarchaeota”).

Apply the following filters:

- Assembly Level: Complete Genome.
- Genome Type: Reference Genomes only.
- Release Date: Optional, to ensure recent high-quality data.
- Download Genome Metadata:

On the NCBI search results page:

Click the box on the first row to select all the genomes;

Select Column “RefSeq” or "GenBank"

Choose "Download Table" to save the list of genome accession numbers as genome_list.txt.


### Step 3: Download Genomes from NCBI

Once the list of accession numbers (genome_accession_list.txt) has been prepared, download the corresponding genome files from NCBI’s FTP server:

#### Recommended Method: 

Use this python script downloading both the sequence fna file and the map file.
[genome_download_map.py](https://github.com/RogerLab/Eukfinder/blob/main/Building_custom_DB/genome_download_map.py)
See [Using a Custom Python Script](https://github.com/RogerLab/Eukfinder/tree/main/Building_custom_DB#method-3-using-a-custom-python-script) for more detailed instructions.

```sh
python3 genome_download_map.py assembly_summary_genbank.txt genome_list.txt
cat *.fna > genome.fasta
cat *genome2taxid.txt > genome2taxid.map
```

#### Alternative Method:
Use NCBI Entrez Direct (Command-Line Method):

Install Entrez Direct if not already available:

```sh
sudo apt-get install -y entrez-direct
```

Run the following command to download all genomes:

```sh
cat genome_list.txt | xargs -n 1 -I {} efetch -db nuccore -id {} -format fasta > bacterial_archaeal_genomes.fna
```

Explanation:
genome_list.txt: File containing one genome accession number per line.
efetch: Retrieves genome sequences in FASTA format.
bacterial_archaeal_genomes.fna: Combined output file containing all bacterial and archaeal genomes.


### Step 4: Combine All Genomes into a Single Database

Merge the mitochondrial genome and the bacterial/archaeal genomes into one file:

```sh
cat mitochondrion.1.1.genomic.fna bacterial_archaeal_genomes.fna > mitochondrial_database.fna
```

### Step 5: Make a BLAST Database

After combining all mitochondrial, bacterial, and archaeal genomes into a single FASTA file, you can create a BLAST database using the makeblastdb command.

```sh
makeblastdb -in mitochondrial_database.fna -dbtype nucl -parse_seqids -taxid_map genome2taxid.map -out mito_blast_db
```

## 3. Final Output
The final database file mitochondrial_database.fna contains:

The mitochondrial genome (mitochondrion.1.1.genomic.fna).
Bacterial genomes from Alphaproteobacteria, Deltaproteobacteria, and Gammaproteobacteria.
Archaeal genomes from Lokiarchaeota, Thorarchaeota, Odinarchaeota, Methanobacteria, Halobacteria, Sulfolobus, and Thermoproteus.
All sequences included are high-quality Reference Genomes at the Complete Assembly level.
