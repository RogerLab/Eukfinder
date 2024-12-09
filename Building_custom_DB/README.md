# Building Custom Databases

## Overview

Here are the instructions on how to use different methods to build custom Centrifuge and Plast databases

## Building Centrifuge Databases

This guide describes several approaches to create custom Centrifuge databases.  

### Method 1: Using the Centrifuge Manual Instructions

See [Centrifuge manual](https://ccb.jhu.edu/software/centrifuge/manual.shtml#database-download-and-index-building) for details.

1.1 Downloading Genomes Directly with centrifuge-download

You can retrieve reference genomes and taxonomy files directly from NCBI using the centrifuge-download utility. For example, to build a database containing archaeal, bacterial, and viral genomes:

Obtain the NCBI taxonomy files:

```sh
centrifuge-download -o taxonomy taxonomy
```

Download all complete archaeal, bacterial, and viral genomes and mask low-complexity regions:

```sh
centrifuge-download -o library -m -d "archaea,bacteria,viral" refseq > seqid2taxid.map
```

The seqid2taxid.map file maps each sequence ID to its corresponding taxonomy ID.

Merge all downloaded sequences into one file:

```sh
cat library/*/*.fna > input-sequences.fna
```

Build the Centrifuge index with multiple threads:

```sh
centrifuge-build -p 4 --conversion-table seqid2taxid.map \
                 --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
                 input-sequences.fna abv
```

1.2 Building an Index from the NCBI nt Database

The NCBI nt database is a comprehensive, non-redundant collection of nucleotide sequences. After downloading the nt FASTA file and the corresponding GI-to-TaxID map, you can construct a Centrifuge index as follows:

Download and decompress the nt database:

```sh
wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz
gunzip nt.gz && mv nt nt.fa
```

Retrieve the GI-to-TaxID mapping:

```sh
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
gunzip -c gi_taxid_nucl.dmp.gz | sed 's/^/gi|/' > gi_taxid_nucl.map
Build the Centrifuge index using more threads and a custom --bmax parameter to manage memory:

```sh
centrifuge-build -p 16 --bmax 1342177280 --conversion-table gi_taxid_nucl.map \
                 --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
                 nt.fa nt
```

(See Centrifuge manual for details)


### Method 2: Using ncbi-genome-download

Download the desired genomes by specifying a genus or taxonomic group. For instance, to download Blastocystis genomes:

```sh
ncbi-genome-download --genera Blastocystis -p 4 -r 10 --flat-output --progress-bar --formats fasta,assembly-report protozoa
```

Combine all downloaded FASTA files into a single file:

```sh
cat *.fna > genome.fasta
```

Generate the sequence-to-taxid mapping file using a custom script:

```sh
# Example script available here:[Build_Centrifuge_map_from_assembly_report.py](https://github.com/RogerLab/Eukfinder/blob/main/Building_custom_DB/Build_Centrifuge_map_from_assembly_report.py)
# https://github.com/RogerLab/Eukfinder/blob/main/Building_custom_DB/Build_Centrifuge_map_from_assembly_report.py
python3 Build_Centrifuge_map_from_assembly_report.py
cat *_genome2taxid.txt > genome2taxid.map
```

Download taxonomy files:

```sh
centrifuge-download -o taxonomy taxonomy
# This creates taxonomy/nodes.dmp and taxonomy/names.dmp
```

Build the Centrifuge index:

```sh
centrifuge-build -p 16 --bmax 1342177280 --conversion-table genome2taxid.map \
                 --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
                 genome.fasta genome
```

### Method 3: Using a Custom Python Script

1. Identify Genome Accession Numbers:
Search the NCBI Nucleotide database (https://www.ncbi.nlm.nih.gov/nuccore) for the organism of interest. 

- Go to [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore)
- Choose "Genome" as search criteria
- Filter results (reference genomes, MAGs, assembly level, release date)
- Download the table and extract accession numbers from the first column	

Assume resulted list of geome accession number is in the file genome_list.txt

2. Download the Assembly Summary File:

```sh
# For GenBank:
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt

#F or RefSeq:
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
```

3. Use the provided Python script to download genomes and produce the mapping file:

```sh
# Example script:[genome_download_map.py](https://github.com/RogerLab/Eukfinder/blob/main/Building_custom_DB/genome_download_map.py)
# https://github.com/RogerLab/Eukfinder/blob/main/Building_custom_DB/genome_download_map.py
python3 genome_download_map.py assembly_summary_genbank.txt genome_list.txt
```

```sh
cat *.fna > genome.fasta
cat *_genome2taxid.txt > genome2taxid.map
```

4. Build the Centrifuge index:

```sh
centrifuge-build -p 16 --bmax 1342177280 --conversion-table genome2taxid.map \
                 --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
                 genome.fasta genome
```

Additional Notes:

1. Always ensure that sequence headers and mapping files align with the taxonomy files (nodes.dmp and names.dmp).
2. Test small datasets first to confirm that the pipeline is set up correctly.
3. For more detailed instructions and troubleshooting, consult the official Centrifuge manual and other provided scripts’ documentation.



## Building Plast databases
