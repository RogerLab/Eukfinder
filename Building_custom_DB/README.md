# Building Custom Databases

## Overview

In metagenomic and genomic analyses, the ability to accurately classify sequences—whether contigs or raw reads—is crucial for understanding the composition of complex microbial and eukaryotic communities. Eukfinder leverages classification tools like Centrifuge to determine which sequences originate from eukaryotic organisms. To perform these classifications, Centrifuge relies on specialized databases that link each sequence to its corresponding taxonomy. By default, users may rely on pre-built databases; however, constructing custom Centrifuge databases allows for greater precision and adaptability to a project’s unique requirements.

Building a custom database ensures that you target the exact taxa, strains, or groups of interest, streamline your analyses, and improve classification accuracy. Whether you need a broad reference set (e.g., all archaeal, bacterial, and viral genomes), a focused reference for a particular genus or species, or a hand-curated set of genomes from a given clade or timeframe, the methods described here provide multiple approaches to assembling the necessary data and indexing it for Centrifuge.

## Building Centrifuge Databases

This document outlines three primary strategies for creating custom Centrifuge databases:

1. **Using Centrifuge’s Built-In Tools**
Leverage the centrifuge-download utility to directly fetch and assemble complete datasets (like archaeal, bacterial, and viral references or the NCBI nt database) from official sources. This method is ideal when working with well-established, large-scale reference sets.

2. **Using ncbi-genome-download**
Target specific organisms, taxa, or genera by employing the ncbi-genome-download command-line tool. This approach lets you define your dataset more narrowly, ensuring that only genomes fitting your scientific question are included.

3. **Using a Custom Python Script and Assembly Summaries**
Gain full control over database composition by manually selecting accession numbers from NCBI, applying filtering criteria (like assembly level or release date), and downloading sequences using a Python script. This method grants maximum customization for specialized studies that require curated sets of genomes.

For each approach, you’ll generate or retrieve:

- **Genomic FASTA files**: The raw nucleotide sequences that form the building blocks of your database.
- **Taxonomy files**: nodes.dmp and names.dmp, which are essential for linking sequences to their respective taxonomic identifiers.
- **Mapping files**: Sequence-to-TaxID maps that align each genome or sequence header to a taxonomic ID, enabling Centrifuge to place sequences within the correct taxonomic branch.
  
After constructing the database, you’ll run centrifuge-build to create the index files that Centrifuge needs for rapid classification. Once indexed, these databases can be integrated into Eukfinder’s workflow, enabling refined eukaryotic sequence detection and improved downstream analyses such as eukaryotic binning, phylogenetic studies, and ecological assessments.

By following the steps in this guide, you’ll have the flexibility to customize your reference databases, improve classification accuracy, and adapt Eukfinder to the evolving demands of your research projects.  

- Method 1: Ideal for large, well-established databases like complete bacterial/archaeal/viral sets or nt.
- Method 2: Good for targeted downloads from NCBI based on genus or taxonomic group.
- Method 3: Fully customizable approach using accession numbers and assembly summaries for highly specific datasets.

### Method 1: Using the Centrifuge Manual Instructions

See [Centrifuge manual](https://ccb.jhu.edu/software/centrifuge/manual.shtml#database-download-and-index-building) for details.

#### 1.1 Downloading Genomes Directly with centrifuge-download

You can retrieve reference genomes and taxonomy files directly from NCBI using the centrifuge-download utility. For example, to build a database containing archaeal, bacterial, and viral genomes:

1. Obtain the NCBI taxonomy files:

```sh
centrifuge-download -o taxonomy taxonomy
```

2. Download all complete archaeal, bacterial, and viral genomes and mask low-complexity regions:

```sh
centrifuge-download -o library -m -d "archaea,bacteria,viral" refseq > seqid2taxid.map
```

The seqid2taxid.map file maps each sequence ID to its corresponding taxonomy ID.

3. Merge all downloaded sequences into one file:

```sh
cat library/*/*.fna > input-sequences.fna
```

4. Build the Centrifuge index:

```sh
centrifuge-build -p 4 --conversion-table seqid2taxid.map \
                 --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
                 input-sequences.fna abv
```

#### 1.2 Building an Index from the NCBI nt Database

The NCBI nt database is a comprehensive, non-redundant collection of nucleotide sequences. After downloading the nt FASTA file and the corresponding GI-to-TaxID map, you can construct a Centrifuge index as follows:

1. Download and decompress the nt database:

```sh
wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz
gunzip nt.gz && mv nt nt.fa
```

2. Retrieve the GI-to-TaxID mapping:

```sh
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
gunzip -c gi_taxid_nucl.dmp.gz | sed 's/^/gi|/' > gi_taxid_nucl.map
```

3. Build the Centrifuge index:
Build the Centrifuge index using more threads and a custom --bmax parameter to manage memory:

```sh
centrifuge-build -p 16 --bmax 1342177280 --conversion-table gi_taxid_nucl.map \
                 --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
                 nt.fa nt
```

(See Centrifuge manual for details)


### Method 2: Using ncbi-genome-download

1. Download the desired genomes by specifying a genus or taxonomic group. 

**Prerequisites**: [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download)

For instance, to download *Blastocystis* genomes:

```sh
source activate ncbi-genome-download
ncbi-genome-download --genera Blastocystis -p 4 -r 10 --flat-output --progress-bar --formats fasta,assembly-report protozoa
```

2. Combine all downloaded FASTA files into a single file:

```sh
cat *.fna > genome.fasta
```

3. Generate the sequence-to-taxid mapping file using a custom script:

Example script available here: [Build_Centrifuge_map_from_assembly_report.py](https://github.com/RogerLab/Eukfinder/blob/main/Building_custom_DB/Build_Centrifuge_map_from_assembly_report.py)

This script will extract the organism name, taxid, and a list of sequence IDs from an assembly report file (XXX_assembly_report.txt), and generate a corresponding XXX_genome2taxid.txt for the genome file XXX_genomic.fna. XXX here represent the genome accession number and assemble name, for example: 

**genome file**: GCF_000743755.1_ASM74375v1_genomic.fna               (Download in step 1)

**assembly report**: GCF_000743755.1_ASM74375v1_assembly_report.txt    (Download in step 1)

**map file**: GCF_000743755.1_ASM74375v1_genomic_seq2taxid_map.txt    (output from script: Build_Centrifuge_map_from_assembly_report.py)


```sh
python3 Build_Centrifuge_map_from_assembly_report.py
cat *_genome2taxid.txt > genome2taxid.map
```

**Note**:
This script will also generate a file name "Genome_list.tsv" that list all the genomes included with the species name and taxoID. 

Heres is an example of the file:
| Genome                     | Organism Name                           |  Taxid |
| :--------------------------|:----------------------------------------|:-------|
| GCF_000151665.1_ASM15166v1 | Blastocystis hominis (eukaryotes)       | 12968  |
| GCF_000743755.1_ASM74375v1 | Blastocystis sp. subtype 4 (eukaryotes) | 944170 |



4. Download taxonomy files:

```sh
centrifuge-download -o taxonomy taxonomy
# This creates taxonomy/nodes.dmp and taxonomy/names.dmp
```

5. Build the Centrifuge index:

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

Assume resulted list of geome accession number is in the file **genome_list.txt**.

2. Download the Assembly Summary File:

```sh
# For GenBank:
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt

#F or RefSeq:
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
```

3. Use the provided Python script to download genomes and produce the mapping file:
   
Example script: [genome_download_map.py](https://github.com/RogerLab/Eukfinder/blob/main/Building_custom_DB/genome_download_map.py)

```sh
python3 genome_download_map.py assembly_summary_genbank.txt genome_list.txt
```

```sh
cat *.fna > genome.fasta
cat *genome2taxid.txt > genome2taxid.map
```

4. Build the Centrifuge index:

```sh
centrifuge-build -p 16 --bmax 1342177280 --conversion-table genome2taxid.map \
                 --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
                 genome.fasta genome
```

**Additional Notes**

- Always ensure that sequence headers and mapping files align with the taxonomy files (nodes.dmp and names.dmp).
- Test small datasets first to confirm that the pipeline is set up correctly.
- For more detailed instructions and troubleshooting, consult the official Centrifuge manual and other provided scripts’ documentation.


## Building Plast databases
