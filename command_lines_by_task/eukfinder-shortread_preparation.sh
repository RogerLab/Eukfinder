#!/bin/bash

# Activate Conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate Eukfinder

# Input data (Update these paths based on your setup)
R1=/path/to/your/R.1.fastq
R2=/path/to/your/R.2.fastq
host_genome=/path/to/your/host_genome.fasta

# Centrifuge Database
# If the folder is "Centrifuge_DB" and contains files:
# Centrifuge_DB.1.cf, Centrifuge_DB.2.cf, Centrifuge_DB.3.cf, Centrifuge_DB.4.cf
# Then specify the location for the Centrifuge database as:
# centrifuge=/path/to/your/Centrifuge_DB/Centrifuge_DB
# Note: The software requires the base path without ".1.cf"
centrifuge=/path/to/your/Centrifuge_DB/Centrifuge_DB

# Adapter file for trimming
adapters=TrueSeq2_NexteraSE-PE.fa

# Create output directory if it doesn't exist
output_dir="read_prep"
mkdir -p $output_dir

# Command line for read preparation
Eukfinder.py read_prep \
    --r1 $R1 \
    --r2 $R2 \
    -n 20 \
    --hcrop 5 \
    -l 10 \
    -t 7 \
    --wsize 40 \
    --qscore 25 \
    --mlen 30 \
    --mhlen 30 \
    --hg $host_genome \
    -o $output_dir \
    --cdb $centrifuge \
    -i $adapters

# Deactivate Conda environment
conda deactivate
