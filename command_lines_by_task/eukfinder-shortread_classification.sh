#!/bin/bash

# Activate Conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate Eukfinder

# Input data
plastdb=Plast_DB.fasta            # Path to the PLAST database
plastmap=Plast_DB.map             # Path to the PLAST map file
centrifuge=/path/to/Centrifuge_DB # Path to the Centrifuge database (base name without .1.cf, .2.cf, etc.)
acc2tax=/path/to/Acc2Tax_DB       # Path to the Accession-to-Taxonomy database

# Previously prepared reads
R1=read_prep_p.1.fastq
R2=read_prep_p.2.fastq
uR1R2=read_prep_un.fastq

# Centrifuge read classification output files
p_class=read_prep_centrifuge_P
up_class=read_prep_centrifuge_UP

# Run Eukfinder classification
eukfinder.py short_seqs \
    --r1 $R1 \
    --r2 $R2 \
    --un $uR1R2 \
    -n 20 \
    -z 6 \
    -t False \
    --max_mem 100 \
    -p $plastdb \
    -m $plastmap \
    -a $acc2tax \
    -e 0.01 \
    --pid 60 \
    --cov 30 \
    --pclass $p_class \
    --uclass $up_class \
    -o first_round \
    --cdb $centrifuge \
    --mhlen 30 \
    -k 21,33,55

# Deactivate Conda environment
conda deactivate
