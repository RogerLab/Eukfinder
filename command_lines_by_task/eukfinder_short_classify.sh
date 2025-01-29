#!/bin/bash
# -------------------------------------------------------------------
# eukfinder_short_classify.sh
#
# Description:
#   This script classifies short reads using Eukfinder. It automatically
#   activates the Eukfinder conda environment, runs the classification,
#   and then deactivates the environment.
#
# Instructions:
#   1. Update the input file paths (R1, R2, uR1R2, etc.) to match your data.
#   2. Update the database paths as needed.
#   3. Run: ./eukfinder_short_classify.sh
#
# Notes:
#   - For a first-time run or when updating taxonomy, set "-t True" in eukfinder.
#   - See `eukfinder short_seqs -h` for additional flags.
# -------------------------------------------------------------------

# 1. Define input files and output prefix
R1="read_prep_p.1.fastq"
R2="read_prep_p.2.fastq"
uR1R2="read_prep_un.fastq"
p_class="read_prep_centrifuge_P"
up_class="read_prep_centrifuge_UP"
prefix="Test_short"

# 2. Full paths to databases (edit these as needed)
plastdb="/path/to/Plast_DB.fasta"       # Path to the PLAST database
plastmap="/path/to/Plast_DB.map"        # Path to the PLAST map file
centrifuge="/path/to/Centrifuge_DB"     # Base name (no .1.cf, .2.cf, etc.)
acc2tax="/path/to/Acc2Tax_DB"           # Accession-to-Taxonomy DB

# 3. Activate the Eukfinder Conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate eukfinder

# 4. Run Eukfinder classification (short reads)
eukfinder short_seqs \
    --r1 "$R1" \
    --r2 "$R2" \
    --un "$uR1R2" \
    -n 20 \
    -z 6 \
    -t False \
    --max_mem 100 \
    -p "$plastdb" \
    -m "$plastmap" \
    -a "$acc2tax" \
    -e 0.01 \
    --pid 60 \
    --cov 30 \
    --pclass "$p_class" \
    --uclass "$up_class" \
    -o "$prefix" \
    --cdb "$centrifuge" \
    --mhlen 30 \
    -k 21,33,55

# 5. Deactivate Conda environment
conda deactivate
