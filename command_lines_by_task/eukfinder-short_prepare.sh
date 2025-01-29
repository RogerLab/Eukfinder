#!/bin/bash
# -------------------------------------------------------------------
# eukfinder_shortreads_prepare.sh
#
# Description:
#   This script runs the short-read preparation step in Eukfinder.
#   It automatically activates the Eukfinder conda environment,
#   runs the read_prep command, and then deactivates the environment.
#
# Instructions:
#   1. Update the input file paths (R1, R2, host_genome, adapters, etc.)
#      to match your system.
#   2. Update pathes for the Centrifuge DB, host genome and adapter file.
#   3. Run: ./eukfinder_shortreads_prepare.sh
#
# Notes:
#   - For additional flags, see: eukfinder read_prep -h
#   - If running for the first time or updating taxonomy, consider using "-t True".
# -------------------------------------------------------------------

# 1. Define input data (Update these paths based on your setup)
R1="/path/to/your/R.1.fastq"
R2="/path/to/your/R.2.fastq"
host_genome="/path/to/your/host_genome.fasta"
centrifuge="/path/to/your/Centrifuge_DB/Centrifuge_DB"  # base name without .1.cf, etc.
adapters="/path/to/adapter_file"

# 2. Activate the Eukfinder Conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate eukfinder

# 3. Run Eukfinder (short-read preparation)
eukfinder read_prep \
    --r1 "$R1" \
    --r2 "$R2" \
    -n 20 \
    --hcrop 5 \
    -l 10 \
    -t 7 \
    --wsize 40 \
    --qscore 25 \
    --mlen 30 \
    --mhlen 30 \
    --hg "$host_genome" \
    -o "$output_dir" \
    --cdb "$centrifuge" \
    -i "$adapters"

# 5. Deactivate Conda environment
conda deactivate
