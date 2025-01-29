#!/bin/bash
# -------------------------------------------------------------------
# eukfinder_longreads_classify.sh
# 
# Description:
#   This script classifies long reads using Eukfinder. It automatically
#   activates the Eukfinder conda environment, runs the classification,
#   and then deactivates the environment.
#
# Instructions:
#   1. Update the input file paths (input, prefix) to match your data.
#   2. Update the database paths as needed.
#   3. Run: ./eukfinder_short_classify.sh
#
# Notes:
#   - Make sure your environment is correctly named ("eukfinder" here).
#   - Update the database paths before running.
#   - For first-time or taxonomy updates, set "-t True" in eukfinder.
# -------------------------------------------------------------------

# 1. Define input files and output prefix
input="read_fastq_or_assembly_fasta"
prefix="Test_long"

# 2. Full paths to databases (edit these as needed)
plastdb="/path/to/Plast_DB.fasta"       # Path to the PLAST database
plastmap="/path/to/Plast_DB.map"        # Path to the PLAST map file
centrifuge="/path/to/Centrifuge_DB"     # base name (no .1.cf, .2.cf, etc.)
acc2tax="/path/to/Acc2Tax_DB"           # Path to folder of acc2tax database

# 3. Activate Conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate Eukfinder

# 4. Run Eukfinder for long reads
#    - See `eukfinder long_seqs -h` for additional flags.
eukfinder long_seqs \
    -l "$input" \
    -n 18 \
    -z 3 \
    -t False \
    -p "$plastdb" \
    -m "$plastmap" \
    -a "$acc2tax" \
    -e 0.01 \
    --pid 80 \
    --cov 30 \
    -o "$prefix" \
    --cdb "$centrifuge" \
    --mhlen 40

# 5. Deactivate Conda environment
conda deactivate
