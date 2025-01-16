#!/bin/bash

# Activate Conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate Eukfinder

# Input file and prefix for output
input=20.fastq
prefix=mytest_lr

# Path to databases (update these paths as needed)
plastdb=/path/to/Plast_DB.fasta         # Path to the PLAST database
plastmap=/path/to/Plast_DB.map          # Path to the PLAST map file
centrifuge=/path/to/Centrifuge_DB       # Path to the Centrifuge database (base name without .1.cf, .2.cf, etc.)
acc2tax=/path/to/Acc2Tax_DB             # Path to the Accession-to-Taxonomy database

# Information about usage
# Use `eukfinder.py -h` for details about the available submenus: read_prep, short_seqs, long_seqs
# Use `eukfinder.py long_seqs -h` to see flags for the `long_seqs` submenu
#
# WARNING: Set `-t True` only when the Eukfinder is used for the first time or when updating the taxonomy.
# Updating the taxonomy may take a long time (up to one hour) depending on your internet connection.

# Run Eukfinder for long sequences
eukfinder.py long_seqs \
    -l $input \
    -n 18 \
    -z 3 \
    -t False \
    -p $plastdb \
    -m $plastmap \
    -a $acc2tax \
    -e 0.01 \
    --pid 80 \
    --cov 30 \
    -o $prefix \
    --cdb $centrifuge \
    --mhlen 40

# Deactivate Conda environment
conda deactivate
