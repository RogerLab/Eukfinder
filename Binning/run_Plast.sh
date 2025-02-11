#!/bin/bash

# Activate the conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate eukfinder

# Input query and database
query=Eukfinder_long.fasta
DB=/scratch5/db/Eukfinder/nt2021/nt.fasta

# Run PLAST
plast -e 1E-5 \
      -max-hit-per-query 1 \
      -outfmt 1 \
      -a 48 \
      -p plastn \
      -i $query \
      -d $DB \
      -force-query-order 1000 \
      -o ${query::-6}.PLAST_nt.tsv

# Deactivate the conda environment
conda deactivate
