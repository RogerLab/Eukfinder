#!/bin/bash

# Activate Conda environment for BLAST
source $(conda info --base)/etc/profile.d/conda.sh
conda activate blast

# Define input query file
query=Eukfinder_long.fasta

# Define the BLAST database path and name
# Users must modify the `DB` variable to point to their local BLAST database.
DB=/path/to/your/blast_database/mito_blast_db

# Check if the database exists
if [ ! -f "${DB}.nin" ]; then
  echo "Error: BLAST database '${DB}' not found. Please specify the correct path."
  exit 1
fi

# Run BLAST
blastn -db $DB \
       -query $query \
       -out ${query::-6}_BLAST4Mit.out \
       -num_threads 30 \
       -outfmt "6 qseqid sseqid stitle evalue pident qcovhsp nident mismatch length slen qlen qstart qend sstart send staxids sscinames sskingdoms" \
       -evalue 1E-5 \
       -max_hsps 1

# Deactivate Conda environment
conda deactivate
