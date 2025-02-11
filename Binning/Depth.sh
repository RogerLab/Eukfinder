#!/bin/bash

# Activate Conda environment for Bowtie2
source $(conda info --base)/etc/profile.d/conda.sh
conda activate eukfinder

# Read files (Paired-end reads and unpaired reads) and assembly file
R1=read_prep_p.1.fastq
R2=read_prep_p.2.fastq
uR1R2=read_prep_un.fastq
assembly=Eukfinder_long.fasta
index=Eukfinder_long_EUnk

# Build the Bowtie2 index
bowtie2-build -f $assembly $index

# Mapping reads using Bowtie2
basename=Eukfinder_long_EUnk
bowtie2 --threads 20 -x $index -1 $R1 -2 $R2 -U $uR1R2 -S $basename.sam --no-unal

# Activate Conda environment for SAMtools
conda activate samtools

# Convert SAM to BAM
samtools view --threads 20 -bS $basename.sam > $basename.bam

# Sort the BAM file
samtools sort --threads 20 -o ${basename}_sorted.bam $basename.bam

# Index the sorted BAM file
samtools index ${basename}_sorted.bam

# Activate Conda environment for MetaBAT2
conda activate metabat2

# Calculate depth using jgi_summarize_bam_contig_depths
jgi_summarize_bam_contig_depths --outputDepth ${basename}.depth.txt ${basename}_sorted.bam

# Deactivate Conda
conda deactivate
