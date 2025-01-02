#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 20
cd $PWD
source activate eukfinder
# Read files (Paired-end reads and unpaired reads) and assembly file
R1=read_prep_p.1.fastq
R2=read_prep_p.2.fastq
uR1R2=read_prep_un.fastq
assembly=Eukfinder_long_EUnk.fasta
index=Eukfinder_long_EUnk

# Build the index 
bowtie2-build -f $assembly $index
# Mapping reads using Bowtie2
basename=Eukfinder_long_EUnk
bowtie2 --threads 20 -x $index -1 $R1 -2 $R2 -U $uR1R2 -S $basename.sam  --no-unal
# Convert SAM to BAM
source activate samtools
samtools view --threads 20 -bS $basename.sam > $basename.bam

# Sort the BAM file
samtools sort --threads 20 -o $basename_sorted.bam $basename.bam

# Index the sorted BAM file
samtools index $basename_sorted.bam


# Calculate depth using jgi_summarize_bam_contig_depths
source activate metabat2
jgi_summarize_bam_contig_depths --outputDepth $basename.depth.txt $basename_sorted.bam

conda deactivate
