#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 10
cd $PWD
input=Read2M_Test1_scaffolds_Eukfinder_long_EUnk.fasta

source activate maxbin2
# Run maxbin2
echo "======== Run maxbin2 ========"
input=Read2M_Test1_scaffolds_Eukfinder_long_EUnk.fasta
run_MaxBin.pl -contig $input -abund depth.txt -out maxbin2_binning -thread 10

source activate metabat2
# Run metabat2
echo "======== Run metabat2 ========"
metabat2 -i $input -o metabat2_binning -m 1500 -t 10  --unbinned
