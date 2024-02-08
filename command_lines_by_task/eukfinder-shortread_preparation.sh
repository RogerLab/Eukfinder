#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 42
cd $PWD

# input data
R1=EUnkBact_overlap.al-con.ht2.R.1.fastq
R2=EUnkBact_overlap.al-con.ht2.R.2.fastq
host_genome=test.host.fasta
centrifuge=/scratch2/software/centrifuge-1.0.3/Centrifuge_Apr7_DB
adapters=/home/dsalas/anvioWrap/TrueSeq2_NexteraSE-PE.fa

# command line for read preparation
source activate Metagenomics-Scavenger
eukfinder.py read_prep --r1 $R1 --r2 $R2 -n 42 --hcrop 5 -l 10 -t 7 \
                       --wsize 40 --qscore 25 --mlen 30 --mhlen 30 \
                       --hg $host_genome -o read_prep \
                       --cdb $centrifuge -i $adapters
conda deactivate
