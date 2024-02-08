#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 42

cd $PWD
source activate Metagenomics-Scavenger
# Classification
plastdb=/scratch4/dsalas/Shared/EukfinderLabDB/GenericDB.lab.fasta
plastmap=/scratch4/dsalas/Shared/EukfinderLabDB/GenericDB.lab.map
centrifuge=/scratch2/software/centrifuge-1.0.3/Centrifuge_Apr7_DB
acc2tax=/db1/extra-data-sets/Acc2tax/Acc2Tax_04_01_2024/

# previously prepared reads
R1=/misc/scratch4/dsalas/Eukfinder_2024/test1/read_prep_p.1.fastq
R2=/misc/scratch4/dsalas/Eukfinder_2024/test1/read_prep_p.2.fastq
uR1R2=/misc/scratch4/dsalas/Eukfinder_2024/test1/read_prep_un.fastq

# centrifuge read classification
p_class=/misc/scratch4/dsalas/Eukfinder_2024/test1/read_prep_centrifuge_P
up_class=/misc/scratch4/dsalas/Eukfinder_2024/test1/read_prep_centrifuge_UP

eukfinder.py short_seqs --r1 $R1 --r2 $R2 --un $uR1R2 -n 42 -z 6 -t False --max_mem 100 \
                        -p $plastdb -m $plastmap -a $acc2tax -e 0.01 --pid 60 --cov 30 \
                        --pclass $p_class --uclass $up_class -o first_round --cdb $centrifuge \
                        --mhlen 30 -k 21,33,55
conda deactivate
