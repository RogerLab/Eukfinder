#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 20

cd $PWD
source activate Eukfinder
# Classification
plastdb=Plast_DB.fasta
plastmap=Plast_DB.map
centrifuge=Centrifuge_DB
acc2tax=Acc2Tax_DB

# previously prepared reads
R1=read_prep_p.1.fastq
R2=read_prep_p.2.fastq
uR1R2=read_prep_un.fastq

# centrifuge read classification
p_class=read_prep_centrifuge_P
up_class=read_prep_centrifuge_UP

eukfinder.py short_seqs --r1 $R1 --r2 $R2 --un $uR1R2 -n 20 -z 6 -t False --max_mem 100 \
                        -p $plastdb -m $plastmap -a $acc2tax -e 0.01 --pid 60 --cov 30 \
                        --pclass $p_class --uclass $up_class -o first_round --cdb $centrifuge \
                        --mhlen 30 -k 21,33,55
conda deactivate
