#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 20

cd $PWD
source activate Metagenomics-Scavenger
input=20.fastq
prefix=mytest_lr
#path to the eukfinder script
# path to databases
acc2tax=/db1/extra-data-sets/Acc2tax/Acc2Tax_04_01_2024/
# this can be customized similarly as in ReadClassifier script
plastdb=/scratch4/dsalas/Shared/EukfinderLabDB/GenericDB.lab.fasta
plastmap=/scratch4/dsalas/Shared/EukfinderLabDB/GenericDB.lab.map
centrifuge=/scratch2/software/centrifuge-1.0.3/Centrifuge_Apr7_DB

# information about flags:
#     $home/eukfinder_pre-class.py -h  --> will indicate the available submenus: read_prep,short_seqs,long_seqs)
#     $home/eukfinder_pre-class.py long_seqs -h   --> will provide information about the flags for that particular submenu
#  *** WARNING: set -t True only when the eukfinder is used for first time or if you really need to update the taxonomy. If the latter 
#               is the case, the update will take a long time (up to one hour depending on a steady internet connection. #
#

eukfinder.py long_seqs -l $input -n 18 -z 3 -t False -p $plastdb -m $plastmap -a $acc2tax \
                       -e 0.01 --pid 80 --cov 30 -o $prefix --cdb $centrifuge --mhlen 40
conda deactivate
