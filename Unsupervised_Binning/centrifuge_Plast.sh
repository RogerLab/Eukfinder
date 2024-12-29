#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 20
cd $PWD
source activate eukfinder
# Run centrifuge
centDB=/scratch3/Eukfinder/DB/Centrifuge_DB/Centrifuge_NewDB_Sept2020
centrifuge -q -- phred33 -x $centDB -U Eukfinder_long_EUnk.fasta -p 20 -S Eukfinder_long_EUnk_centrifuge --report-file Eukfinder_long_EUnk_centrifuge_report.tsv
# Run plast
DB=/scratch5/db/Eukfinder/nt2022/nt
plast -e 1E-5 -max-hit-per-query 1 -outfmt 1 -a 20 -p plastn -i Eukfinder_long_EUnk.fasta -d $DB -force-query-order 1000 -o Eukfinder_long_EUnk_PLAST.tsv

