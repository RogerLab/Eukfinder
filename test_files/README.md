 Test files can be downloaded from [test files](https://perun.biochem.dal.ca/Metagenomics-Scavenger/test_files/)

# Examples

## 1. Preparation steps for short reads: run Trimmomatics, Bowtie2, Centrifuge
###   Input files: paired fastq files and host genome [short_reads](https://perun.biochem.dal.ca/Metagenomics-Scavenger/test_files/short_reads_workflow/)
   EUnkBact_overlap.al-con.ht2.R.1.fastq <br>
   EUnkBact_overlap.al-con.ht2.R.2.fastq	<br>
   test.host.fasta 
   adaptor.fa
   <br>
###   Output files:
Two paired fastq files and one unpaired fastq file
Two centrifuge result files for paired and unpaired reads

## 2. Short reads classification
###   Input files: cleaned reads 
(2 paired fastq files and 1 unpaired fastq file, centrifuge results for paired and unpaired reads)  [short_reads](https://perun.biochem.dal.ca/Metagenomics-Scavenger/test_files/short_reads_workflow/)
<br>
###   Output files:
Directory Name: TempEukfinder/Classified_contigs
<br>
Five fasta files containing bacterial, archeal, eukaryotic, viral, and eukaryotic+unknown contigs:
<br>
Arch.fasta, Bact.fasta, Euk.fasta, Misc.fasta, EUnk.fasta 
<br>

## 3. Long reads classification
###   Input files: long reads fasta file [long reads](https://perun.biochem.dal.ca/Metagenomics-Scavenger/test_files/long_reads_workflow/)
   longreads.fastq <br>
   
###   Output files:

Directory Name: TempEukfinder
<br>
Five fasta files containing bacterial, archeal, eukaryotic, viral, and eukaryotic+unknown contigs:
<br>
Arch.fasta, Bact.fasta, Euk.fasta, Misc.fasta, EUnk.fasta 
<br>

## 4. Assembled genome classification
   Input files: assembly (fasta file) <br>



###   Output files:

Directory Name: TempEukfinder
<br>
Five fasta files containing bacterial, archeal, eukaryotic, viral, and eukaryotic+unknown contigs:
<br>
Arch.fasta, Bact.fasta, Euk.fasta, Misc.fasta, EUnk.fasta    
<br>
