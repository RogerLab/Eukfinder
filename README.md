# eukfinder
[![Bioconda platforms](https://img.shields.io/conda/pn/bioconda/eukfinder?style=flag)](https://anaconda.org/bioconda/eukfinder)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/eukfinder.svg?style=flag&label=Bioconda%20install)](https://anaconda.org/bioconda/eukfinder)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/RogerLab/Eukfinder/blob/main/LICENSE.txt)
[![CI](https://github.com/RogerLab/Eukfinder/workflows/Build/badge.svg)](https://github.com/RogerLab/Eukfinder/actions)

A bioinformatics tool designed to enable rapid recovery of high-quality draft genomes of microbial eukaryotes from various environmental metagenomic samples.

- reference-independent and cultivation-free
- separate reads or contigs into five groups: Bacteria, Archaeal, Eukaryotic, Viral, and Unknown. 
- generate a fasta file with Euk and unknown contigs for binning


## Overview
Schematic representation of eukfinder workflows. eukfinder is a taxonomic classification-based bioinformatics approach to retrieve microbial eukaryotic nuclear and mitochondrial genomes from WGS metagenomic sequencing data. eukfinder has two different workflows based on the input files:

- (a) using Illumina short reads, it first classifies Illumina reads into 5 distinct taxonomic categories (Archaeal, Bacterial, Viral, Eukaryotic, and Unknown), after assembling the Eukaryotic and Unknown reads together, second round of classification and supervised binning, it will output Eukaryotic nuclear and mitochondrial genomes.
- (b) using MAG assembled contigs or long-read sequencing data generated by Nanopore or Pacbio platforms. It goes through one round of classification to select Eukaryotic and Unknown contigs, and after supervised binning generates Eukaryotic nuclear and mitochondrial genomes.

![Graphical_abstract](https://github.com/RogerLab/Eukfinder/assets/39600837/1d9e690e-be40-4255-b00a-07742219d92e)

## Installation
Anaconda or miniconda required*
### 1. Created environment and install eukfinder

```sh
conda create -n eukfinder -c bioconda eukfinder
```

### 2. Download or build databases
 
  Default reference databases can be downloaded from [Eukfinder Databases](https://perun.biochem.dal.ca/Eukfinder/)
- Plast Database
- Centrifuge Database
- acc2tax Database 
- Human Genome for read decontamination
- Read Adapters for Illumina sequencing

```shell
./download_db.sh
```

 Users can flexibly customize the reference data (see [here](https://github.com/dzhao2019/eukfindertest/wiki/Build-a-customized-reference-database))
 
### 3. (!) Activate eukfinder environment before running the command

```shell
source activate eukfinder
```

<!-- USAGE EXAMPLES -->
## Usage
See [Wiki](https://github.com/dzhao2019/eukfindertest/wiki) for detailed description

**eukfinder read_prep**

Run Trimmomatic to remove low-quality reads, and adaptors
Run Bowtie2 to remove host reads
Run Centrifuge for the first round of classification

    eukfinder read_prep [-h] --r1 R1 --r2 R2 -n THREADS -i ILLUMINA_CLIP
                               --hcrop HCROP -l LEADING_TRIM -t TRAIL_TRIM --wsize
                               WSIZE --qscore QSCORE --mlen MLEN --hg HG -o
                               OUT_NAME --cdb CDB

  
**eukfinder short_seqs**

    eukfinder short_seqs [-h] --r1 R1 --r2 R2 --un UN -o OUT_NAME -n
                            NUMBER_OF_THREADS -z NUMBER_OF_CHUNKS -t
                            TAXONOMY_UPDATE -p PLAST_DATABASE -m PLAST_ID_MAP
                            [-p2 ANCILLARY_PLAST_DATABASE]
                            [-m2 ANCILLARY_PLAST_ID_MAP]
                            [--force-pdb FORCE_PDB] -a ACC2TAX_DATABASE --cdb
                            CDB -e E_VALUE --pid PID --cov COV --max_m MAX_M
                            --mhlen MHLEN --pclass PCLASS --uclass UCLASS


**eukfinder long_seqs**

    eukfinder long_seqs [-h] -l LONG_SEQS -o OUT_NAME --mhlen MHLEN --cdb
                           CDB -n NUMBER_OF_THREADS -z NUMBER_OF_CHUNKS -t
                           TAXONOMY_UPDATE -p PLAST_DATABASE -m PLAST_ID_MAP
                           -a ACC2TAX_DATABASE -e E_VALUE --pid PID --cov COV

<!-- CONTRIBUTING -->
## Contributing

Contributions are what makes the open-source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

<!-- Publication -->
## Publication/Citation

Zhao, D., Salas-Leiva, D.E., Williams, S.K., Dunn, K.A. and Roger, A.J., 2023. Eukfinder: a pipeline to retrieve microbial eukaryote genomes from metagenomic sequencing data. bioRxiv, pp.2023-12.


<!-- CONTACT -->
## Contact

d.zhao@dal.ca

ds2000@cam.ac.uk

