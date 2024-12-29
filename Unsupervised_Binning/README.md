# Supervised Binning for Recovering Eukaryotic Genomes from Metagenomes

## Overview

This workflow is designed to guide you through recovering eukaryotic genomes from metagenomic datasets using supervised binning, based on the contigs classified as eukaryotic or unknown (EUnk) by Eukfinder.

### Summary of Approach:
Input Selection: Contigs classified as eukaryotic or unknown (EUnk) are extracted from Eukfinder results for further analysis.
Multiple MyCC Analyses: Three separate MyCC analyses are conducted using different k-mers (4-mer, 5-mer, and 5-6-mer).
Data Mapping: Read coverage depth, identified rRNA sequences, and taxonomic identities of contigs are mapped to the corresponding MyCC bins.
Filtering and Classification: Contigs are filtered based on multiple criteria and assigned to final eukaryotic bins or marked as mitochondrial genomes.


1. run Plast, then using acc2tax to translate results into taxonomy
2. run Maxbin2, metabat2, MyCC (4mer, 5mer, 5/6mer)
3. run Metaxa2, Centrifuge. Parse all the results
4. Extract Mitochondiral genomes
5. Extract Eukaryotic genomes
