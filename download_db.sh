#!/bin/bash
mkdir "db"
cd "db/"
set -e  # Exit immediately if a command exits with a non-zero status.
set -u  # Treat unset variables as an error and exit immediately.

# URLs for the databases
ACC2TAX_DIR_URL="https://perun.biochem.dal.ca/Eukfinder/Acc2tax/"
CENTRIFUGE_DIR_URL="https://perun.biochem.dal.ca/Eukfinder/Centrifuge_DB/"
PLAST_DIR_URL="https://perun.biochem.dal.ca/Eukfinder/PLAST_DB/"
HUMAN_GENOME_URL="https://perun.biochem.dal.ca/Eukfinder/GCF_000001405.39_GRCh38.p13_human_genome.fna"
READ_ADAPTERS_URL="https://perun.biochem.dal.ca/Eukfinder/TrueSeq2_NexteraSE-PE.fa"

# Sizes of the databases
ACC2TAX_SIZE="221.4 GB"
CENTRIFUGE_SIZE="91.4 GB"
PLAST_SIZE="10.2 GB"
HUMAN_GENOME_SIZE="3.1 GB"
READ_ADAPTERS_SIZE="16.0 KB"

# Display the available databases and their sizes
echo ""
echo "1. acc2tax database - $ACC2TAX_SIZE"
echo "2. centrifuge database - $CENTRIFUGE_SIZE"
echo "3. PLAST database - $PLAST_SIZE"
echo "4. Human Genome for read decontamination - $HUMAN_GENOME_SIZE"
echo "5. Read Adapters for Illumina sequencing - $READ_ADAPTERS_SIZE"

echo ""
# Prompt the user to enter the number of databases to install
echo "Please enter the number of database(s) which you would like to install, separated by spaces (e.g., 1 2)."
echo "Or 0, if you would like to skip for now:"
read -r file_selection
echo ""

# Function to download a database
download_database() {
  local db_url=$1
  local db_name=$2
  echo "Downloading $db_name..."
  mkdir $db_name
  cd $db_name
  wget -r -np -nH $db_url
  echo "$db_name downloaded successfully."
  cd ..
}

# Loop through the user input and download the selected databases
for file in $file_selection; {
  case $file in
    0)
      echo "Exiting..."
      cd ..
      rmdir "db/"
      exit 0
      ;;
    1)
      download_database $ACC2TAX_DIR_URL "acc2tax"
      ;;
    2)
      download_database $CENTRIFUGE_DIR_URL "centrifuge"
      ;;
    3)
      download_database $PLAST_DIR_URL "PLAST"
      ;;
    4)
      wget $HUMAN_GENOME_URL
      ;;
    5)
      wget $READ_ADAPTERS_URL
      ;;
    *)
      echo "Invalid selection: $file"
      exit -1
      ;;
  esac
}

echo "All selected databases have been downloaded."
