#!/bin/bash
mkdir "db"
cd "db/"
set -e  # Exit immediately if a command exits with a non-zero status.
set -u  # Treat unset variables as an error and exit immediately.

# URLs for the databases
ACC2TAX_URL="https://perun.biochem.dal.ca/Eukfinder/compressed_db/acc2tax_db.tar.gz/"
CENTRIFUGE_URL="https://perun.biochem.dal.ca/Eukfinder/compressed_db/centrifuge_db.tar.gz/"
PLAST_URL="https://perun.biochem.dal.ca/Eukfinder/compressed_db/plast_db.tar.gz/"
HUMAN_GENOME_URL="https://perun.biochem.dal.ca/Eukfinder/compressed_db/GCF_000001405.39_GRCh38.p13_human_genome.fna.tar.gz"
READ_ADAPTERS_URL="https://perun.biochem.dal.ca/Eukfinder/TrueSeq2_NexteraSE-PE.fa"

# Sizes of the databases
ACC2TAX_SIZE="36 GB"
CENTRIFUGE_SIZE="70 GB"
PLAST_SIZE="3.7 GB"
HUMAN_GENOME_SIZE="0.92 GB"
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
echo "Please enter the number of compressed database(s) which you would like to install, separated by spaces (e.g., 1 2)."
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
for file in $file_selection; do
  case $file in
    0)
      echo "Exiting..."
      cd ..
      rmdir "db/"
      exit 0
      ;;
    1)
      download_database $ACC2TAX_URL "acc2tax"
      ;;
    2)
      download_database $CENTRIFUGE_URL "centrifuge"
      ;;
    3)
      download_database $PLAST_URL "PLAST"
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
done

echo "All selected databases have been downloaded."
echo "Decompressing..."
tar -xzvf *.tar.gz
echo "Cleaning up..."
rm -f *.tar.gz
echo "Download completed."