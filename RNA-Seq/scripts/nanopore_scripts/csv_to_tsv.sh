#!/bin/bash
#SBATCH --job-name=csv_to_tsv.sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00

## csv_to_tsv.sh
## Created by Ella Chee on 10/23/24.
## This script converts CSV files to TSV (comma-separated to tab-separated)
## Usage: ./csv_to_tsv.sh /path/to/csv_directory

## Check if the directory argument is provided
if [ -z "$1" ]; then
  echo "Usage: $0 /path/to/csv_directory"
  exit 1
fi

## Define the path to raw files
bedPath=$1

## Check if the path exists
if [ ! -d "$bedPath" ]; then
  echo "Error: Directory $bedPath does not exist."
  exit 1
fi

## Create output directory if it doesn't exist
mkdir -p $bedPath/tsv

cd $bedPath

## Loop over all the CSV files in the folder
for file in *.csv; do
    ## Extract filename without extension
    OUTNAME=$(basename "$file" .csv)
    
    ## Convert CSV to TSV and write to the output directory
    awk -F ',' 'BEGIN {OFS="\t"} {$1=$1}1' "$file" > "$bedPath/tsv/${OUTNAME}.tsv"
done

echo "All files converted to TSVs :)"
