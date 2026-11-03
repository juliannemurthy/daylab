#!/bin/bash

# Launch_Trim_SE.sh
# Created by Julianne (Yulia) Murthy on 06/24/2025

# Safety error catch
set -euo pipefail

# Provide input and output paths
projPath="/projects/daylab/yulia/RNA-Seq"
fastqPath="$projPath/raw_data/RoachEtAl_fastq"
trimmedPath="$projPath/trimmed/RoachEtAl_trimmed"

# error handling: make output trimmed fastqc path directory if it does not exist already
mkdir -p $trimmedPath

# Go to the fastq data path
cd $fastqPath

# Loop that loops through every fastq.gz file in the path and strips its base name,
# designates outFile with new name using the base name,
# sends the trimming script to sbatch using each input and output file, running in parallel
for file in *.fastq.gz; do
  base=$(basename $file .fastq.gz)
  outFile="$trimmedPath/${base}_trimmed.fastq.gz"
  sbatch $projPath/scripts/modified_scripts/Trim_Job_SE.sh $fastqPath/$file $outFile
done
