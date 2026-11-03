#!/bin/bash

# Launch_Align_PE.sh
# Created by Julianne (Yulia) Murthy on 07/09/2025

# Safety error catch
set -euo pipefail

# Create paths
projPath="/projects/daylab/yulia/RNA-Seq"
pairedPath="$projPath/trimmed/L1vsRA_trimmed/paired"

# Go to the trimmed fastq path 
cd $pairedPath

# Loop that strips base name from each paired file, identifies read 1 and read 2 files, submits batch job so reads run in parallel
for R1 in *_R1_paired.fastq.gz; do
  base=$(basename "$R1" _R1_paired.fastq.gz)
  R2="${base}_R2_paired.fastq.gz"
  sbatch $projPath/scripts/illumina_scripts/Align_Job_PE.sh "$pairedPath/$R1" "$pairedPath/$R2" "$base"
done
