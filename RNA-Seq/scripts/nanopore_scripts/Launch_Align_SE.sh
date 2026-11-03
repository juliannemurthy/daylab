#!/bin/bash

# Launch_Align_SE.sh
# Created by Julianne (Yulia) Murthy on 06/24/2025

# Safety error catch
set -euo pipefail

# Create paths
projPath="/projects/daylab/yulia/RNA-Seq"
trimmedPath="$projPath/trimmed/RoachEtAl_trimmed"

# Go to the trimmed fastq path 
cd $trimmedPath

for file in *.fastq.gz; do
  sbatch $projPath/scripts/modified_scripts/Align_Job_SE.sh "$trimmedPath/$file"
done
