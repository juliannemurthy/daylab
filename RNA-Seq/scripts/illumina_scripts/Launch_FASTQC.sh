#!/bin/bash

## Launch_FASTQC.sh
## Modifed by Julianne (Yulia) Murthy on 07/09/2025
## Created by Adrianna Vandeuren on 08/02/2023
## This script serves to launch the FASTQC jobs in parallel on all files.
## It will significantly reduce the runtime needed for this analysis. 

# Safety error catch
set -euo pipefail

# Provide input and output paths
projPath="/projects/daylab/yulia/RNA-Seq"
fastqPath="$projPath/raw_data/L1vsRA_fastq"
qcPath="$projPath/QC/L1vsRA_QC"

# error handling: make output fastqc path directory if it does not exist already
mkdir -p $qcPath

# Go to the path 
cd $fastqPath

# Loop that sends each file for FASTQC analysis using sbatch
# View results by looking at sbatch err and log files 
for file in *.fastq.gz; do
  sbatch $projPath/scripts/illumina_scripts/FASTQC_Job.sh $fastqPath $qcPath $file
done
