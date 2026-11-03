#!/bin/bash

## Launch_NANO_QC.sh
## Modifed by Julianne (Yulia) Murthy on 06/24/2025
## Created by Adrianna Vandeuren on 08/02/2023
## This script serves to launch the Oxford Nanopore FASTQ QC jobs in parallel on all files.
## It will significantly reduce the runtime needed for this analysis. 

# Provide input and output paths
projPath="/projects/daylab/yulia/RNA-Seq"
fastqPath="$projPath/raw_data/RoachEtAl_fastq"
qcPath="$projPath/QC"

# error handling: make output fastqc path directory if it does not exist already
mkdir -p $qcPath

# Go to the path 
cd $fastqPath

# Loop that sends each file for QC analysis using sbatch
# View results by looking at sbatch err and log files 
for file in *.fastq.gz; do
  sbatch $projPath/scripts/modified_scripts/NANO_QC.sh $fastqPath $qcPath $file
done
