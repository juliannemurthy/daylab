#!/bin/bash

# Launch_Trim_PE.sh
# Modified by Julianne (Yulia) Murthy on 07/09/2025
# Created by Adrianna Vandeuren 08/02/2023

# Safety error catch
set -euo pipefail

# Provide input and output paths
projPath="/projects/daylab/yulia/RNA-Seq"
fastqPath="$projPath/raw_data/L1vsRA_fastq"
pairedPath="$projPath/trimmed/L1vsRA_trimmed/paired"
unpairedPath="$projPath/trimmed/L1vsRA_trimmed/unpaired"
trimmomaticPath="/shared/centos7/anaconda3/2021.11/envs/BINF-12-2021/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2"

# error handling: make output trimmed fastqc path directory if it does not exist already
mkdir -p $pairedPath
mkdir -p $unpairedPath

# Go to the fastq data path
cd $fastqPath

# Loop over all the files in the folder defined in your fastqPath
for file in *_R1_001.fastq.gz
do
    # Define the output prefix; -F is where the sample name will be cut
    OUTNAME=$(basename "$file" _R1_001.fastq.gz)

    # Define the read suffixes
    READ1=${OUTNAME}_R1_001.fastq.gz
    READ2=${OUTNAME}_R2_001.fastq.gz

    # --export means that you pass on the defined arguments to the script you call here
    sbatch $projPath/scripts/illumina_scripts/Trim_Job_PE.sh \
     "$READ1" "$READ2" "$OUTNAME" "$trimmomaticPath" "$pairedPath" "$unpairedPath"
done
