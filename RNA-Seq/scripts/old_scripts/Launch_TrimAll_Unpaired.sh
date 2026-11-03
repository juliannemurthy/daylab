#!/bin/bash

## Launch_TrimAll.sh
## Created by Adrianna Vandeuren on 01/24/2024
## This script serves to call TrimAll_Job.sh.
## It will significantly reduce the runtime needed for this analysis. 


## Define the path to your sequencing results
projPath="/work/daylab/Shannon/RNA-seq/Roach"
fastqPath="/work/daylab/Shannon/RNA-seq/Roach/samples/"

## Create the output directories
mkdir -p $fastqPath/trimmed
cd $fastqPath

## Loop over all the files in the folder defined in your fastqPath
for file in *.fastq
do
    ## Define the output prefix; -F is where the sample name will be cut
    OUTNAME=`ls ${file} | awk -F'.fastq' '{print $1}'`

    ## Checkpoint to test code so far 
    echo "$OUTNAME"

    ## --export means that you pass on the defined arguments to the script you call here
    sbatch --export=NAME=${OUTNAME} /work/daylab/Shannon/RNA-seq/scripts/TrimAll_Job_Unpaired.sh

done
