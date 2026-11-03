#!/bin/bash

## Launch_Unzip.sh
## Created by Adrianna Vandeuren on 02/12/2024
## This script serves to call Unzip_Job.sh.
## It will significantly reduce the runtime needed for this analysis. 


## Define the path to your sequencing results
projPath="/work/daylab/adrianna/RNA-seq_analysis/GSE85443"
pairedPath="/work/daylab/adrianna/RNA-seq_analysis/GSE85443/trimmed"

cd $pairedPath/

## Loop over all the files in the folder defined in your fastqPath
for file in *_trim.fastq.gz
do
    ## Define the output prefix; -F is where the sample name will be cut
    OUTNAME=`ls ${file} | awk -F'_trim' '{print $1}'`

    ## Checkpoint to test code so far 
    echo "$OUTNAME"


    ## --export means that you pass on the defined arguments to the script you call here
    sbatch --export=NAME=${OUTNAME} /work/daylab/adrianna/RNA-seq_analysis/analysis/Unzip_Job_SE.sh

done
