#!/bin/bash

## Launch_Unzip.sh
## Created by Adrianna Vandeuren on 02/12/2024
## This script serves to call Unzip_Job.sh.
## It will significantly reduce the runtime needed for this analysis. 


## Define the path to your sequencing results
projPath="/work/daylab/adrianna/RNA-seq_analysis/GSE101788"
pairedPath="/work/daylab/adrianna/RNA-seq_analysis/GSE101788/trimmed/Paired"

cd $pairedPath/

## Loop over all the files in the folder defined in your fastqPath
for file in *_R1_paired.fastq.gz
do
    ## Define the output prefix; -F is where the sample name will be cut
    OUTNAME=`ls ${file} | awk -F'_R1' '{print $1}'`

    ## Checkpoint to test code so far 
    echo "$OUTNAME"

    ## Define the read suffixes
    READ1=`ls ${OUTNAME}_R1_paired.fastq.gz`
    READ1=`ls ${OUTNAME}_R2_paired.fastq.gz`
    
    ## Checkpoint to test code so far
    ## echo "$READ1 $READ2"

    ## --export means that you pass on the defined arguments to the script you call here
    sbatch --export=R1=${READ1},R2=${READ2},NAME=${OUTNAME} /work/daylab/adrianna/RNA-seq_analysis/analysis/Unzip_Job.sh

done
