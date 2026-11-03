#!/bin/bash

## Launch_Align.sh
## Created by Adrianna Vandeuren on 01/25/2024
## This script serves to call Align_Job.sh.
## It will significantly reduce the runtime needed for this analysis. 


## Define the path to your sequencing results
projPath="/work/daylab/Shannon/RNA-seq/30-1172669254"
pairedPath="/work/daylab/Shannon/RNA-seq/30-1172669254/samples/trimmed/paired"

mkdir -p $projPath/STAR_Aligned
mkdir -p $projPath/analysis/counts

cd $pairedPath/

## Loop over all the files in the folder defined in your fastqPath
for file in *_R1_paired.fastq
do
    ## Define the output prefix; -F is where the sample name will be cut
    OUTNAME=`ls ${file} | awk -F'_R1_paired.fastq' '{print $1}'`

    ## Checkpoint to test code so far 
    echo "$OUTNAME"

    ## Define the read suffixes
    READ1=`ls ${OUTNAME}_R1_paired.fastq`
    READ2=`ls ${OUTNAME}_R2_paired.fastq`
    
    ## Checkpoint to test code so far
    echo "$READ1 $READ2"

    ## --export means that you pass on the defined arguments to the script you call here
    sbatch --export=R1=${READ1},R2=${READ2},NAME=${OUTNAME} /work/daylab/Shannon/RNA-seq/scripts/Align_Job.sh

done
