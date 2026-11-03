#!/bin/bash

## Launch_Align.sh
## Created by Adrianna Vandeuren on 01/25/2024
## This script serves to call Align_Job.sh.
## It will significantly reduce the runtime needed for this analysis. 


## Define the path to your sequencing results
projPath="work/daylab/Shannon/RNA-seq/Roach"
#trimmedPath="/work/daylab/Shannon/RNA-seq/Roach/samples/trimmed"
samplePath="/work/daylab/Shannon/RNA-seq/Roach/samples/bam"

#mkdir -p $trimmedPath/STAR_Aligned
mkdir -p $samplePath/STAR_Aligned
mkdir -p $projPath/counts

#cd $trimmedPath/
cd $samplePath/

## Loop over all the files in the folder defined in your fastqPath
for file in *.bam
do
    ## Define the output prefix; -F is where the sample name will be cut
    OUTNAME=`ls ${file} | awk -F'.bam' '{print $1}'`

    ## Checkpoint to test code so far 
    echo "$OUTNAME"

    ## --export means that you pass on the defined arguments to the script you call here
    sbatch --export=NAME=${OUTNAME} /work/daylab/Shannon/RNA-seq/scripts/Align_Job_SE.sh

done
