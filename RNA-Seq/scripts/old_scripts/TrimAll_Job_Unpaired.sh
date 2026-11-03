#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=TrimAll
#SBATCH --time=4:00:00
#SBATCH --mem=10G
#SBATCH --partition=short

## TrimAll_Job.sh 
## Created by Adrianna Vandeuren 01/24/2024
## This script will serve to remove the adapter sequences from your reads as well as to remove the poor-quality regions. 

source /etc/profile.d/modules.sh

## Define necessary paths 
projPath="/work/daylab/Shannon/RNA-seq/Roach/samples"
trimmedPath="/work/daylab/Shannon/RNA-seq/Roach/samples/trimmed"

## Activate the environment that containes the tools you need
module load anaconda3
source activate RNA-seq

## Since trimmomatic is installed in the RNA-seq environment we do not have to call the jar file
## The LEADING, TRAILING and MINLEN are recommended through the Azenta analysis guide, since we are sequencing with them it seems wise to follow their recommendations. 

trimmomatic SE -phred33 ${NAME}.fastq ${trimmedPath}/${NAME}_trimmed.fastq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36 
#trimmomatic-0.35.jar SE -phred33 input.fq.gz output.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#trimmomatic SE ${NAME}.fastq.gz ${trimPath}/${NAME}_trim.fastq.gz -phred33 ILLUMINACLIP:TruSeq3-SE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36 



