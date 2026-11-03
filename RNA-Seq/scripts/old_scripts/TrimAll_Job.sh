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
fastqPath="/work/daylab/Shannon/RNA-seq/30-1172669254/samples"
pairedPath="/work/daylab/Shannon/RNA-seq/30-1172669254/samples/trimmed/paired"
unpairedPath="/work/daylab/Shannon/RNA-seq/30-1172669254/samples/trimmed/unpaired"
PATH_TO_TRIMMOMATIC="/shared/centos7/anaconda3/2021.11/envs/BINF-12-2021/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2"

## Activate the environment that contains the tools you need
module load anaconda3
source activate RNA-seq

cd $fastqPath
## Since trimmomatic is installed in the RNA-seq environment we do not have to call the jar file
## The LEADING, TRAILING and MINLEN are recommended through the Azenta analysis guide, since we are sequencing with them it seems wise to follow their recommendations. 

#trimmomatic PE -phred33 ${R1} ${R2} ${pairedPath}/${NAME}_R1_paired.fastq.gz ${unpairedPath}/${NAME}_R1_unpaired.fastq.gz ${pairedPath}/${NAME}_R2_paired.fastq.gz ${unpairedPath}/${NAME}_R2_unpaired.fastq.gz ILLUMINACLIP:$PATH_TO_TRIMMOMATIC/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

java -jar $PATH_TO_TRIMMOMATIC/trimmomatic.jar PE ${R1} ${R2} ${pairedPath}/${NAME}_R1_paired.fastq.gz ${unpairedPath}/${NAME}_R1_unpaired.fastq.gz ${pairedPath}/${NAME}_R2_paired.fastq.gz ${unpairedPath}/${NAME}_R2_unpaired.fastq.gz ILLUMINACLIP:$PATH_TO_TRIMMOMATIC/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

