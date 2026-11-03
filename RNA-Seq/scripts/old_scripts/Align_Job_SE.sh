#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=Align
#SBATCH --time=20:00:00
#SBATCH --mem=50G
#SBATCH --partition=short

## Align_Job.sh 
## Created by Adrianna Vandeuren 01/25/2024
## This script will serve to align your paired reads to the indexed reference genome

## Define necessary paths 
projPath="/work/daylab/Shannon/RNA-seq/Roach"
trimmedPath="/work/daylab/Shannon/RNA-seq/Roach/samples"
alPath="/work/daylab/Shannon/RNA-seq/Roach/samples/bam/STAR_Aligned"
bamPath="/work/daylab/Shannon/RNA-seq/Roach/samples/bam"
#ce11Path="/work/daylab/Shannon/CT/referenceGenomes/ce11"
cel235Path="/work/daylab/Shannon/RNA-seq/geneAnnotations"
cel235PathGA="/work/daylab/Shannon/RNA-seq/geneAnnotations/WBcel235"

## Activate the environment that contains the tools you need
source /etc/profile.d/modules.sh
module load anaconda3
conda init
source activate RNA-seq

## Align the paired reads to the WBcel235 genome 
#STAR --genomeDir $cel235Path/WBcel235 --readFilesIn ${NAME}_trimmed.fastq --outFileNamePrefix $alPath/${NAME} --outSAMtype BAM Unsorted 
#STAR --genomeDir $cel235Path/WBcel235 --readFilesIn ${NAME}.fastq --outFileNamePrefix $alPath/${NAME} --outSAMtype BAM Unsorted 
## Sort the bam files
#samtools sort $alPath/${NAME}Aligned.out.bam -o $alPath/${NAME}.sorted.bam
samtools sort $bamPath/${NAME}.bam -o $bamPath/sorted/${NAME}.sorted.bam
## Index the sorted bam files
samtools index $bamPath/sorted/${NAME}.sorted.bam

## Calculate the gene hit counts with FeatureCounts
featureCounts -a ${cel235PathGA}/WBcel235genomic.gtf -o $projPath/counts/${NAME}.counts.txt --largestOverlap -t exon -g gene_id $bamPath/sorted/${NAME}.sorted.bam
