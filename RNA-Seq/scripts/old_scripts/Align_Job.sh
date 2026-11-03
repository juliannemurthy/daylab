#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=Align
#SBATCH --time=20:00:00
#SBATCH --mem=50G
#SBATCH --partition=short
#SBATCH --cpus-per-task=20

## Align_Job.sh 
## Created by Adrianna Vandeuren 01/25/2024
## This script will serve to align your paired reads to the indexed reference genome

## Define necessary paths 
projPath="/work/daylab/Shannon/RNA-seq/30-1172669254/samples"
pairedPath="/work/daylab/Shannon/RNA-seq/30-1172669254/samples/trimmed/paired"
unpairedPath="/work/daylab/Shannon/RNA-seq/30-1172669254/samples/trimmed/unpaired"
alignedPath="/work/daylab/Shannon/RNA-seq/30-1172669254/STAR_Aligned"
WBcel235Path="/work/daylab/Shannon/RNA-seq/geneAnnotations/WBcel235"

## Activate the environment that containes the tools you need
source /etc/profile.d/modules.sh
module load anaconda3
conda init
source activate RNA-seq

## Align the paired reads to the hg19 genome 
STAR --runThreadN 20 --genomeDir $WBcel235Path --readFilesIn ${R1} ${R2} --outFileNamePrefix $alignedPath/${NAME} --outSAMtype BAM Unsorted 

## Sort the bam files
samtools sort $alignedPath/${NAME}Aligned.out.bam -o $alignedPath/${NAME}.sorted.bam

## Index the sorted bam files
samtools index $alignedPath/${NAME}.sorted.bam

## Calculate the gene hit counts with FeatureCounts
##featureCounts -a ${hg19Path}/gencode.v19.chr_patch_hapl_scaff.annotation.gtf -o $projPath/Counts/${NAME}.counts.txt --largestOverlap -t exon -g gene_id $alPath/${NAME}.sorted.bam

featureCounts -a /work/daylab/Shannon/RNA-seq/geneAnnotations/WBcel235genomic.gtf -o /work/daylab/Shannon/RNA-seq/30-1172669254/analysis/counts/${NAME}.counts.txt --largestOverlap -t exon -g gene_id /work/daylab/Shannon/RNA-seq/30-1172669254/STAR_Aligned/${NAME}.sorted.bam
