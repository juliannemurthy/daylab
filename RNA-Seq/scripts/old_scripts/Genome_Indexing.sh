#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=Genome_Index
#SBATCH --time=15:00:00
#SBATCH --mem=50G
#SBATCH --partition=short
#SBATCH --cpus-per-task=8

## Genome_Indexing.sh
## Created by Adrianna Vandeuren on 01/25/2024
## This script serves to builf the reference genome you need 

## Activate necessary environment
source /etc/profile.d/modules.sh
module load anaconda3
conda init
source activate RNA-seq

## Define the paths to the reference genome
projPath="/work/daylab/Shannon"
WBcel235Path="${projPath}/CT/referenceGenomes/WBcel235"
WBcel235PathGA="/work/daylab/Shannon/RNA-seq/geneAnnotations"
#refSeqPath="/work/daylab/Shannon/RNA-seq/geneAnnotations"

cd ${WBcel235PathGA}/

## We are going to use the hg19 genome like we do for the CUT&Tag analysis 
## Step 1: get the data from all the chromosomes from UCSC
#wget hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz

#echo "downloaded reference genome"

## Step 2: unzip the files
#tar -zxvf chromFa.tar.gz

echo "unzipped files"

## Step 3: create one big file containing the whole genome, to save space you can delete all the generated files except hg19.fa after this point. 
#cat chr*.fa > hg19.fa

#echo "condense genome to one fasta file"

## Download the GTF file with the gene annotations for your reference genome
## You can find some here: https://www.gencodegenes.org/human/release_19.html
## You need to download the one that matches the origin of your genome file, so here it needs to come from UCSC
#wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.ensGene.gtf.gz

echo "downloaded GTF file for our genome of interest"

## unzip the file 
#gunzip hg19.ensGene.gtf.gz

echo "unzipped GTF file" 

## Use STAR to index the reference genome since that is the tool we are going to use for the alignment of our reads too
STAR --runMode genomeGenerate --genomeDir ${WBcel235Path} --genomeFastaFiles ${WBcel235Path}/GCF_000002985.6_WBcel235_genomic.fna --sjdbGTFfile ${WBcel235Path}/WBcel235genomic.gtf


