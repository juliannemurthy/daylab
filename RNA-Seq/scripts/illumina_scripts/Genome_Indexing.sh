#!/bin/bash
#SBATCH --job-name=Genome_Indexing
#SBATCH --time=15:00:00
#SBATCH --mem=50G
#SBATCH --partition=short
#SBATCH --cpus-per-task=10
#SBATCH --output=/projects/daylab/yulia/RNA-Seq/slurm_logs/%x_%j.out
#SBATCH --error=/projects/daylab/yulia/RNA-Seq/slurm_logs/%x_%j.err

# NOTE : Modify sbatch error and log file paths above for you

## Genome_Indexing.sh
## Modifed by Julianne (Yulia) Murthy on 07/09/2025
## Created by Adrianna Vandeuren on 01/25/2024
## This script serves to index C.elegans WBcel235 genome with STAR 

# Safety error catch
set -euo pipefail

## Activate conda environment
source /etc/profile.d/modules.sh
module load anaconda3
eval "$(conda shell.bash hook)"
conda activate RNA-seq

# Run below commands on your terminal to make sure GTF + ref genome file formats have matching FASTA chromosome names
	# i.e. are they both RefSeq or both UCSC or different? If different change to same
# Check FASTA chromosome names
# grep ">" genome.fna | cut -d " " -f1
# Check GTF sequence names
# cut -f1 annotation.gtf | grep -v "^#" | sort | uniq

## Define the paths to the reference genome
projPath="/projects/daylab/yulia/RNA-Seq"
genomeDir="${projPath}/reference/index_UCSC"
# For UCSC (compatible with BETA), use ce11.fa and UCSC GTF file 
genomeFasta="${projPath}/reference/ref/UCSC/ce11.fa"
gtfFile="${projPath}/reference/ref/UCSC/ce11.ncbiRefSeq.gtf"

# error handling: make output genomeDir if it does not exist already
mkdir -p $genomeDir

# run STAR alignment
	# STAR uses 8 threads
	# runMode genomeGenerate tells STAR to build genome index
	# genomeDir is where index will be saved
	# genomeFastaFiles is the genome FASTA file to index
	# sjdbGTFfile is the GTF annotation for improved splice junction detection
	# sjdbOverhang is Read Length - 1, specifies the length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database
		# run this command in terminal on .fastq.gz files to determine read length: zcat sample.fastq.gz | awk 'NR%4==2 {print length($1)}' | sort -n | uniq -c
	
STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir $genomeDir \
--genomeFastaFiles $genomeFasta \
--sjdbGTFfile $gtfFile \
--sjdbOverhang 149
