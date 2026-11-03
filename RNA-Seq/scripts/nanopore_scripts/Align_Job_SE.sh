#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=Align_Job_SE
#SBATCH --time=12:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=18
#SBATCH --partition=short
#SBATCH --output=/projects/daylab/yulia/RNA-Seq/slurm_logs/%x_%j.out
#SBATCH --error=/projects/daylab/yulia/RNA-Seq/slurm_logs/%x_%j.err

# NOTE : Modify sbatch error and log file paths above for you

# Created by Julianne (Yulia) Murthy on 06/24/2025

# Safety error catch
set -euo pipefail

# Activate conda environment 
source /etc/profile.d/modules.sh
# use these step if you do not have minimap2 in your RNA-seq env 
# conda activate RNA-seq
# conda install -c bioconda minimap2
module load anaconda3
eval "$(conda shell.bash hook)"
conda activate RNA-seq

# Provide input and output paths
projPath="/projects/daylab/yulia/RNA-Seq"
genomeDir="$projPath/reference/ref/RefSeq/GCF_000002985.6_WBcel235_genomic.fna"
trimmedPath="$projPath/trimmed/RoachEtAl_trimmed"
sortedBAMPath="$projPath/aligned/minimap2_Align/RoachEtAl_sortedbam"

# error handling: Make output path if it doesnt not already exist
mkdir -p $sortedBAMPath

# $1 refers to the first input written in the Launch script -> $fastqPath/$file
fastqFile="$1"
# Obtain only the base of the file names by stripping .fastq.gz
base=$(basename $fastqFile .fastq.gz)
# Output BAM path
outputBAM="$sortedBAMPath/${base}.bam"

# Log process
echo "Aligning: $fastqFile → $outputBAM"

# Performs minimap2 alignment with parameters:
	# -ax splice : preset for spliced alignment
	# -k14 : minimizer k-mer size, smaller value is better for noisy reads like ONT
	# -u : long read preset
	# --secondary=no : only output primary alignments, suitable for featureCounts
	# -G 25000: max intron length, common for C. elegans
	# -t 8 : use 8 threads 
# convert SAM to BAM
# sort the BAM files
minimap2 -ax splice -uf -k14 --secondary=no -G 25000 -t 8 "$genomeDir" "$fastqFile" | \
samtools view -F 2048 -@ 4 -b | \
samtools sort -@ 4 -o "$outputBAM"

# Index the sorted BAM files for formatted use in featurecounts later
samtools index "$outputBAM"   

