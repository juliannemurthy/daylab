#!/bin/bash
#SBATCH --job-name=BAM_Indexing
#SBATCH --time=15:00:00
#SBATCH --mem=50G
#SBATCH --partition=short
#SBATCH --cpus-per-task=10
#SBATCH --output=/projects/daylab/yulia/G4_Peak_Diff_Analysis/experiments/L1vsRA/slurm_logs/%x_%j.out
#SBATCH --error=/projects/daylab/yulia/G4_Peak_Diff_Analysis/experiments/L1vsRA/slurm_logs/%x_%j.err

# NOTE: adjust the slurm log paths above for your usage 

# Written by Julianne (Yulia) Murthy on 08/21/2025
# Script creates indexes of aligned paired-end BAM files for required use in DiffBind chip-seq peak differential analysis

# Load module
module load samtools/1.21

# Set path to bams
bamPath="/projects/daylab/yulia/G4_Peak_Diff_Analysis/experiments/L1vsRA/infiles/bamfiles"

# Index the bam files
for bam in "$bamPath"/*.bam; do
	samtools index -@ "${SLURM_CPUS_PER_TASK:-1}" $bam
done
