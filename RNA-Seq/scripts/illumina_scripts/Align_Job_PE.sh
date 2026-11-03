#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=Align_Job_PE
#SBATCH --time=12:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=8
#SBATCH --partition=short
#SBATCH --output=/projects/daylab/yulia/RNA-Seq/slurm_logs/%x_%j.out
#SBATCH --error=/projects/daylab/yulia/RNA-Seq/slurm_logs/%x_%j.err

# NOTE : Modify sbatch error and log file paths above for you

# Modified by Julianne (Yulia) Murthy on 07/09/2025
# Created by Adrianna Vandeuren 08/02/2023

# Safety error catch
set -euo pipefail

# Activate conda environment 
source /etc/profile.d/modules.sh
module load anaconda3
eval "$(conda shell.bash hook)"
conda activate RNA-seq

# Provide input and output paths
projPath="/projects/daylab/yulia/RNA-Seq"
genomeDir="$projPath/reference/index_UCSC"
alignPath="$projPath/aligned/STAR_Aligned/L1vsRA_sortedbam"

# Make aligned path if it doesn't exist
mkdir -p "$alignPath"

# Refer to inputs defined in Launch script
R1="$1"
R2="$2"
base="$3"

# Run STAR
	# uses indexed genome files
	# handles .gz files
	# The output files will be SORTED BAM files by coordinate compatible for featureCounts later
STAR --runThreadN 8 \
  --genomeDir "$genomeDir" \
  --readFilesIn "$R1" "$R2" \
  --readFilesCommand zcat \
  --outFileNamePrefix "$alignPath/${base}_" \
  --outSAMtype BAM SortedByCoordinate 

  
