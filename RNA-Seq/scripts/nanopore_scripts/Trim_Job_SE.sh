#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=Trim_Job_SE
#SBATCH --time=02:00:00
#SBATCH --mem=10G
#SBATCH --partition=short
#SBATCH --output=/projects/daylab/yulia/RNA-Seq/slurm_logs/%x_%j.out
#SBATCH --error=/projects/daylab/yulia/RNA-Seq/slurm_logs/%x_%j.err

# NOTE : Modify sbatch error and log file paths above for you

# Created by Julianne (Yulia) Murthy on 06/24/2025
# Do quality trimming of poor quality bases with filtlong

# So cutadapt can catch errors if fails
set -euo pipefail

# Activate conda environment 
source /etc/profile.d/modules.sh
# use this step if you do not have filtlong in your env 
module load anaconda3
eval "$(conda shell.bash hook)"
conda activate RNA-seq

# The paths are numbered since the actual paths are given in the Launch_Trim_SE.sh in this order 
fastqFile="$1"
outFile="$2"

# error handlng: make sure output directory exists 
mkdir -p "$(dirname "$outFile")"

# Use filtlong to filter the fastq files
filtlong \
  --min_length 50 \
  --min_mean_q 7 \
  --mean_q_weight 2.0 \
  "$fastqFile" | gzip > "$outFile"
