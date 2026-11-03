#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=FASTQC_Job
#SBATCH --time=10:00:00
#SBATCH --mem=10G
#SBATCH --partition=short
#SBATCH --output=/projects/daylab/yulia/RNA-Seq/slurm_logs/%x_%j.out
#SBATCH --error=/projects/daylab/yulia/RNA-Seq/slurm_logs/%x_%j.err

# NOTE : Modify sbatch error and log file paths above for you

## FASTQC_Job.sh 
## Modifed by Julianne (Yulia) Murthy on 07/09/2025
## Created by Adrianna Vandeuren 08/02/2023

# Safety error catch
set -euo pipefail

# Activate conda environment
source /etc/profile.d/modules.sh
module load fastqc/0.12.1 OpenJDK/22.0.2 perl/5.40.0

# The paths are numbered since the actual paths are given in the Launch_FASTQC.sh in this order 
fastqPath="$1"
outputPath="$2"
inputFile="$3"

# Run FASTQC, it will output results to output path 
# Check the quality by clicking on the output html in OOD. Also copy paste the adapter sequence there for later
fastqc -o $outputPath -f fastq $fastqPath/$inputFile

