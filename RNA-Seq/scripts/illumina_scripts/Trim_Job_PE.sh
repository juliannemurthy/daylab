#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=Trim_Job_PE
#SBATCH --time=02:00:00
#SBATCH --mem=8G
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

# The paths are numbered since the actual paths are given in the Launch_Trim_PE.sh in this order 
R1="$1"
R2="$2"
NAME="$3"
trimmomaticPath="$4"
pairedPath="$5"
unpairedPath="$6"

# Since trimmomatic is installed in the RNA-seq environment we do not have to call the jar file
# The LEADING, TRAILING and MINLEN are recommended through the Azenta analysis guide 
	# ILLUMINACLIP does the adapter trimming
	# LEADING and TRAILING removes bases below quality score 3
	# MINLEN discards reads shorter than 36 nucleotides
java -jar "$trimmomaticPath/trimmomatic.jar" \
PE ${R1} ${R2} \
${pairedPath}/${NAME}_R1_paired.fastq.gz \
${unpairedPath}/${NAME}_R1_unpaired.fastq.gz \
${pairedPath}/${NAME}_R2_paired.fastq.gz \
${unpairedPath}/${NAME}_R2_unpaired.fastq.gz \
ILLUMINACLIP:$trimmomaticPath/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36 
