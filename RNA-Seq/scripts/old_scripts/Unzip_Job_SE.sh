#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=Unzip
#SBATCH --time=20:00:00
#SBATCH --mem=50G
#SBATCH --partition=short

## Unzip_Job.sh 
## Created by Adrianna Vandeuren 02/12/2024
## This script will serve to unzip your fastq files


## Activate the environment that containes the tools you need
source /etc/profile.d/modules.sh
module load anaconda3
conda init
source activate RNA-seq

## Unzip files
gunzip ${NAME}*
