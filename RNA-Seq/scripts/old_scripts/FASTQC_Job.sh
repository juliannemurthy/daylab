#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=FASTQC
#SBATCH --time=10:00:00
#SBATCH --mem=10G
#SBATCH --partition=short

## FASTQC_Job.sh 
## Created by Adrianna Vandeuren 08/02/2023
source /etc/profile.d/modules.sh
projPath="/work/daylab/Shannon/RNA-seq/30-1172669254"

## Depending on how you load picard and your server environment, the picardCMD can be different. Adjust accordingly.
## picardCMD="java -jar picard.jar"
## We have picard in the conda environment called Java
module load anaconda3
conda init bash
source activate RNA-seq


fastqc -o ${projPath}/samples/FASTQC/${NAME} -f fastq ${R1}
fastqc -o ${projPath}/samples/FASTQC/${NAME} -f fastq ${R2}

