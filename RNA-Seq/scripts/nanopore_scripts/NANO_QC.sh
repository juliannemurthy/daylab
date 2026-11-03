#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=NANO_QC
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=14
#SBATCH --mem=10G
#SBATCH --partition=short
#SBATCH --output=/projects/daylab/yulia/RNA-Seq/slurm_logs/%x_%j.out
#SBATCH --error=/projects/daylab/yulia/RNA-Seq/slurm_logs/%x_%j.err

# NOTE : Modify sbatch error and log file paths above for you

## NANO_QC.sh 
## Modifed by Julianne (Yulia) Murthy on 06/24/2025
## Created by Adrianna Vandeuren 08/02/2023
# Designed for use with Oxford Nanopore RNA-seq data

# Activate conda environment, download new packages as needed
source /etc/profile.d/modules.sh
# module load anaconda3
# source activate RNA-seq
# pip install nanoplot
# pip install nanostat
# pip install "plotly>=6.1.1"
# pip install "kaleido>=0.2.1"

# The paths are numbered since the actual paths are given in the Launch_NANO_QC.sh in this order 
fastqPath="$1"
outputPath="$2"
inputFile="$3"

# Extract base name from fastq
base=$(basename "$inputFile" .fastq.gz)

# error handlng: make sure output directory exists 
mkdir -p "$outputPath/$base"

# Run NanoPlot to determine read lengths, and get QC plots. Add quality filters simialr to what will be run with Filtlong later
	 # --plots kde hex: Include KDE and hexbin plots for quality vs length, etc.
	 # --format: Output plot formats
	 # --N50 : Display N50 on the length histogram
	 # Show log-scaled length histogram, recommended if long tails skew plots
	 # Filter out reads shorter than 100 bp
	 # Filter out reads with average Q-score below 7
	 # 	Higher resolution plots
	 # Custom plot/report title
	 # Embed command info and settings in HTML report
NanoPlot \
  --fastq "$fastqPath/$inputFile" \
  -o "$outputPath/$base" \
  -t 8 \
  --prefix "${base}_nanoplot" \
  --plots kde hex \
  --N50 \
  --loglength \
  --minlength 100 \
  --minqual 7 \
  --dpi 150 \
  --title "NanoPlot QC: ${base}" \
  --info_in_report

# Run Nanostat to get QC summary stats. Note will not have quality filters like above so stats can vary from above slightly
NanoStat \
  --fastq "$fastqPath/$inputFile" \
  --name "$outputPath/$base/${base}_nanostat.tsv" \
  --prefix "${base}_nanostat" \
  --threads 4 \
  --tsv
