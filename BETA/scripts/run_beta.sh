#!/bin/bash

## run_beta.sh
## Created by Ella Chee on 10/23/24.
## This script batches BETA analysis with RNA-seq/BED files corersponding to the given timepoint
## Usage: ./run-beta.sh 01hr (timepoint)

module load python/2.7.15
module load R/4.4.0

timepoint=$1

mkdir -p ./beta-plus-0.05-unique/${timepoint}
sbatch --job-name="run_beta_${fasta_name}.sh" --nodes="1" --time="03:00:00" --mem="10G" --wrap="BETA plus -p cut-and-tag-unique-peaks/converted-breaks/A549_unique_${timepoint}_all.bed -e rna-seq-data/tsv/A549_UV_timecourse_RNA_seq_${timepoint}.tsv --info 3,5,8 -g hg19 -c .05 -n ${timepoint} --gname2 -k O -o ./beta-plus-0.05-unique/${timepoint} --gs hg19.fa"

echo 'timepoint:' $timepoint 'batched'