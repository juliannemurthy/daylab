#!/bin/bash

## run_beta_ce11.sh
## Created by Ella Chee on 10/23/24.
## Modified by Julianne (Yulia) Murthy on 07/22/25
## This script batches BETA analysis with RNA-seq/BED files corersponding to the given timepoint
## Usage: bash ./run_beta_ce11.sh RA (timepoint)

# Load R and beta_env conda environment required to run BETA
module load R/4.4.1
# source /etc/profile.d/modules.sh
# module load anaconda3
# To replicate this conda env run: conda env create -f /projects/daylab/yulia/BETA/env/beta_env.yml
# Contains python 2.7, argparse, numpy, setuptools, etc
# conda activate beta_env 

# Takes a timepoint argument so it can locate correct input files and name output files
# run example: bash run_beta_ce11.sh RA
timepoint=$1

# Define project path
projPath="/projects/daylab/yulia/BETA"

# Makes output dir for BETA run with timepoint name if not already made
mkdir -p $projPath/betaplus_outfiles/RAvsAged_beta/${timepoint}

# sbatch slurm job submission
# BETA execution -> --wrap:
	# BETA plus : runs in plus mode, TF active and repressive function prediction, direct targets detecting and motif analysis in target regions
	# -p : path to BED peaks file of G4 binding sites 
		# do not contain the header in the bed file, and make sure it is tab delimited
	# -e : RNA-seq differential expression file TSV format
		# differential expression file should contain all the genes in the genome, BETA will use all the info to get the static genes, and isolate the up regulated genes and down regulated genes based on the threshold you input
		# Make sure your differential expression file do not have the header or add the ‘#’ in the front of your header line
	# --info 3,5,8 : which columns are gene ID, log2FC, and p-value respectively
	# -g : specifies the genome of your data, hg19 for human and mm9 for mouse, others, ignore this one
	# -c .05 : Input a value 0–1 as a threshold of one-tail Kolmogorov-Smirnov test to determine significant difference
	# -da : example, 500 means get top 500 most significant expression changed genes of up and down, default 0.5
	# -n : specifies the prefix of the output files
	# --gname2: add if your gene ID is official gene symbol (not Entrez IDs), ensures correct gene name mapping
	# -k O : specifies the format of the differential expression data—LIM, CUF, BSF or O for Other
	# -d DISTANCE, --distance DISTANCE : Set a number which unit is 'base'. It will get peaks within this distance from gene TSS. DEFAULT=100000 (100kb)
	# -o : output dir for the run 
	# -r : UCSC annotation file refGene.txt !!! input if you have genome not hg19 nor mm9
	# --gs : MUST BE UCSC -> FASTA file of whole genome sequence for motif scan

sbatch \
--job-name="run_beta_ce11_${timepoint}.sh" \
--output=$projPath/slurm_logs/run_beta_ce11_${timepoint}_--daALL.out \
--error=$projPath/slurm_logs/run_beta_ce11_${timepoint}_--daALL.err \
--nodes=1 \
--time="03:00:00" \
--mem="10G" \
--wrap="source /shared/EL9/explorer/anaconda3/2024.06/etc/profile.d/conda.sh && conda activate beta_env && \
BETA plus \
	-p $projPath/betaplus_infiles/RAvsAged_beta/${timepoint}.bed \
	-e $projPath/betaplus_infiles/RAvsAged_beta/RAvsAged_DESeq.tsv \
	--info 1,2,3 \
	-c 0.1 \
	-n ${timepoint} \
	--gname2 \
	-k O \
	-o $projPath/betaplus_outfiles/RAvsAged_beta/${timepoint}_--daALL \
	-r $projPath/reference/refGene_ce11_BETA.bed \
	--gs $projPath/reference/ce11_fixed.fa \
	--da 5619 "
	
# BETA output contents include
	# Target Gene Associated Motif Analysis .html sumamry report
	# up/downtargeted genes in a TXT file containing location, refseqID, rank product,strand, gene symbol 
	# BED file containing same cols and score
	# PDF of activating/repressive function prediction
	# R script to generate PDF figure
	
# BETA Optional Parameters
	# -n NAME, --name NAME : This argument is used to name the result file.
	# -o OUTPUT, --output OUTPUT : The directory to store all the output files.
	# --gname2 : If this switch is on, gene or transcript IDs in files given through -e will be considered as official gene symbols, DEFAULT=FALSE
	# --info EXPREINFO : Specify the geneID, up/down status and statistical values column of your expression data. DEFAULT:2,5,7 for LIMMA; 2,10,13 for Cuffdiff and 1,2,3 for BETA specific format
	# --pn PEAKNUMBER : The number of peaks will be contribute to the regulatory potential score. DEFAULT=10000
	# -d DISTANCE, --distance DISTANCE : Set a number which unit is 'base'. It will get peaks within this distance from gene TSS. DEFAULT=100000 (100kb)
	# --df DIFF_FDR : Input a number 0~1 as a threshold to pick out the most significant differential expressed genes by FDR or other statistical values, DEFAULT = 1,that is select all the genes
	# --da DIFF_AMOUNT : Get the most significant differential expressed genes by the percentage(0-1) or number (larger than 1). Genes will be ranked by the statistical values. For example, 2000, BETA will set top 2000 genes of up and down as the differentially expressed genes. DEFAULT = 0.5, that is to get top 50 percent genes of up and down separately
	# -c CUTOFF, --cutoff CUTOFF : Input a number between 0~1 as a threshold of One tail KS test to determine if two datasets differ significantly, default=1e-3
	# -r REFERENCE, --reference REFERENCE : the refgene info file downloaded from UCSC genome browser.input this file only if your genome is neither hg19 nor mm9
	# --bl BOUNDARY LIMIT : Boolean Value. Whether or not use CTCF boundary to get a peak’s associated gene, DEFAULT=FALSE
	# --bf BOUNDARYFILE : Some BED format boundary file, use this only when You set --bl and the genome is neither hg19 nor mm9
	
echo 'timepoint:' $timepoint 'batched'
