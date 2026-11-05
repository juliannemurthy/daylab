# Script for G4 peaks differential analysis between L1 and RA groups
# Goal is to quantify peak differences not just as presence/absence of of a peak in both conditions, but increases or decreases of a peak
# Written by Julianne (Yulia) Murthy on 10/28/2025

# Install packages
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DiffBind")
BiocManager::install("DESeq2", dependencies = TRUE, force = TRUE)
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("pheatmap") 
BiocManager::install("edgeR")

# Load necessary libraries
library(BiocParallel)
register(MulticoreParam(workers=4))
library(DiffBind)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(edgeR)

# set paths for input files
# BAM files must be coordinate sorted and indexed using samtools index .bam
bam_dir <- "/projects/daylab/yulia/G4_Peak_Diff_Analysis/experiments/L1vsRA/infiles/bamfiles"
L1_bamfiles <- file.path(bam_dir, c("L1-1_Aligned.sortedByCoord.out.bam", "L1-2_Aligned.sortedByCoord.out.bam", "L1-3_Aligned.sortedByCoord.out.bam"))
RA_bamfiles <- file.path(bam_dir, c("RA-1_Aligned.sortedByCoord.out.bam", "RA-2_Aligned.sortedByCoord.out.bam", "RA-3_Aligned.sortedByCoord.out.bam"))
# Individual narrowPeak bed files per rep are included here
L1_peak_beds   <- c("/projects/daylab/yulia/G4_Peak_Diff_Analysis/experiments/L1vsRA/infiles/peakfiles/L1_BG4_2.narrowPeak",
                    "/projects/daylab/yulia/G4_Peak_Diff_Analysis/experiments/L1vsRA/infiles/peakfiles/L1_BG4_3.narrowPeak")
RA_peak_beds <- c("/projects/daylab/yulia/G4_Peak_Diff_Analysis/experiments/L1vsRA/infiles/peakfiles/RA_BG4_5.narrowPeak",
                    "/projects/daylab/yulia/G4_Peak_Diff_Analysis/experiments/L1vsRA/infiles/peakfiles/RA_BG4_6.narrowPeak")

# Create an output directory if it doesn't already exist + any other subfolders within
outdir <- "DiffBind_L1vsRA"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

# Create a temporary sample sheet in csv format for only 4 sample-peak combos
# NOTE: if BED file numeric chrom score is not in column below, must add using PeakFormat, ScoreCol as detailed in the DiffBind manual
samples_csv <- file.path(outdir,"L1vsRA_samples.csv")
make_sample_line <- function(id, cond, rep, bam, peaks,
tissue="C.elegans",
factor="G4", peakcaller="macs") {
  paste(id, tissue, factor, cond, rep, bam, peaks, peakcaller, sep=",")
}

csv_lines <- c(
  "SampleID,Tissue,Factor,Condition,Replicate,bamReads,Peaks,PeakCaller",
  make_sample_line("L1_1","L1",1,L1_bamfiles[1],L1_peak_beds[1]),
  make_sample_line("L1_2","L1",2,L1_bamfiles[2],L1_peak_beds[2]),
  make_sample_line("RA_1","RA",1,RA_bamfiles[1],RA_peak_beds[1]),
  make_sample_line("RA_2","RA",2,RA_bamfiles[2],RA_peak_beds[2])
)
writeLines(csv_lines, samples_csv)
# Load in that temporary sample sheet
db_peaks <- dba(sampleSheet = samples_csv)

# Build a consensus peakset from these four samples across conditions
db_peaks <- dba.peakset(db_peaks, consensus=DBA_CONDITION, minOverlap=1)

# Retrieve the consensus as GRanges
consensus <- dba.peakset(db_peaks, bRetrieve=TRUE)

# Now build full sample sheet with all 6 BAMs, leave peaks blank
samples_all_csv <- file.path(outdir, "L1vsRA_allBAMs.csv")
make_sample_line <- function(id, cond, rep, bam, peaks,
                             tissue="C.elegans",
                             factor="G4", peakcaller="macs") {
  paste(id, tissue, factor, cond, rep, bam, peaks, peakcaller, sep=",")
}

csv_lines_all <- c(
  "SampleID,Tissue,Factor,Condition,Replicate,bamReads,Peaks,PeakCaller",
  make_sample_line("L1_1","L1",1,L1_bamfiles[1],""),  
  make_sample_line("L1_2","L1",2,L1_bamfiles[2],""),  
  make_sample_line("L1_3","L1",3,L1_bamfiles[3],""), 
  make_sample_line("RA_1","RA",1,RA_bamfiles[1],""),  
  make_sample_line("RA_2","RA",2,RA_bamfiles[2],""), 
  make_sample_line("RA_3","RA",3,RA_bamfiles[3],"")   
)
writeLines(csv_lines_all, samples_all_csv)
# Load in the sample sheet
db <- dba(sampleSheet = samples_all_csv)

# Optional blacklist/greylist for the species if BED available??
# db <- dba.blacklist(db, blacklist= "/path/to/ce11_blacklist.bed", greylist=FALSE)

head(db$peaks[[1]])


# Count reads across all samples against fixed consensus peakset
# Set standardized widths via summits with narrow marks like H3K4me3
# Note: Diffbind 3.0 will automatically re-center peaks around the point of maximum read pileup (summit=)
# If not specifying summit default is 200, using 250 since it was used in one of the manuals and I wanted more wiggle room
#### For G4 structures are peaks sharp like TFs or broader???
# Requires a peak to show up in at least one rep, runs in parallel to speed up

# NOTE: THIS STEP REQUIRES A LOT OF MEMORY! RUN ON OOD RSTUDIO ON CLUSTER AND REQUEST 75 GB OF MEMORY AND BE PREPARED TO WAIT
db <- dba.count(db, summit=250, peaks=consensus, bParallel=TRUE) 
db

# Save progress checkpoint
saveRDS(db, file.path(outdir,"db_counted.rds"))
db <- readRDS(file.path(outdir, "db_counted.rds"))

# Run correlation and pca plots to QC before modeling data
png(file.path(outdir,"pre_correlation_L1vsRA.png"));  plot(db); dev.off()
png(file.path(outdir,"pre_pca_L1vsRA.png")); dba.plotPCA(db, DBA_CONDITION, label=DBA_ID); dev.off()

# Define the contrasts (RA vs Aged), only set up where each condition has at least 2 replicates for statistical significance
db <- dba.contrast(db, categories=DBA_CONDITION, minMembers=2)

# Normalize with TMM normalization
tmm <- dba.normalize(db, method=DBA_EDGER, normalize=DBA_NORM_TMM)

# Run Differential Analysis on the peaks
db_tmm <- dba.analyze(tmm, method=DBA_EDGER, bBlacklist=FALSE, bGreylist=FALSE)
db_tmm

# MA plots visualize which data points are differentially bound. Each of these points represents a binding site and red points indicate differentailly bound ones
png(file.path(outdir,"MA_DESeq2_L1vsRA_TMMNORM.png"));     dba.plotMA(db_tmm, method=DBA_EDGER); dev.off()

# Volcano plot to show the significantly differentially bound sites and show their fold enrichment.
png(file.path(outdir,"Volcano_DESeq2_L1vsRA_TMMNORM.png")); dba.plotVolcano(db_tmm, method=DBA_EDGER); dev.off()

# Write differential peaks to csv file data frame
# Saves threshold FDR/adjusted p-value ≤ 0.05 and includes normalized read counts per sample (binding matrix values) for each reported peak
res_DE <- dba.report(db_tmm, method=DBA_EDGER, th=0.05, bCounts=TRUE)
write.csv(as.data.frame(res_DE), file.path(outdir,"DE_peaks_L1vsRA_TMMNORM.csv"), row.names=FALSE)

# Run correlation and pca plots QC post-normalization
png(file.path(outdir,"post_correlation_L1vsRA.png")); plot(db); dev.off()
png(file.path(outdir,"post_pca_L1vsRA.png")); dba.plotPCA(db, DBA_CONDITION, label=DBA_ID); dev.off()

