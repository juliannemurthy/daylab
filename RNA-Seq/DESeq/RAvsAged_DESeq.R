## Script for featureCounts and DESeq2 analysis for glp-1 Reproductive Adult (RA) vs Aged Illumina RNA-seq alignments
## Written by Julianne (Yulia) Murthy on 7/18/2025

# Install BiocManager and subread
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsubread", dependencies = TRUE, force = TRUE)
BiocManager::install("DESeq2", dependencies = TRUE, force = TRUE)
BiocManager::install("vsn", dependencies = TRUE, force = TRUE)
BiocManager::install("EnhancedVolcano", dependencies = TRUE, force = TRUE)
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("pheatmap") 
install.packages("stringr")

# Load necessary libraries
library(Rsubread)
library(DESeq2)
library(ggplot2)
library(vsn)
library(ggrepel)
library(pheatmap)
library(stringr)
library(EnhancedVolcano)

# set paths for BAM and gtf files
bam.path <- "/Users/juliannemurthy/RAvsAged_sortedbam"
bam.files <- list.files(bam.path, pattern = "\\.bam$", full.names = TRUE)
gtf.file <- "/Users/juliannemurthy/ce11.ncbiRefSeq.gtf"

# Run featureCounts, assign to fc object
# specify path to BAM files with read alignments to be counted
# path to GTF annotation file
# Indicates GTF file is in GTF not SAF format
# exon -> count reads mapped to exons
# "gene_id" -> extract gene id attribute to group exon-level counts under each gene
# useMetaFeatures = TRUE -> summarize counts at gene level than exon level
# allowMultiOverlap = FALSE -> discard reads that overlap more than one feature, such as overlapping genes or exons
# isPairedEnd = TRUE -> reads are paired end
# strandSpecific = 0 -> strandedness, these are unstranded as confirmed by 50/50 strandednress split confirmed by rseqc 
# primaryOnly = TRUE -> use only primary alignments, not secondary or supplementary in counts
# 8 threads for parallel processing
fc <- featureCounts(
  files = bam.files,
  annot.ext = gtf.file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE,
  allowMultiOverlap = FALSE,
  isPairedEnd = TRUE,
  strandSpecific = 0,
  primaryOnly = TRUE,
  nthreads = 8,
)
# Assign the gene level raw counts matrix to "counts"
counts <- fc$counts
# Rename columns to bam file names
colnames(counts) <- basename(bam.files)
# Inspect the counts data -> genes are rows and BAM file samples columns
head(counts)

# Write results to csv file, ensures gene names (row names) go into the first column. first row are sample names
write.csv(counts, file = "RAvsAged_counts.csv")

# Read in the counts data, row.names = 1 says genes are in first column
countsdata <- read.csv("RAvsAged_counts.csv", header = TRUE, row.names = 1)
# Check the counts data
head(countsdata)

# Filter counts to reduce noise and improve statistical power for DE
# Use <=50 for Illumina, nanopore data has lower read depth, higher variability, greater length per read
# keep genes with at least 10 counts in at least 2 samples, better than using overall sum cutoffs
countsdata <- countsdata[rowSums(countsdata >= 10) >= 2, ]
# Check the filtered counts table
head(countsdata)

# Set up the conditions of the experiment -> L1 and adult, MUST MATCH order of columns in countsdata
condition <- factor(c("Aged", "Aged", "Aged", "RA", "RA", "RA"))
# Check the condition factor
condition  

# Generate a table that matches the samples to the conditions for DESeq2
coldata <- data.frame(row.names = colnames(countsdata), condition = condition)
# Check the coldata table
coldata  

# Generate a DESeq object for analysis, design = ~condition models gene expression as a function of life stage
dds <- DESeqDataSetFromMatrix(countData = countsdata, colData = coldata, design = ~condition)

# Run DESeq analysis
# estimates size factors, normalizes counts
# estimates dispersion, fit negative binomial GLMs, perform Wald tests
dds <- DESeq(dds)

# VST normalization (Variance Stabilizing Transformation) to normalize counts for visualization later
# blind = FALSE preserves experimental design structure in transformation
vsdata <- vst(dds, blind = FALSE)

# PCA (Principal Component Analysis) _> visualize sample clustering by condition using top genes
# Should show separation between two groups of samples if strong signal
pca <- plotPCA(vsdata, intgroup = "condition")  + ggtitle("PCA of glp-1 Reproductive Adult vs Aged (VST-normalized)")
pca
ggsave("RAvsAged_PCA_plot.png", pca, width = 6, height = 5)
# Dispersion vs mean expression plot -> check fit and detect over dispersion or outlier for QC
png("RAvsAged_Dispersion_plot.png", width = 6, height = 5, units = "in", res = 300)
plotDispEsts(dds)
dev.off()

# Extract differential expression results 
res <- results(dds, contrast = c("condition", "Aged", "RA"), test = "Wald", alpha = 0.05)
res

# Extract log 2 fold changes (LFC) for better volcano plots and heatmaps and ranking purposes
# If the shrinkage estimator apeglm is used in published research, please cite: Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics. 10.1093/bioinformatics/bty895
# Run resultsNames(dds) to get the coef. For example, "condition_RA_vs_Aged" means log2 fold change of RA relative to Aged. 
resLFC <- lfcShrink(dds, coef = "condition_RA_vs_Aged", type="apeglm")
resLFC

# Filter significant genes by adjusted p-value
sigs <- na.omit(resLFC) # remove genes without valid p-values, low count
sigs <- sigs[sigs$padj < 0.05, ]
# Check the significant results
sigs  

# Save only significant genes sigs to CSV 
write.csv(sigs, file = "DESeq2_RAvsAged_significant_genes.csv")

# Plot heatmap based on the 30 most significant (p<0.05) DEGs. Very statistically significant but may have small effect size
topgenes <- head(order(sigs$padj),30)
mat <- assay(vsdata)[rownames(sigs)[topgenes], ]
# Average replicates for each condition
condition_labels <- colData(vsdata)$condition
RA_mean <- rowMeans(mat[, condition_labels == "RA"])
Aged_mean <- rowMeans(mat[, condition_labels == "Aged"])
mat_avg <- cbind(RA = RA_mean, Aged = Aged_mean)
# Center the data by subtracting row means, tells you how much each condition deviates from gene's avg expression
mat_avg <- mat_avg - rowMeans(mat_avg)
#Annotation for the heatmap
annotation_avg <- data.frame(condition = c("RA", "Aged"))
rownames(annotation_avg) <- c("RA", "Aged")
#Plot heatmap  
heatmap <- pheatmap(mat_avg,
         annotation_col = annotation_avg, 
         main = "Top 30 DE Genes Heatmap (averaged): glp-1 RA vs Aged",
         fontsize_col = 10
         )
heatmap
ggsave("RAvsAged_Heatmap.png", heatmap, width = 8, height = 6)

# Volcano plot creation
# Select the 30 genes with smallest padj, regardless of significance threshold -> Top 30 by |log2FC| × –log10(padj)
top30resLFC=head(order(resLFC$padj),30)
# x-axis: large absolute value of log 2 fold change -> genes with strong biological effect size
# y-axis: –log10(pvalue or padj) -> measures statistical significance
volcano <- EnhancedVolcano(resLFC,
                           lab = rownames(resLFC),
                           selectLab = rownames(resLFC)[top30resLFC],
                           x = "log2FoldChange",
                           y = "pvalue",
                           title = "Volcano Plot: glp-1 RA vs Aged",
                           pointSize = 1.5,
                           labSize = 3,
                           colAlpha = 0.8,
                           max.overlaps = Inf 
)
volcano
ggsave("RAvsAged_Volcano.png", volcano, width = 8, height = 6)

### Important to note that the top 30 genes shown in heatmap vs volcano plot are different


# Save plots
# one way; save dat.long to a variable ex: p<- dat.long
# use function ggsave, 1st parameter is the file variable p, and then define filename + specify dimensions
#ggsave(p, filename = 'heatmap_save1.pdf', width = 10, height = 8)

# second way to save plot is typing this before the command dat.long etc
# pdf("heatmap_save2.pdf", width = 10, height = 8)
# must write the below after the command is done as well
# dev.off()

# Save session infor at the end
writeLines(capture.output(sessionInfo()), "RAvsAged_session_info.txt")

nrow(sigs)
nrow(resLFC)

# Save DESeq results in TSV file for BETA analysis
resLFC_beta <- data.frame(
  genename = rownames(resLFC),
  log2FC = resLFC$log2FoldChange,
  pval = resLFC$pvalue
)
resLFC_beta <- na.omit(resLFC_beta)
resLFC_beta
write.table(resLFC_beta,
            file = "RAvsAged_DESeq.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)


