## Script for featureCounts and DESeq2 analysis for RoachetAl nanopore long read RNA-seq alignments -> L1 vs adult groups
## Written by Julianne (Yulia) Murthy on 7/2/2025

# Install BiocManager and subread
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsubread", dependencies = TRUE, force = TRUE)
BiocManager::install("DESeq2", dependencies = TRUE, force = TRUE)
install.packages("ggplot2")
install.packages("tidyverse")

# Load necessary libraries
library(Rsubread)
library(DESeq2)
library(ggplot2)
library(tidyverse)


# set paths for BAM and gtf files
bam.path <- "/Users/juliannemurthy/RoachEtAl_sortedbam"
bam.files <- list.files(bam.path, pattern = "\\.bam$", full.names = TRUE)
gtf.file <- "/Users/juliannemurthy/GCF_000002985.6_WBcel235_genomic.gtf"

# Run featureCounts, assign to fc object
  # specify path to BAM files with read alignments to be counted
  # path to GTF annotation file
  # Indicates GTF file is in GTF not SAF format
  # exon -> count reads mapped to exons
  # "gene_id" -> extract gene id attribute to group exon-level counts under each gene
  # useMetaFeatures = TRUE -> summarize counts at gene level than exon level
  # allowMultiOverlap = FALSE -> discard reads that overlap more than one feature, such as overlapping genes or exons
  # isPairedEnd = FALSE -> reads are single end not paired
  # strandSpecific = 1 -> strandedness, these are forward stranded as confirmed by rseqc and RoachetAl
  # primaryOnly = TRUE -> use only primary alignments, not secondary or supplementary in counts
  # 8 threads for parallel processing
  # largestOverlap = TRUE -> used with allowMultiOverlap = FALSE, 5read is assigned to the feature with largest overlap when read overlaps multiple features
  # isLongRead = TRUE -> for longread Nanopore Pacbio data, sets no maximum to number of M (match/mismatch) operations in CIGAR string of BAM
fc <- featureCounts(
  files = bam.files,
  annot.ext = gtf.file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE,
  allowMultiOverlap = FALSE,
  isPairedEnd = FALSE,
  strandSpecific = 1,
  primaryOnly = TRUE,
  nthreads = 8,
  isLongRead = TRUE
)
# Assign the gene level raw counts matrix to "counts"
counts <- fc$counts
# Rename columns to bam file names
colnames(counts) <- basename(bam.files)
# Inspect the counts data -> genes are rows and BAM file samples columns
head(counts)

# Write results to csv file, ensures gene names (row names) go into the first column. first row are sample names
write.csv(counts, file = "RoachEtAl_counts1.csv")

# Read in the counts data, row.names = 1 says genes are in first column
countsdata <- read.csv("RoachEtAl_counts1.csv", header = TRUE, row.names = 1)
# Check the counts data
head(countsdata)

# Filter counts to reduce noise and improve statistical power for DE
  # Use <=50 for Illumina, nanopore data has lower read depth, higher variability, greater length per read
  # keep genes with at least 10 counts in at least 2 samples, better than using overall sum cutoffs
countsdata <- countsdata[rowSums(countsdata >= 10) >= 2, ]
# Check the filtered counts table
head(countsdata)

# Set up the conditions of the experiment -> L1 and adult, MUST MATCH order of columns in countsdata
condition <- factor(c("adult", "adult", "L1", "L1"))
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
  # Should show separation between L1 and adult samples if strong signal
plotPCA(vsdata, intgroup = "condition")  + ggtitle("PCA of RoachetAl Samples (VST-normalized)")
# Dispersion vs mean expression plot -> check fit and detect over dispersion or outlier for QC
plotDispEsts(dds)

5# Extract differential expression results l1 vs adult
res <- results(dds, contrast = c("condition", "L1", "adult"))
res

# Set a more stringent p-value for res, filter significant genes with p-value
sigs <- na.omit(res) # remove genes without valid p-values, low count
sigs <- sigs[sigs$padj < 0.05, ]
# Check the significant results
sigs  

# Save only significant genes sigs to CSV 
write.csv(sigs, file = "DESeq2_RoachEtAl_significant_genes1.csv")

# Make heatmap -> use dds object (not sigs or res) since you need normalized expression values for visualization

# Get top 30 gene names from significant genes in sigs
top_genes <- rownames(sigs)[order(sigs$padj)][1:30]

# Extract VST-normalized expression values from vsdata for those genes
top_expr <- assay(vsdata)[top_genes, ]

# Convert to long format for ggplot use
top_expr_long <- as.data.frame(top_expr) %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "vst") %>%
  mutate(condition = rep(coldata$condition, each = length(top_genes)))
top_expr_long

# Plot heatmap using geom_tile
ggplot(top_expr_long, aes(x = sample, y = gene, fill = vst)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(title = "RoachetAl Top 30 DE Genes (VST-normalized)", x = "Sample", y = "Gene") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Make volcano plot -> use res object (not sigs or dds) since you want to show all genes and highlight significant ones

# Convert DESeq results to data.frame and add significance column to distinguish DE vs non-DE genes
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  mutate(
    significance = case_when(
      padj < 0.05 & abs(log2FoldChange) > 1 ~ "Significant",
      TRUE ~ "Not significant"
    )
  )

# Plot volcano plot that highlights significant hits
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme(legend.position = "bottom")

# Save plots
# one way; save dat.long to a variable ex: p<- dat.long
# use function ggsave, 1st parameter is the file variable p, and then define filename + specify dimensions
#ggsave(p, filename = 'heatmap_save1.pdf', width = 10, height = 8)

# second way to save plot is typing this before the command dat.long etc
# pdf("heatmap_save2.pdf", width = 10, height = 8)
# must write the below after the command is done as well
# dev.off()

# Save session infor at the end
writeLines(capture.output(sessionInfo()), "RoachetAl_session_info.txt")

nrow(sigs)
table(res_df$significance)

```