## Script for featureCounts and DESeq2 analysis for L1 (Larvae) vs RA (Reproductive Adult) Illumina RNA-Seq alignments
## Written by Julianne (Yulia) Murthy on 10/06/2025

# Install BiocManager and subread
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

renv::install("bioc::DESeq2")
renv::install("bioc::apeglm")
renv::install("bioc::Rsubread")
renv::install("bioc::vsn")
renv::install("bioc::EnhancedVolcano")
renv::install("ggplot2@3.5.2")
renv::snapshot()

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
bam.path <- "/projects/daylab/yulia/RNA-Seq/aligned/STAR_Aligned/L1vsRA_sortedbam"
bam.files <- list.files(bam.path, pattern = "\\.bam$", full.names = TRUE)
gtf.file <- "/projects/daylab/yulia/RNA-Seq/reference/ref/UCSC/ce11.ncbiRefSeq.gtf"

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
write.csv(counts, file = "L1vsRA_counts.csv")

# Read in the counts data, row.names = 1 says genes are in first column
countsdata <- read.csv("/projects/daylab/yulia/RNA-Seq/DESeq/DESeq_Results_L1vsRA/L1vsRA_counts.csv", header = TRUE, row.names = 1)
# Check the counts data
head(countsdata)

# Filter counts to reduce noise and improve statistical power for DE
# Use <=50 for Illumina, nanopore data has lower read depth, higher variability, greater length per read
# keep genes with at least 10 counts in at least 3 samples (smallestGroupsize), better than using overall sum cutoffs
countsdata <- countsdata[rowSums(countsdata >= 10) >= 3, ]
# Check the filtered counts table
head(countsdata)

# Set up the conditions of the experiment -> L1 and RA, MUST MATCH order of columns in countsdata
condition <- factor(c("L1", "L1", "L1", "RA", "RA", "RA"))
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
pca <- plotPCA(vsdata, intgroup = "condition")  + ggtitle("PCA of L1 vs Reproductive Adult (VST-normalized)")
pca
ggsave("L1vsRA_PCA_plot.png", pca, width = 6, height = 5)
# Dispersion vs mean expression plot -> check fit and detect over dispersion or outlier for QC
png("L1vsRA_Dispersion_plot.png", width = 6, height = 5, units = "in", res = 300)
plotDispEsts(dds)
dev.off()

# Extract differential expression results 
res <- results(dds, contrast = c("condition", "L1", "RA"), test = "Wald", alpha = 0.05)
res

# Extract log 2 fold changes (LFC) for better volcano plots and heatmaps and ranking purposes
# If the shrinkage estimator apeglm is used in published research, please cite: Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics. 10.1093/bioinformatics/bty895
# Run resultsNames(dds) to get the coef. For example, "condition_RA_vs_L1" means log2 fold change of RA relative to L1. 
resLFC <- lfcShrink(dds, coef = "condition_RA_vs_L1", type="apeglm")
#resLFC <- results(dds, contrast = c("condition", "L1", "RA"))
resLFC

#Save resLFC
res_LFC1 <- na.omit(resLFC)
write.csv(res_LFC1, file = "/projects/daylab/yulia/RNA-Seq/DESeq/DESeq_Results_L1vsRA/DESeq2_L1vsRA.csv")

# Filter significant genes by adjusted p-value
sigs <- na.omit(resLFC) # remove genes without valid p-values, low count
sigs <- sigs[sigs$padj < 0.05, ]
# Check the significant results
sigs  

# Save only significant genes sigs to CSV 
write.csv(sigs, file = "DESeq2_L1vsRA_significant_genes.csv")

# Plot heatmap based on the 30 most significant (p<0.05) DEGs. Very statistically significant but may have small effect size
topgenes <- head(order(sigs$padj),30)
mat <- assay(vsdata)[rownames(sigs)[topgenes], ]
# Average replicates for each condition
condition_labels <- colData(vsdata)$condition
L1_mean <- rowMeans(mat[, condition_labels == "L1"])
RA_mean <- rowMeans(mat[, condition_labels == "RA"])
mat_avg <- cbind(L1 = L1_mean, RA = RA_mean)
# Center the data by subtracting row means, tells you how much each condition deviates from gene's avg expression
mat_avg <- mat_avg - rowMeans(mat_avg)
#Annotation for the heatmap
annotation_avg <- data.frame(condition = c("L1", "RA"))
rownames(annotation_avg) <- c("L1", "RA")
#Plot heatmap  
heatmap <- pheatmap(mat_avg,
                    annotation_col = annotation_avg, 
                    main = "Top 30 DE Genes Heatmap (averaged): L1 vs RA",
                    fontsize_col = 10,
                    filename = "L1vsRA_Heatmap.png",   # writes directly
                    width = 8, height = 6
)
heatmap


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
                           title = "Volcano Plot: L1 vs RA",
                           pointSize = 1.5,
                           labSize = 3,
                           colAlpha = 0.8,
                           max.overlaps = Inf 
)
volcano
ggsave("L1vsRA_Volcano.png", volcano, width = 8, height = 6)

### Important to note that the top 30 genes shown in heatmap vs volcano plot are different


# Save plots
# one way; save dat.long to a variable ex: p<- dat.long
# use function ggsave, 1st parameter is the file variable p, and then define filename + specify dimensions
#ggsave(p, filename = 'heatmap_save1.pdf', width = 10, height = 8)

# second way to save plot is typing this before the command dat.long etc
# pdf("heatmap_save2.pdf", width = 10, height = 8)
# must write the below after the command is done as well
# dev.off()

# Save DESeq results in TSV file for BETA analysis
resLFC_beta <- data.frame(
  genename = rownames(resLFC),
  log2FC = resLFC$log2FoldChange,
  pval = resLFC$pvalue
)
resLFC_beta <- na.omit(resLFC_beta)
resLFC_beta
write.table(resLFC_beta,
            file = "L1vsRA_DESeq.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)


nrow(sigs)
nrow(resLFC)

# Count significant upregulated genes
sum(res$log2FoldChange > 0 & res$padj < 0.05, na.rm = TRUE)

# Count significant downregulated genes
sum(res$log2FoldChange < 0 & res$padj < 0.05, na.rm = TRUE)

####
# Creating chr X highlights in volcano plot and getting proportions

# Read in BETA hits
beta_hits <- read.csv("/projects/daylab/yulia/RNA-Seq/DESeq/DESeq_Results_L1vsRA/RepAd_BG4_HQfinal_uptarget_da1000_3kb.csv")
beta_genes <- beta_hits$GeneSymbol

# Read in feature counts object annotations
annotation <- as.data.frame(fc$annotation)
rownames(annotation) <- annotation$GeneID

# Create resLFC based results data frame 
res_df <- as.data.frame(resLFC)

# Map chromosome using the same IDs that became countData rownames
res_df$Chr <- annotation[rownames(res_df), "Chr"]

#Strip semicolon-separated exon duplicates ("chrX;chrX;chrX" -> "chrX") 
res_df$Chr <- sub(";.*", "", res_df$Chr)

# Flag chromosome X
x_aliases <- c("X", "chrX")
res_df$on_chrX <- !is.na(res_df$Chr) & res_df$Chr %in% x_aliases

# Flag BETA hits
res_df$is_beta <- rownames(res_df) %in% beta_genes
cat("BETA hits found in res_df:", sum(res_df$is_beta), "of", length(beta_genes), "\n")
cat("Missing BETA hits:", beta_genes[!beta_genes %in% rownames(res_df)], "\n")

# Calculate Proportions !
n_all_deg  <- sum(!is.na(res_df$pvalue))
n_sig_deg  <- sum(res_df$padj < 0.05, na.rm = TRUE)
n_chrX_all <- sum(res_df$on_chrX, na.rm = TRUE)
n_chrX_sig <- sum(res_df$on_chrX & res_df$padj < 0.05, na.rm = TRUE)

n_all_deg  
n_sig_deg 
n_chrX_all
n_chrX_sig 

# Create a secondary volcano plot with separate colors for chrX
# EnhancedVolcano base layer colors
ev_cols <- c(
  "NS"                 = "grey75",   # non-significant, no FC threshold
  "Log2 FC"            = "grey50",   # FC only
  "p-value"            = "#4393C3",  # p only  (mid-blue)
  "p-value and Log2 FC"= "#2166AC"   # significant (dark blue)
)

# Overlay colors (named for legend)
col_chrX_ns  <- "#C994C7"   # light purple  – chrX, not significant
col_chrX_sig <- "#6A3D9A"   # dark purple   – chrX, significant
col_beta     <- "#FF7F00"   # vivid orange  – BETA hits

# Build base volcano (all default labels suppressed)
volcano_chrX <- EnhancedVolcano(
  res_df,
  lab       = rownames(res_df),
  selectLab = character(0),
  x         = "log2FoldChange",
  y         = "padj",
  title     = "Volcano: L1 vs RA (with chr X DEGs highlighted)",
  subtitle  = sprintf(
    "chrX ALL DEGs: %d/%.0f (%.2f%%) | chrX SIG DEGs: %d/%d (%.2f%%)",
    n_chrX_all, n_all_deg, 100*n_chrX_all/n_all_deg,
    n_chrX_sig, n_sig_deg, 100*n_chrX_sig/n_sig_deg),
  col          = ev_cols,
  pointSize    = 1.5, labSize = 3, colAlpha = 0.6,
  max.overlaps = Inf,
  pCutoff      = 0.05, FCcutoff = 0,
  legendPosition = "none"   # we build our own legend below
)

# Add chrX and BETA layers with color mapped to legend labels
volcano_chrX <- volcano_chrX +
  # chrX non-significant
  geom_point(
    data = subset(res_df, on_chrX & (is.na(padj) | padj >= 0.05)),
    aes(x = log2FoldChange, y = -log10(padj), color = "chrX (NS)"),
    inherit.aes = FALSE, size = 1.8, alpha = 0.85
  ) +
  # chrX significant
  geom_point(
    data = subset(res_df, on_chrX & !is.na(padj) & padj < 0.05),
    aes(x = log2FoldChange, y = -log10(padj), color = "chrX (sig.)"),
    inherit.aes = FALSE, size = 2.0, alpha = 0.95
  ) +
  # BETA hits (drawn last so always on top)
  geom_point(
    data = subset(res_df, is_beta),
    aes(x = log2FoldChange, y = -log10(padj), color = "BETA hit"),
    inherit.aes = FALSE, size = 3.0, alpha = 1, shape = 18
  ) +
  # BETA labels
  geom_text_repel(
    data = subset(res_df, is_beta),
    aes(x     = log2FoldChange,
        y     = -log10(padj),
        label = rownames(subset(res_df, is_beta))),
    inherit.aes   = FALSE,
    size          = 3,
    color         = col_beta,
    fontface      = "bold",
    box.padding   = 0.45,
    point.padding = 0.3,
    max.overlaps  = Inf,
    segment.color = col_beta,
    segment.alpha = 0.6
  ) +
  # Manual color scale for the three overlay categories
  scale_color_manual(
    name   = NULL,
    values = c(
      "chrX (NS)"   = col_chrX_ns,
      "chrX (sig.)" = col_chrX_sig,
      "BETA hit"    = col_beta
    ),
    guide = guide_legend(
      override.aes = list(
        shape = c(16, 16, 18),
        size  = c(3, 3, 4),
        alpha = 1
      )
    )
  ) +
  # Append base-layer color legend entries manually
  annotate("point", x = Inf, y = Inf, color = NA) +   # spacer
  theme(
    legend.position  = "right",
    legend.text      = element_text(size = 9),
    legend.key       = element_blank(),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3)
  )


ggsave("/projects/daylab/yulia/RNA-Seq/DESeq/DESeq_Results_L1vsRA/L1vsRA_Volcano_chrX.png", volcano_chrX, width = 10, height = 8, dpi = 300)

# Non-significant genes - grey75
# Significant non-chrX genes - dark blue
# chrX non-significant - light purple 
# chrX significant - dark purple 
# BETA hits - orange diamond



