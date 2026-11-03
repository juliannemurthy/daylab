---
  title: "30-1172669254 C. elegans RNA-seq RA vs Aged
output: html_document
author: lindsay julio modified by Shannon MacLeod
date: "2025-05-21"
editor_options: 
  chunk_output_type: console
---
  
  This R Notebook denotes the steps to perform following the steps performed in terminal to trim, align and count features of your RNA-seq output. 

## References
Most of this guide is based on the instructions taught in "HarvardX: Case Studies in Functional Genomics" offered on EdX in April 2024. The bioconductor document detailing RNA-seq data analysis with DESeq2 was a useful resource for data shrinkage: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering
Azenta's guidelines were also followed: https://www.azenta.com/blog/quick-start-guide-rna-seq-data-analysis 

## Necessary packages to install

```{r}
#install.packages(c("foreign", "KernSmooth", "lattice", "Matrix", "mgcv", "nlme", "rpart", "spatial", "survival", "stringr"), lib = "/home/sykes-macleod.s/R/x86_64-pc-linux-gnu-library/4.4")

#BiocManager::install(c("BiocGenerics", "BiocVersion", "rhdf5", "rhdf5filters", "Rhdf5lib", "S4Vectors", "Rsubread", "DESeq2", "vsn", "apeglm", "ashr", "pheatmap", "org.Hs.eg.db", "EnsDb.Hsapiens.v86", "EnhancedVolcano"), force = TRUE, lib = "/home/sykes-macleod.s/R/x86_64-pc-linux-gnu-library/4.4")
```

```{r}

#BiocManager::install(c("BiocGenerics","S4Vectors", "DESeq2", "vsn", "apeglm", "ashr", "org.Hs.eg.db", "EnsDb.Hsapiens.v86", "EnhancedVolcano"), force = TRUE, lib = "/home/sykes-macleod.s/R/x86_64-pc-linux-gnu-library/4.4")

#BiocManager::install("clusterProfiler", force = TRUE, lib = "/home/sykes-macleod.s/R/x86_64-pc-linux-gnu-library/4.4")

```



## Feature Counts

You can either count the features using a bash command (see RNA-seq analysis guide document) or you can do it here in R. 

First we need to make a csv file containing the sample information. With in the first column the file name and in the second column the condition that sample is assigned to (example: file name: A549-siBAZ2A-1_S1_L001.sorted.bam;	condition: A549-siBAZ2A). 

```{r}
dir <- "/work/daylab/Shannon/RNA-seq/30-1172669254/STAR_Aligned"
setwd(dir)
## List the bam files you want to get the counts for
fileName <- list.files(dir, ".sorted.bam$")
## Extract the sample ID as well as the o=borader condition information 
SampleID <- sub('\\.sorted.bam$', '', fileName)
Condition <- substr(SampleID, 1, nchar(SampleID)-2)
## Create your CSV file 
df <- data.frame(unlist(fileName), unlist(Condition))
names(df) = c("FileName", "Condition")
write.csv(df, "sample_table.csv")
df
```

Then we need to define the location of the bam files we want to analyze as well as the gtf file, where-after we will generate a feature count list using the options recommended by Azenta. 


AT THIS POINT IF THIS CODE ISN'T RUNNING, YOU CAN USE Launch_SubRead.sh AND SubRead_Job.sh IN /work/daylab/lindsay/rna_seq/analysis!
  ```{r}
## define the necessary directories 
dir <- "/work/daylab/Shannon/RNA-seq/30-1172669254/STAR_Aligned"
dirgtf <- "/work/daylab/Shannon/RNA-seq/geneAnnotations/WBcel235"

## Read the csv file we created before and list all the bam files as well as the gtf file
csv.file <- file.path(dir, "sample_table.csv")
sample.table <- read.csv(csv.file, row.names=1)
bam.files <- file.path(dir, paste0(sample.table$FileName))
gtf.file <- file.path(dirgtf, "WBcel235genomic.gtf")

## Count the features and get a dataframe with the counts as result
#library("Rsubread")
#fc <- featureCounts(bam.files, annot.ext=gtf.file,
#                    isGTFAnnotationFile=TRUE, 
#                    countMultiMappingReads=FALSE,
#                    GTF.featureType="exon",
#                    largestOverlap=TRUE,
#                    GTF.attrType="gene_id")
```

rep1_01hr.sorted.bam rep1_24hr.sorted.bam rep1_48hr.sorted.bam rep1_72hr.sorted.bam rep1_Sham.sorted.bam rep2_01hr.sorted.bam rep2_24hr.sorted.bam rep2_48hr.sorted.bam rep2_72hr.sorted.bam rep2_Sham.sorted.bam rep3_01hr.sorted.bam rep3_24hr.sorted.bam rep3_48hr.sorted.bam rep3_72hr.sorted.bam rep3_Sham.sorted.bam rep4_01hr.sorted.bam rep4_24hr.sorted.bam rep4_48hr.sorted.bam rep4_72hr.sorted.bam rep4_Sham.sorted.bam rep5_01hr.sorted.bam rep5_24hr.sorted.bam rep5_48hr.sorted.bam rep5_72hr.sorted.bam rep5_Sham.sorted.bam rep6_01hr.sorted.bam rep6_24hr.sorted.bam rep6_48hr.sorted.bam rep6_72hr.sorted.bam rep6_Sham.sorted.bam


rep1_01hr_counts.txt
rep1_24hr_counts.txt
rep1_48hr_counts.txt
rep1_72hr_counts.txt
rep1_Sham_counts.txt
rep2_01hr_counts.txt
rep2_24hr_counts.txt
rep2_48hr_counts.txt
rep2_72hr_counts.txt
rep2_Sham_counts.txt
rep3_01hr_counts.txt
rep3_24hr_counts.txt
rep3_48hr_counts.txt
rep3_72hr_counts.txt
rep3_Sham_counts.txt
rep4_01hr_counts.txt
rep4_24hr_counts.txt
rep4_48hr_counts.txt
rep4_72hr_counts.txt
rep4_Sham_counts.txt
rep5_01hr_counts.txt
rep5_24hr_counts.txt
rep5_48hr_counts.txt
rep5_72hr_counts.txt
rep5_Sham_counts.txt
rep6_01hr_counts.txt
rep6_24hr_counts.txt
rep6_48hr_counts.txt
rep6_72hr_counts.txt
rep6_Sham_counts.txt

If you ran SubRead_Job_combined.sh in /work/daylab/lindsay/rna_seq/analysis, read the txt file into the df "fc" before moving on
```{r}
dir <- "/work/daylab/Shannon/RNA-seq/30-1172669254/analysis/counts"
setwd(dir)
fc <- read.table("testRA_A_All.counts.txt", sep = "\t", header=TRUE)

row.names(fc) = fc[,1]

fc <- fc[,c(7:12)]

head(fc)
```





## Compare hit counts between groups with DESeq2 

Generate a DESeq object to be able to normalize the sample for sequencing depth
```{r}
library("DESeq2")
dds.fc <- DESeqDataSetFromMatrix(fc, df, design= ~Condition)
## Make sure that the reference level is siControl
dds.fc$Condition <- relevel(dds.fc$Condition, "RA")
## Check if the levels are in the order of comparison you want
levels(dds.fc$Condition)

## Run the DESeq2 model
dds <- DESeq(dds.fc)
```



If you want to look at the size factors (these will give you a robust estimate of the differences in sequencing depth), you can run the following line: 
  ```{r}
## This will give you the numbers
sizeFactors(dds)
## This will give you a visual representation
plot(sizeFactors(dds), colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))
```
Samples with a size factor > 1 had a higher sequencing depth, samples with a size vactor < 1 had a lower sequencing depth. 


## [OPTIONAL] Preliminary data visualisation

We will want to transform the data in a way which minimizes differences between samples for rows with small counts and which normalizes with respect to library size. To do this in a preliminary way to be able to inspect the data we will use the rlog function of DESeq2. 

```{r}
rld <- rlog(dds)
```

[OPTIONAL] If you would like to visually see what it does to your data by comparing the 2 first samples you can run this code: 
  ```{r}

log.norm.counts <- log2(counts(dds, normalize = TRUE) + 1)
## plot of only the log2 without additional normalization
plot(log.norm.counts[,1:2], cex=.1)
## Plot of the rlog normalized data
plot(assay(rld)[,1], assay(rld)[,2], cex=.1)

## Different representation
library(vsn)
## plot of only the log2 without additional normalization
meanSdPlot(log.norm.counts, rank=FALSE)
## Plot of the rlog normalized data
meanSdPlot(assay(rld), rank=FALSE)

## To get a quick idea of the variance within your samples based on conditions, make a PCA plot:
library(ggplot2)
library(ggrepel)
z <- plotPCA(rld, intgroup="Condition")
nudge <- position_nudge(y = 1)
z + theme_bw()
z + geom_label_repel(aes(label = name), position = nudge)
```


## Compare hit counts between groups with DESeq2 (continued)

You will need to generate  results table to be able to make visual representations of the 
data as well as to analyze the differential gene expression. We will want to transform the
data in a way which minimizes differences between samples for rows with small counts and 
which normalizes with respect to library size. To do this we will use the log fold change 
shrinkage analysis using the lfcShrink function of DESeq2. There are three different types 
of shrinkage analysis you can use here: apeglm , ashr and normal. It seems like apeglm is the 
most recommended one online, however this yielded the following message for most of my data: 
"the line search routine failed, unable to sufficiently decrease the function value". 
In that case ashr seems to be the 2nd most recommended type. 

```{r}
## Log fold change shrinkage analysis 
resLFC <- lfcShrink(dds, coef = "Condition_Aged_vs_RA", type="apeglm")

## Visual representation 
plotMA(resLFC)
```


## Heatmap creation 

```{r}
library(pheatmap)
library(stringr)

### Plot based on top 20 siBAZ2B padjusted genes
resSort48hr <- resLFC.48hr[order(resLFC.48hr$padj),]
topgenes <- head(rownames(resSort48hr),25)
mat <- assay(rld)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,"Condition"])
colnames(mat) <- str_sub(colnames(mat), 1, -12)
rownames(df) <- colnames(mat)
pheatmap(mat, annotation_col=df)
```
```{r}

library(EnsDb.Hsapiens.v86)
library(pheatmap)
library(stringr)
library(viridis)

### Plot based on top 20 siBAZ2B padjusted genes
resSort48hr <- resLFC.48hr[order(resLFC.48hr$padj),]
topgenes <- head(rownames(resSort48hr),25)
mat <- assay(rld)[topgenes,]
mat <- mat - rowMeans(mat)

rownames(mat) <- sub("\\_.*", "", rownames(mat))
rownames(mat) <- sub("\\..*", "", rownames(mat))
EnsID <- row.names(mat)
geneID <- ensembldb::select(EnsDb.Hsapiens.v86, keys=EnsID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
rownames(mat)[match(geneID[,2], row.names(mat))] <- geneID[,1]


df <- as.data.frame(colData(dds)[,"Condition"])
colnames(mat) <- str_sub(colnames(mat), 1, -12)
rownames(df) <- colnames(mat)
pheatmap(mat, annotation_col=df, color = inferno(20), main = "48h DEGs")


```




### Volcano plot creation

```{r}
library(EnhancedVolcano)

table(resLFC.48hr$padj < 0.05)
EnhancedVolcano(resLFC.48hr,
                lab = rownames(resLFC.48hr),
                x = 'log2FoldChange',
                y = 'pvalue')

```



```{r}
library(EnsDb.Hsapiens.v86)
library(EnhancedVolcano)
library(stringr)

rownames(resLFC.48hr) <- sub("\\_.*", "", rownames(resLFC.48hr))
rownames(resLFC.48hr) <- sub("\\..*", "", rownames(resLFC.48hr))
EnsID <- row.names(resLFC.48hr)
geneID <- ensembldb::select(EnsDb.Hsapiens.v86, keys=EnsID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
rownames(resLFC.48hr)[match(geneID[,2], row.names(resLFC.48hr))] <- geneID[,1]

table(resLFC.48hr$padj < 0.05)
EnhancedVolcano(resLFC.48hr,
                lab = rownames(resLFC.48hr),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 0.5,
                title = 'Differentially Expressed Genes at 48h post-UV',
                labSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)

```









01hr

```{r}
## Log fold change shrinkage analysis 
resLFC.01hr <- lfcShrink(dds, coef = "Condition_01hr_vs_Sham", type="apeglm")

## Visual representation 
plotMA(resLFC.01hr) 
```

```{r}
library(pheatmap)
library(stringr)

### Plot based on top 20 padjusted genes
resSort01hr <- resLFC.01hr[order(resLFC.01hr$padj),]
topgenes <- head(rownames(resSort01hr),25)
mat <- assay(rld)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,"Condition"])
colnames(mat) <- str_sub(colnames(mat), 1, -3)
rownames(df) <- colnames(mat)
pheatmap(mat, annotation_col=df)
```

```{r}

library(EnsDb.Hsapiens.v86)
library(pheatmap)
library(stringr)

### Plot based on top 20 padjusted genes
resSort01hr <- resLFC.01hr[order(resLFC.01hr$padj),]
topgenes <- head(rownames(resSort01hr),25)
mat <- assay(rld)[topgenes,]
mat <- mat - rowMeans(mat)

rownames(mat) <- sub("\\_.*", "", rownames(mat))
rownames(mat) <- sub("\\..*", "", rownames(mat))
EnsID <- row.names(mat)
geneID <- ensembldb::select(EnsDb.Hsapiens.v86, keys=EnsID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
rownames(mat)[match(geneID[,2], row.names(mat))] <- geneID[,1]


df <- as.data.frame(colData(dds)[,"Condition"])
colnames(mat) <- str_sub(colnames(mat), 1, -12)
rownames(df) <- colnames(mat)
pheatmap(mat, annotation_col=df, color = inferno(20), main = "01h DEGs")


```





```{r}
library(EnhancedVolcano)

table(resLFC.01hr$padj < 0.05)
EnhancedVolcano(resLFC.01hr,
                lab = rownames(resLFC.01hr),
                x = 'log2FoldChange',
                y = 'pvalue')

```

```{r}
library(EnsDb.Hsapiens.v86)
library(EnhancedVolcano)
library(stringr)

rownames(resLFC.01hr) <- sub("\\_.*", "", rownames(resLFC.01hr))
rownames(resLFC.01hr) <- sub("\\..*", "", rownames(resLFC.01hr))
EnsID <- row.names(resLFC.01hr)
geneID <- ensembldb::select(EnsDb.Hsapiens.v86, keys=EnsID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
rownames(resLFC.01hr)[match(geneID[,2], row.names(resLFC.01hr))] <- geneID[,1]

table(resLFC.01hr$padj < 0.05)
EnhancedVolcano(resLFC.01hr,
                lab = rownames(resLFC.01hr),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 0.5,
                title = 'Differentially Expressed Genes at 01h post-UV',
                labSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)

```





24hr

```{r}
## Log fold change shrinkage analysis 
resLFC.24hr <- lfcShrink(dds, coef = "Condition_24hr_vs_Sham", type="apeglm")

## Visual representation 
plotMA(resLFC.24hr)
```

```{r}


require(gridExtra)
plot1 <- plotMA(resLFC.01hr)
plot2 <- plotMA(resLFC.24hr)
plot3 <- plotMA(resLFC.48hr)
#grid.arrange(plot1, plot2, plot3, ncol=3)

par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-0.4,0.4)
plotMA(resLFC.01hr, xlim=xlim, ylim=ylim, main="01hr vs Mock")
plotMA(resLFC.24hr, xlim=xlim, ylim=ylim, main="24hr vs Mock")
plotMA(resLFC.48hr, xlim=xlim, ylim=ylim, main="48hr vs Mock")

```


```{r}
library(pheatmap)
library(stringr)

### Plot based on top 20 siBAZ2B padjusted genes
resSort24hr <- resLFC.24hr[order(resLFC.24hr$padj),]
topgenes <- head(rownames(resSort24hr),25)
mat <- assay(rld)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,"Condition"])
colnames(mat) <- str_sub(colnames(mat), 1, -3)
rownames(df) <- colnames(mat)
pheatmap(mat, annotation_col=df)
```

```{r}

library(EnsDb.Hsapiens.v86)
library(pheatmap)
library(stringr)

### Plot based on top 20 padjusted genes
resSort24hr <- resLFC.24hr[order(resLFC.24hr$padj),]
topgenes <- head(rownames(resSort24hr),25)
mat <- assay(rld)[topgenes,]
mat <- mat - rowMeans(mat)

rownames(mat) <- sub("\\_.*", "", rownames(mat))
rownames(mat) <- sub("\\..*", "", rownames(mat))
EnsID <- row.names(mat)
geneID <- ensembldb::select(EnsDb.Hsapiens.v86, keys=EnsID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
rownames(mat)[match(geneID[,2], row.names(mat))] <- geneID[,1]


df <- as.data.frame(colData(dds)[,"Condition"])
colnames(mat) <- str_sub(colnames(mat), 1, -12)
rownames(df) <- colnames(mat)
pheatmap(mat, annotation_col=df, color = inferno(20), main = "24h DEGs")


```



```{r}
library(EnhancedVolcano)

table(resLFC.24hr$padj < 0.05)
EnhancedVolcano(resLFC.24hr,
                lab = rownames(resLFC.24hr),
                x = 'log2FoldChange',
                y = 'pvalue')

```

```{r}
library(EnsDb.Hsapiens.v86)
library(EnhancedVolcano)
library(stringr)

rownames(resLFC.24hr) <- sub("\\_.*", "", rownames(resLFC.24hr))
rownames(resLFC.24hr) <- sub("\\..*", "", rownames(resLFC.24hr))
EnsID <- row.names(resLFC.24hr)
geneID <- ensembldb::select(EnsDb.Hsapiens.v86, keys=EnsID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
rownames(resLFC.24hr)[match(geneID[,2], row.names(resLFC.24hr))] <- geneID[,1]

table(resLFC.24hr$padj < 0.05)
EnhancedVolcano(resLFC.24hr,
                lab = rownames(resLFC.24hr),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 0.5,
                title = 'Differentially Expressed Genes at 24h post-UV',
                labSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)

```








72hr
```{r}
## Log fold change shrinkage analysis 
resLFC.72hr <- lfcShrink(dds, coef = "Condition_01hr_vs_Sham", type="apeglm")

## Visual representation 
plotMA(resLFC.72hr)
```

```{r}
library(pheatmap)
library(stringr)

### Plot based on top 20 siBAZ2B padjusted genes
resSort72hr <- resLFC.72hr[order(resLFC.72hr$padj),]
topgenes <- head(rownames(resSort72hr),25)
mat <- assay(rld)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,"Condition"])
colnames(mat) <- str_sub(colnames(mat), 1, -3)
rownames(df) <- colnames(mat)
pheatmap(mat, annotation_col=df)
```


```{r}

library(EnsDb.Hsapiens.v86)
library(pheatmap)
library(stringr)

### Plot based on top 20 padjusted genes
resSort72hr <- resLFC.72hr[order(resLFC.72hr$padj),]
topgenes <- head(rownames(resSort72hr),25)
mat <- assay(rld)[topgenes,]
mat <- mat - rowMeans(mat)

rownames(mat) <- sub("\\_.*", "", rownames(mat))
rownames(mat) <- sub("\\..*", "", rownames(mat))
EnsID <- row.names(mat)
geneID <- ensembldb::select(EnsDb.Hsapiens.v86, keys=EnsID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
rownames(mat)[match(geneID[,2], row.names(mat))] <- geneID[,1]


df <- as.data.frame(colData(dds)[,"Condition"])
colnames(mat) <- str_sub(colnames(mat), 1, -12)
rownames(df) <- colnames(mat)
pheatmap(mat, annotation_col=df, color = inferno(20))


```




```{r}
library(EnhancedVolcano)

table(resLFC.72hr$padj < 0.05)
EnhancedVolcano(resLFC.72hr,
                lab = rownames(resLFC.72hr),
                x = 'log2FoldChange',
                y = 'pvalue')

```


```{r}
library(EnsDb.Hsapiens.v86)
library(EnhancedVolcano)
library(stringr)

rownames(resLFC.72hr) <- sub("\\_.*", "", rownames(resLFC.72hr))
rownames(resLFC.72hr) <- sub("\\..*", "", rownames(resLFC.72hr))
EnsID <- row.names(resLFC.72hr)
geneID <- ensembldb::select(EnsDb.Hsapiens.v86, keys=EnsID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
rownames(resLFC.72hr)[match(geneID[,2], row.names(resLFC.72hr))] <- geneID[,1]

table(resLFC.72hr$padj < 0.05)
EnhancedVolcano(resLFC.72hr,
                lab = rownames(resLFC.72hr),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 0.5,
                title = 'Differentially Expressed Genes at 72h post-UV',
                labSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)

```





### Export results table

```{r}
library(EnsDb.Hsapiens.v86)
clean.reSort72hr <- na.omit(resSort72hr)
rownames(clean.reSort72hr) <- sub("\\_.*", "", rownames(clean.reSort72hr))
rownames(clean.reSort72hr) <- sub("\\..*", "", rownames(clean.reSort72hr))
EnsID <- row.names(clean.reSort72hr)
geneID <- ensembldb::select(EnsDb.Hsapiens.v86, keys=EnsID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
new72hr <- data.frame(GENEID = row.names(clean.reSort72hr), clean.reSort72hr)
gene72hr <- merge(geneID, new72hr, by="GENEID")

library(EnsDb.Hsapiens.v86)
clean.reSort48hr <- na.omit(resSort48hr)
rownames(clean.reSort48hr) <- sub("\\_.*", "", rownames(clean.reSort48hr))
rownames(clean.reSort48hr) <- sub("\\..*", "", rownames(clean.reSort48hr))
EnsID <- row.names(clean.reSort48hr)
geneID <- ensembldb::select(EnsDb.Hsapiens.v86, keys=EnsID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
new48hr <- data.frame(GENEID = row.names(clean.reSort48hr), clean.reSort48hr)
gene48hr <- merge(geneID, new48hr, by="GENEID")

library(EnsDb.Hsapiens.v86)
clean.reSort24hr <- na.omit(resSort24hr)
rownames(clean.reSort24hr) <- sub("\\_.*", "", rownames(clean.reSort24hr))
rownames(clean.reSort24hr) <- sub("\\..*", "", rownames(clean.reSort24hr))
EnsID <- row.names(clean.reSort24hr)
geneID <- ensembldb::select(EnsDb.Hsapiens.v86, keys=EnsID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
new24hr <- data.frame(GENEID = row.names(clean.reSort24hr), clean.reSort24hr)
gene24hr <- merge(geneID, new24hr, by="GENEID")

library(EnsDb.Hsapiens.v86)
clean.reSort01hr <- na.omit(resSort01hr)
rownames(clean.reSort01hr) <- sub("\\_.*", "", rownames(clean.reSort01hr))
rownames(clean.reSort01hr) <- sub("\\..*", "", rownames(clean.reSort01hr))
EnsID <- row.names(clean.reSort01hr)
geneID <- ensembldb::select(EnsDb.Hsapiens.v86, keys=EnsID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
new01hr <- data.frame(GENEID = row.names(clean.reSort01hr), clean.reSort01hr)
gene01hr <- merge(geneID, new01hr, by="GENEID")


```


```{r}
write.csv(gene72hr, file="A549_UV_timecourse_RNA_seq_72hr.csv")
write.csv(gene48hr, file="A549_UV_timecourse_RNA_seq_48hr.csv")
write.csv(gene24hr, file="A549_UV_timecourse_RNA_seq_24hr.csv")
write.csv(gene01hr, file="A549_UV_timecourse_RNA_seq_01hr.csv")

```






RERUNNING THINGS WITHOUT REP4_48HR REP2_SHAM and REP6_1HR

```{r eval=FALSE, include=FALSE}
# after seeing the clustering of "fc" i am trying again without reps 2,4,and 6, as well as without just rep2_sham, rep4_48hr, and rep6_01hr

fc_no_outliers <- fc[,c(1:9, 11:17, 19:25, 27:30)]
df_no_outliers <- df[c(1:9, 11:17, 19:25, 27:30),]

```

DESeq object for no outliers
```{r eval=FALSE, include=FALSE}
library("DESeq2")
dds.fc_no_outliers <- DESeqDataSetFromMatrix(fc_no_outliers, df_no_outliers, design= ~Condition)
## Make sure that the reference level is siControl
dds.fc_no_outliers$Condition <- relevel(dds.fc_no_outliers$Condition, "Sham")
## Check if the levels are in the order of comparison you want
levels(dds.fc_no_outliers$Condition)

## Run the DESeq2 model
dds_no_outliers <- DESeq(dds.fc_no_outliers)
```


```{r eval=FALSE, include=FALSE}
## This will give you the numbers
sizeFactors(dds_no_outliers)
## This will give you a visual representation
plot(sizeFactors(dds_no_outliers), colSums(counts(dds_no_outliers)))
abline(lm(colSums(counts(dds_no_outliers)) ~ sizeFactors(dds_no_outliers) + 0))
```



```{r eval=FALSE, include=FALSE}
rld_no_outliers <- rlog(dds_no_outliers)

log.norm.counts_no_outliers <- log2(counts(dds_no_outliers, normalize = TRUE) + 1)
## plot of only the log2 without additional normalization
plot(log.norm.counts_no_outliers[,1:2], cex=.1)
## Plot of the rlog normalized data
plot(assay(rld_no_outliers)[,1], assay(rld_no_outliers)[,2], cex=.1)

## Different representation
library(vsn)
## plot of only the log2 without additional normalization
meanSdPlot(log.norm.counts_no_outliers, rank=FALSE)
## Plot of the rlog normalized data
meanSdPlot(assay(rld_no_outliers), rank=FALSE)

## To get a quick idea of the variance within your samples based on conditions, make a PCA plot:
library(ggplot2)
library(ggrepel)
z_no_outliers <- plotPCA(rld_no_outliers, intgroup="Condition")
nudge <- position_nudge(y = 1)
z_no_outliers
z_no_outliers + geom_label_repel(aes(label = name), position = nudge)

```
## FIND DEGS when excluding the outliers

```{r eval=FALSE, include=FALSE}
## Log fold change shrinkage analysis 
resLFC.48hr_no_outliers <- lfcShrink(dds_no_outliers, coef = "Condition_48hr_vs_Sham", type="apeglm")

## Visual representation 
plotMA(resLFC.48hr_no_outliers)
```


## Heatmap creation 

```{r eval=FALSE, include=FALSE}
library(pheatmap)
library(stringr)

### Plot based on top 20 siBAZ2B padjusted genes
resSort48hr_no_outliers <- resLFC.48hr_no_outliers[order(resLFC.48hr_no_outliers$padj),]
topgenes_no_outliers <- head(rownames(resSort48hr_no_outliers),25)
mat_no_outliers <- assay(rld_no_outliers)[topgenes_no_outliers,]
mat_no_outliers <- mat_no_outliers - rowMeans(mat_no_outliers)
df_no_outliers <- as.data.frame(colData(dds_no_outliers)[,"Condition"])
colnames(mat_no_outliers) <- str_sub(colnames(mat_no_outliers), 1, -3)
rownames(df_no_outliers) <- colnames(mat_no_outliers)
pheatmap(mat_no_outliers, annotation_col=df_no_outliers)
```

### Volcano plot creation

```{r eval=FALSE, include=FALSE}
library(EnhancedVolcano)

table(resLFC.48hr_no_outliers$padj < 0.05)
EnhancedVolcano(resLFC.48hr_no_outliers,
                lab = rownames(resLFC.48hr_no_outliers),
                x = 'log2FoldChange',
                y = 'pvalue')

```


01hr

```{r eval=FALSE, include=FALSE}
## Log fold change shrinkage analysis 
resLFC.01hr_no_outliers <- lfcShrink(dds_no_outliers, coef = "Condition_01hr_vs_Sham", type="apeglm")

## Visual representation 
plotMA(resLFC.01hr_no_outliers)
```

```{r eval=FALSE, include=FALSE}
library(pheatmap)
library(stringr)

### Plot based on top 20 siBAZ2B padjusted genes
resSort01hr_no_outliers <- resLFC.01hr_no_outliers[order(resLFC.01hr_no_outliers$padj),]
topgenes_no_outliers <- head(rownames(resSort01hr_no_outliers),25)
mat_no_outliers <- assay(rld_no_outliers)[topgenes_no_outliers,]
mat_no_outliers <- mat_no_outliers - rowMeans(mat_no_outliers)
df_no_outliers <- as.data.frame(colData(dds_no_outliers)[,"Condition"])
colnames(mat_no_outliers) <- str_sub(colnames(mat_no_outliers), 1, -3)
rownames(df_no_outliers) <- colnames(mat_no_outliers)
pheatmap(mat_no_outliers, annotation_col=df_no_outliers)
```

```{r eval=FALSE, include=FALSE}
library(EnhancedVolcano)

table(resLFC.01hr_no_outliers$padj < 0.05)
EnhancedVolcano(resLFC.01hr_no_outliers,
                lab = rownames(resLFC.01hr_no_outliers),
                x = 'log2FoldChange',
                y = 'pvalue')

```


24hr

```{r eval=FALSE, include=FALSE}
## Log fold change shrinkage analysis 
resLFC.24hr_no_outliers <- lfcShrink(dds_no_outliers, coef = "Condition_24hr_vs_Sham", type="apeglm")

## Visual representation 
plotMA(resLFC.24hr_no_outliers)
```

```{r eval=FALSE, include=FALSE}
library(pheatmap)
library(stringr)

### Plot based on top 20 siBAZ2B padjusted genes
resSort24hr_no_outliers <- resLFC.24hr_no_outliers[order(resLFC.24hr_no_outliers$padj),]
topgenes_no_outliers <- head(rownames(resSort24hr_no_outliers),25)
mat_no_outliers <- assay(rld_no_outliers)[topgenes_no_outliers,]
mat_no_outliers <- mat_no_outliers - rowMeans(mat_no_outliers)
df_no_outliers <- as.data.frame(colData(dds_no_outliers)[,"Condition"])
colnames(mat_no_outliers) <- str_sub(colnames(mat_no_outliers), 1, -3)
rownames(df_no_outliers) <- colnames(mat_no_outliers)
pheatmap(mat_no_outliers, annotation_col=df_no_outliers)
```

```{r eval=FALSE, include=FALSE}
library(EnhancedVolcano)

table(resLFC.24hr_no_outliers$padj < 0.05)
EnhancedVolcano(resLFC.24hr_no_outliers,
                lab = rownames(resLFC.24hr_no_outliers),
                x = 'log2FoldChange',
                y = 'pvalue')

```


72hr
```{r eval=FALSE, include=FALSE}
## Log fold change shrinkage analysis 
resLFC.72hr_no_outliers <- lfcShrink(dds_no_outliers, coef = "Condition_72hr_vs_Sham", type="apeglm")

## Visual representation 
plotMA(resLFC.72hr_no_outliers)
```

```{r eval=FALSE, include=FALSE}
library(pheatmap)
library(stringr)

### Plot based on top 20 siBAZ2B padjusted genes
resSort72hr_no_outliers <- resLFC.72hr_no_outliers[order(resLFC.72hr_no_outliers$padj),]
topgenes_no_outliers <- head(rownames(resSort72hr_no_outliers),25)
mat_no_outliers <- assay(rld_no_outliers)[topgenes_no_outliers,]
mat_no_outliers <- mat_no_outliers - rowMeans(mat_no_outliers)
df_no_outliers <- as.data.frame(colData(dds_no_outliers)[,"Condition"])
colnames(mat_no_outliers) <- str_sub(colnames(mat_no_outliers), 1, -3)
rownames(df_no_outliers) <- colnames(mat_no_outliers)
pheatmap(mat_no_outliers, annotation_col=df_no_outliers)
```

```{r eval=FALSE, include=FALSE}
library(EnhancedVolcano)

table(resLFC.72hr_no_outliers$padj < 0.05)
EnhancedVolcano(resLFC.72hr_no_outliers,
                lab = rownames(resLFC.72hr_no_outliers),
                x = 'log2FoldChange',
                y = 'pvalue')

```

### Export results table

```{r eval=FALSE, include=FALSE}
library(EnsDb.Hsapiens.v86)
clean.reSort72hr_no_outliers <- na.omit(resSort72hr_no_outliers)
rownames(clean.reSort72hr_no_outliers) <- sub("\\_.*", "", rownames(clean.reSort72hr_no_outliers))
rownames(clean.reSort72hr_no_outliers) <- sub("\\..*", "", rownames(clean.reSort72hr_no_outliers))
EnsID_no_outliers <- row.names(clean.reSort72hr_no_outliers)
geneID_no_outliers <- ensembldb::select(EnsDb.Hsapiens.v86, keys=EnsID_no_outliers, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
new72hr_no_outliers <- data.frame(GENEID = row.names(clean.reSort72hr_no_outliers), clean.reSort72hr_no_outliers)
gene72hr_no_outliers <- merge(geneID_no_outliers, new72hr_no_outliers, by="GENEID")

library(EnsDb.Hsapiens.v86)
clean.reSort48hr_no_outliers <- na.omit(resSort48hr_no_outliers)
rownames(clean.reSort48hr_no_outliers) <- sub("\\_.*", "", rownames(clean.reSort48hr_no_outliers))
rownames(clean.reSort48hr_no_outliers) <- sub("\\..*", "", rownames(clean.reSort48hr_no_outliers))
EnsID_no_outliers <- row.names(clean.reSort48hr_no_outliers)
geneID_no_outliers <- ensembldb::select(EnsDb.Hsapiens.v86, keys=EnsID_no_outliers, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
new48hr_no_outliers <- data.frame(GENEID = row.names(clean.reSort48hr_no_outliers), clean.reSort48hr_no_outliers)
gene48hr_no_outliers <- merge(geneID_no_outliers, new48hr_no_outliers, by="GENEID")

library(EnsDb.Hsapiens.v86)
clean.reSort24hr_no_outliers <- na.omit(resSort24hr_no_outliers)
rownames(clean.reSort24hr_no_outliers) <- sub("\\_.*", "", rownames(clean.reSort24hr_no_outliers))
rownames(clean.reSort24hr_no_outliers) <- sub("\\..*", "", rownames(clean.reSort24hr_no_outliers))
EnsID_no_outliers <- row.names(clean.reSort24hr_no_outliers)
geneID_no_outliers <- ensembldb::select(EnsDb.Hsapiens.v86, keys=EnsID_no_outliers, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
new24hr_no_outliers <- data.frame(GENEID = row.names(clean.reSort24hr_no_outliers), clean.reSort24hr_no_outliers)
gene24hr_no_outliers <- merge(geneID_no_outliers, new24hr_no_outliers, by="GENEID")

library(EnsDb.Hsapiens.v86)
clean.reSort01hr_no_outliers <- na.omit(resSort01hr_no_outliers)
rownames(clean.reSort01hr_no_outliers) <- sub("\\_.*", "", rownames(clean.reSort01hr_no_outliers))
rownames(clean.reSort01hr_no_outliers) <- sub("\\..*", "", rownames(clean.reSort01hr_no_outliers))
EnsID_no_outliers <- row.names(clean.reSort01hr_no_outliers)
geneID_no_outliers <- ensembldb::select(EnsDb.Hsapiens.v86, keys=EnsID_no_outliers, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
new01hr_no_outliers <- data.frame(GENEID = row.names(clean.reSort01hr_no_outliers), clean.reSort01hr_no_outliers)
gene01hr_no_outliers <- merge(geneID_no_outliers, new01hr_no_outliers, by="GENEID")


```


```{r eval=FALSE, include=FALSE}
write.csv(gene72hr_no_outliers, file="A549_UV_timecourse_RNA_seq_72hr_no_outliers.csv")
write.csv(gene48hr_no_outliers, file="A549_UV_timecourse_RNA_seq_48hr_no_outliers.csv")
write.csv(gene24hr_no_outliers, file="A549_UV_timecourse_RNA_seq_24hr_no_outliers.csv")
write.csv(gene01hr_no_outliers, file="A549_UV_timecourse_RNA_seq_01hr_no_outliers.csv")

```




















TRYING THINGS WITHOUT ALL OF REPS 2 4 and 6!!!!
  
  ```{r eval=FALSE, include=FALSE}
# after seeing the clustering of "fc" i am trying again without reps 2,4,and 6, as well as without just rep2_sham, rep4_48hr, and rep6_01hr

fc_reps_1_3_5  <- fc[,c(1:5, 11:15, 21:25)]
df_reps_1_3_5 <- df[c(1:5, 11:15, 21:25),]
```


DESeq object for only reps 1, 3, and 5
```{r eval=FALSE, include=FALSE}
library("DESeq2")
dds.fc_reps_1_3_5 <- DESeqDataSetFromMatrix(fc_reps_1_3_5, df_reps_1_3_5, design= ~Condition)
## Make sure that the reference level is siControl
dds.fc_reps_1_3_5$Condition <- relevel(dds.fc_reps_1_3_5$Condition, "Sham")
## Check if the levels are in the order of comparison you want
levels(dds.fc_reps_1_3_5$Condition)

## Run the DESeq2 model
dds_reps_1_3_5 <- DESeq(dds.fc_reps_1_3_5)
```



```{r eval=FALSE, include=FALSE}
## This will give you the numbers
sizeFactors(dds_reps_1_3_5)
## This will give you a visual representation
plot(sizeFactors(dds_reps_1_3_5), colSums(counts(dds_reps_1_3_5)))
abline(lm(colSums(counts(dds_reps_1_3_5)) ~ sizeFactors(dds_reps_1_3_5) + 0))
```

```{r eval=FALSE, include=FALSE}
rld_reps_1_3_5 <- rlog(dds_reps_1_3_5)

log.norm.counts_reps_1_3_5 <- log2(counts(dds_reps_1_3_5, normalize = TRUE) + 1)
## plot of only the log2 without additional normalization
plot(log.norm.counts_reps_1_3_5[,1:2], cex=.1)
## Plot of the rlog normalized data
plot(assay(rld_reps_1_3_5)[,1], assay(rld_reps_1_3_5)[,2], cex=.1)

## Different representation
library(vsn)
## plot of only the log2 without additional normalization
meanSdPlot(log.norm.counts_reps_1_3_5, rank=FALSE)
## Plot of the rlog normalized data
meanSdPlot(assay(rld_reps_1_3_5), rank=FALSE)

## To get a quick idea of the variance within your samples based on conditions, make a PCA plot:
library(ggplot2)
library(ggrepel)
z_reps_1_3_5 <- plotPCA(rld_reps_1_3_5, intgroup="Condition")
nudge <- position_nudge(y = 1)
z_reps_1_3_5
z_reps_1_3_5 + geom_label_repel(aes(label = name), position = nudge)

```



## FIND DEGS without reps 2, 4, and 6

```{r eval=FALSE, include=FALSE}
## Log fold change shrinkage analysis 
resLFC.48hr_reps_1_3_5 <- lfcShrink(dds_reps_1_3_5, coef = "Condition_48hr_vs_Sham", type="apeglm")

## Visual representation 
plotMA(resLFC.48hr_reps_1_3_5)
```


## Heatmap creation 

```{r eval=FALSE, include=FALSE}
library(pheatmap)
library(stringr)

### Plot based on top 20 siBAZ2B padjusted genes
resSort48hr_reps_1_3_5 <- resLFC.48hr_reps_1_3_5[order(resLFC.48hr_reps_1_3_5$padj),]
topgenes_reps_1_3_5 <- head(rownames(resSort48hr_reps_1_3_5),25)
mat_reps_1_3_5 <- assay(rld_reps_1_3_5)[topgenes_reps_1_3_5,]
mat_reps_1_3_5 <- mat_reps_1_3_5 - rowMeans(mat_reps_1_3_5)
df_reps_1_3_5 <- as.data.frame(colData(dds_reps_1_3_5)[,"Condition"])
colnames(mat_reps_1_3_5) <- str_sub(colnames(mat_reps_1_3_5), 1, -3)
rownames(df_reps_1_3_5) <- colnames(mat_reps_1_3_5)
pheatmap(mat_reps_1_3_5, annotation_col=df_reps_1_3_5)
```

### Volcano plot creation

```{r eval=FALSE, include=FALSE}
library(EnhancedVolcano)

table(resLFC.48hr_reps_1_3_5$padj < 0.05)
EnhancedVolcano(resLFC.48hr_reps_1_3_5,
                lab = rownames(resLFC.48hr_reps_1_3_5),
                x = 'log2FoldChange',
                y = 'pvalue')

```


01hr

```{r eval=FALSE, include=FALSE}
## Log fold change shrinkage analysis 
resLFC.01hr_reps_1_3_5 <- lfcShrink(dds_reps_1_3_5, coef = "Condition_01hr_vs_Sham", type="apeglm")

## Visual representation 
plotMA(resLFC.01hr_reps_1_3_5)
```

```{r eval=FALSE, include=FALSE}
library(pheatmap)
library(stringr)

### Plot based on top 20 siBAZ2B padjusted genes
resSort01hr_reps_1_3_5 <- resLFC.01hr_reps_1_3_5[order(resLFC.01hr_reps_1_3_5$padj),]
topgenes_reps_1_3_5 <- head(rownames(resSort01hr_reps_1_3_5),25)
mat_reps_1_3_5 <- assay(rld_reps_1_3_5)[topgenes_reps_1_3_5,]
mat_reps_1_3_5 <- mat_reps_1_3_5 - rowMeans(mat_reps_1_3_5)
df_reps_1_3_5 <- as.data.frame(colData(dds_reps_1_3_5)[,"Condition"])
colnames(mat_reps_1_3_5) <- str_sub(colnames(mat_reps_1_3_5), 1, -3)
rownames(df_reps_1_3_5) <- colnames(mat_reps_1_3_5)
pheatmap(mat_reps_1_3_5, annotation_col=df_reps_1_3_5)
```

```{r eval=FALSE, include=FALSE}
library(EnhancedVolcano)

table(resLFC.01hr_reps_1_3_5$padj < 0.05)
EnhancedVolcano(resLFC.01hr_reps_1_3_5,
                lab = rownames(resLFC.01hr_reps_1_3_5),
                x = 'log2FoldChange',
                y = 'pvalue')

```


24hr

```{r eval=FALSE, include=FALSE}
## Log fold change shrinkage analysis 
resLFC.24hr_reps_1_3_5 <- lfcShrink(dds_reps_1_3_5, coef = "Condition_24hr_vs_Sham", type="apeglm")

## Visual representation 
plotMA(resLFC.24hr_reps_1_3_5)
```

```{r eval=FALSE, include=FALSE}
library(pheatmap)
library(stringr)

### Plot based on top 20 siBAZ2B padjusted genes
resSort24hr_reps_1_3_5 <- resLFC.24hr_reps_1_3_5[order(resLFC.24hr_reps_1_3_5$padj),]
topgenes_reps_1_3_5 <- head(rownames(resSort24hr_reps_1_3_5),25)
mat_reps_1_3_5 <- assay(rld_reps_1_3_5)[topgenes_reps_1_3_5,]
mat_reps_1_3_5 <- mat_reps_1_3_5 - rowMeans(mat_reps_1_3_5)
df_reps_1_3_5 <- as.data.frame(colData(dds_reps_1_3_5)[,"Condition"])
colnames(mat_reps_1_3_5) <- str_sub(colnames(mat_reps_1_3_5), 1, -3)
rownames(df_reps_1_3_5) <- colnames(mat_reps_1_3_5)
pheatmap(mat_reps_1_3_5, annotation_col=df_reps_1_3_5)
```

```{r eval=FALSE, include=FALSE}
library(EnhancedVolcano)

table(resLFC.24hr_reps_1_3_5$padj < 0.05)
EnhancedVolcano(resLFC.24hr_reps_1_3_5,
                lab = rownames(resLFC.24hr_reps_1_3_5),
                x = 'log2FoldChange',
                y = 'pvalue')

```


72hr
```{r eval=FALSE, include=FALSE}
## Log fold change shrinkage analysis 
resLFC.72hr_reps_1_3_5 <- lfcShrink(dds_reps_1_3_5, coef = "Condition_01hr_vs_Sham", type="apeglm")

## Visual representation 
plotMA(resLFC.72hr_reps_1_3_5)
```

```{r eval=FALSE, include=FALSE}
library(pheatmap)
library(stringr)

### Plot based on top 20 siBAZ2B padjusted genes
resSort72hr_reps_1_3_5 <- resLFC.72hr_reps_1_3_5[order(resLFC.72hr_reps_1_3_5$padj),]
topgenes_reps_1_3_5 <- head(rownames(resSort72hr_reps_1_3_5),25)
mat_reps_1_3_5 <- assay(rld_reps_1_3_5)[topgenes_reps_1_3_5,]
mat_reps_1_3_5 <- mat_reps_1_3_5 - rowMeans(mat_reps_1_3_5)
df_reps_1_3_5 <- as.data.frame(colData(dds_reps_1_3_5)[,"Condition"])
colnames(mat_reps_1_3_5) <- str_sub(colnames(mat_reps_1_3_5), 1, -3)
rownames(df_reps_1_3_5) <- colnames(mat_reps_1_3_5)
pheatmap(mat_reps_1_3_5, annotation_col=df_reps_1_3_5)
```

```{r eval=FALSE, include=FALSE}
library(EnhancedVolcano)

table(resLFC.72hr_reps_1_3_5$padj < 0.05)
EnhancedVolcano(resLFC.72hr_reps_1_3_5,
                lab = rownames(resLFC.72hr_reps_1_3_5),
                x = 'log2FoldChange',
                y = 'pvalue')

```

### Export results table

```{r eval=FALSE, include=FALSE}
library(EnsDb.Hsapiens.v86)
clean.reSort72hr_reps_1_3_5 <- na.omit(resSort72hr_reps_1_3_5)
rownames(clean.reSort72hr_reps_1_3_5) <- sub("\\_.*", "", rownames(clean.reSort72hr_reps_1_3_5))
rownames(clean.reSort72hr_reps_1_3_5) <- sub("\\..*", "", rownames(clean.reSort72hr_reps_1_3_5))
EnsID_reps_1_3_5 <- row.names(clean.reSort72hr_reps_1_3_5)
geneID_reps_1_3_5 <- ensembldb::select(EnsDb.Hsapiens.v86, keys=EnsID_reps_1_3_5, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
new72hr_reps_1_3_5 <- data.frame(GENEID = row.names(clean.reSort72hr_reps_1_3_5), clean.reSort72hr_reps_1_3_5)
gene72hr_reps_1_3_5 <- merge(geneID_reps_1_3_5, new72hr_reps_1_3_5, by="GENEID")

library(EnsDb.Hsapiens.v86)
clean.reSort48hr_reps_1_3_5 <- na.omit(resSort48hr_reps_1_3_5)
rownames(clean.reSort48hr_reps_1_3_5) <- sub("\\_.*", "", rownames(clean.reSort48hr_reps_1_3_5))
rownames(clean.reSort48hr_reps_1_3_5) <- sub("\\..*", "", rownames(clean.reSort48hr_reps_1_3_5))
EnsID_reps_1_3_5 <- row.names(clean.reSort48hr_reps_1_3_5)
geneID_reps_1_3_5 <- ensembldb::select(EnsDb.Hsapiens.v86, keys=EnsID_reps_1_3_5, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
new48hr_reps_1_3_5 <- data.frame(GENEID = row.names(clean.reSort48hr_reps_1_3_5), clean.reSort48hr_reps_1_3_5)
gene48hr_reps_1_3_5 <- merge(geneID_reps_1_3_5, new48hr_reps_1_3_5, by="GENEID")

library(EnsDb.Hsapiens.v86)
clean.reSort24hr_reps_1_3_5 <- na.omit(resSort24hr_reps_1_3_5)
rownames(clean.reSort24hr_reps_1_3_5) <- sub("\\_.*", "", rownames(clean.reSort24hr_reps_1_3_5))
rownames(clean.reSort24hr_reps_1_3_5) <- sub("\\..*", "", rownames(clean.reSort24hr_reps_1_3_5))
EnsID_reps_1_3_5 <- row.names(clean.reSort24hr_reps_1_3_5)
geneID_reps_1_3_5 <- ensembldb::select(EnsDb.Hsapiens.v86, keys=EnsID_reps_1_3_5, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
new24hr_reps_1_3_5 <- data.frame(GENEID = row.names(clean.reSort24hr_reps_1_3_5), clean.reSort24hr_reps_1_3_5)
gene24hr_reps_1_3_5 <- merge(geneID_reps_1_3_5, new24hr_reps_1_3_5, by="GENEID")

library(EnsDb.Hsapiens.v86)
clean.reSort01hr_reps_1_3_5 <- na.omit(resSort01hr_reps_1_3_5)
rownames(clean.reSort01hr_reps_1_3_5) <- sub("\\_.*", "", rownames(clean.reSort01hr_reps_1_3_5))
rownames(clean.reSort01hr_reps_1_3_5) <- sub("\\..*", "", rownames(clean.reSort01hr_reps_1_3_5))
EnsID_reps_1_3_5 <- row.names(clean.reSort01hr_reps_1_3_5)
geneID_reps_1_3_5 <- ensembldb::select(EnsDb.Hsapiens.v86, keys=EnsID_reps_1_3_5, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
new01hr_reps_1_3_5 <- data.frame(GENEID = row.names(clean.reSort01hr_reps_1_3_5), clean.reSort01hr_reps_1_3_5)
gene01hr_reps_1_3_5 <- merge(geneID_reps_1_3_5, new01hr_reps_1_3_5, by="GENEID")


```


```{r eval=FALSE, include=FALSE}
write.csv(gene72hr_reps_1_3_5, file="A549_UV_timecourse_RNA_seq_72hr_reps_1_3_5.csv")
write.csv(gene48hr_reps_1_3_5, file="A549_UV_timecourse_RNA_seq_48hr_reps_1_3_5.csv")
write.csv(gene24hr_reps_1_3_5, file="A549_UV_timecourse_RNA_seq_24hr_reps_1_3_5.csv")
write.csv(gene01hr_reps_1_3_5, file="A549_UV_timecourse_RNA_seq_01hr_reps_1_3_5.csv")

```


GENE SET ENRICHMENT ANALYSIS

```{r}
library("clusterProfiler")
library("org.Hs.eg.db")

#geneID <- read.delim("/work/daylab/lindsay/rna_seq/20240819_A549/A549_UV_timecourse_RNA_seq_01hr.csv", sep=",")
#geneID_01h <- geneID["GENEID"]
#geneID_01h <- geneID[["GENEID"]]


#GO_01hr <- enrichGO(geneID_01h, org.Hs.eg.db, keyType = "ENSEMBL")

# reading in data from deseq2
df = read.csv("/work/daylab/lindsay/rna_seq/20240819_A549/A549_UV_timecourse_RNA_seq_01hr.csv", sep=",", header = TRUE)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$GENEID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "bonferroni")

require(DOSE)
dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign) + scale_fill_viridis(option = "plasma") + labs(title = "Enriched pathways at 01hr post UV")
#emapplot(gse, showCategory = 10)

cnetplot(gse, categorySize="pvalue", color.params = list(foldChange = gene_list), showCategory = 3) +
  scale_color_viridis(option = "plasma") +
  guides(color=guide_colorbar(title="fold change")) +
  labs(title = "Enriched pathways at 01hr post UV")

library("ggplot2")
library("ggridges")
library("viridis")
ridgeplot(gse) + labs(title = "Enriched pathways at 01hr post UV", x = "enrichment distribution") + scale_fill_viridis(option = "plasma")

```



```{r}
library("clusterProfiler")
library("org.Hs.eg.db")

# reading in data from deseq2
df = read.csv("/work/daylab/lindsay/rna_seq/20240819_A549/A549_UV_timecourse_RNA_seq_24hr.csv", sep=",", header = TRUE)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$GENEID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "bonferroni")

require(DOSE)
dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign) + scale_fill_viridis(option = "plasma") + labs(title = "Enriched pathways at 24hr post UV")
#emapplot(gse, showCategory = 10)

cnetplot(gse, categorySize="pvalue", color.params = list(foldChange = gene_list), showCategory = 3) +
  scale_color_viridis(option = "plasma") +
  guides(color=guide_colorbar(title="fold change")) +
  labs(title = "Enriched pathways at 24hr post UV")

library("ggplot2")
library("ggridges")
library("viridis")
ridgeplot(gse) + labs(title = "Enriched pathways at 24hr post UV", x = "enrichment distribution") + scale_fill_viridis(option = "plasma")



```

```{r}
library("clusterProfiler")
library("org.Hs.eg.db")

# reading in data from deseq2
df = read.csv("/work/daylab/lindsay/rna_seq/20240819_A549/A549_UV_timecourse_RNA_seq_48hr.csv", sep=",", header = TRUE)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$GENEID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "bonferroni")

require(DOSE)
dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign) + scale_fill_viridis(option = "plasma") + labs(title = "Enriched pathways at 48hr post UV")
#emapplot(gse, showCategory = 10)

cnetplot(gse, categorySize="pvalue", color.params = list(foldChange = gene_list), showCategory = 3) +
  scale_color_viridis(option = "plasma") +
  guides(color=guide_colorbar(title="fold change")) +
  labs(title = "Enriched pathways at 48hr post UV")

library("ggplot2")
library("ggridges")
library("viridis")
ridgeplot(gse) + labs(title = "Enriched pathways at 48hr post UV", x = "enrichment distribution") + scale_fill_viridis(option = "plasma")
```


```{r}
library("clusterProfiler")
library("org.Hs.eg.db")

# reading in data from deseq2
df = read.csv("/work/daylab/lindsay/rna_seq/20240819_A549/A549_UV_timecourse_RNA_seq_72hr.csv", sep=",", header = TRUE)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$GENEID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "bonferroni")

require(DOSE)
dotplot(gse, showCategory=20, split=".sign", font.size = 8) + facet_grid(.~.sign) + scale_fill_viridis(option = "plasma") + labs(title = "Enriched pathways at 72hr post UV")
#emapplot(gse, showCategory = 10)

cnetplot(gse, categorySize="pvalue", color.params = list(foldChange = gene_list), showCategory = 3) +
  scale_color_viridis(option = "plasma") +
  guides(color=guide_colorbar(title="fold change")) +
  labs(title = "Enriched pathways at 72hr post UV")

library("ggplot2")
library("ggridges")
library("viridis")
ridgeplot(gse) + labs(title = "Enriched pathways at 72hr post UV", x = "enrichment distribution", font.size = 8) + scale_fill_viridis(option = "plasma")
```


The Gene Ontology is a controlled vocabulary of terms to represent biology in a structured way. The terms are subdivided in three distinct ontologies that represent different biological aspects: Molecular Function (MF), Biological Process (BP), and Cellular Component (CC) (1). - THESE ARE THE OPTIONS FOR ont



```{r}
library("clusterProfiler")
library("org.Hs.eg.db")


# reading in data from deseq2
df = read.csv("/work/daylab/lindsay/rna_seq/20240819_A549/A549_UV_timecourse_RNA_seq_72hr.csv", sep=",", header = TRUE)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$GENEID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# Exctract significant results (padj < 0.05)
sig_genes_df = subset(df, padj < 0.05)

# From significant results, we want to filter on log2fold change
genes <- sig_genes_df$log2FoldChange

# Name the vector
names(genes) <- sig_genes_df$GENEID

# omit NA values
genes <- na.omit(genes)

# filter on min log2fold change (log2FoldChange > 2)
genes <- names(genes)[abs(genes) > 1]


go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "MF",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)



barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 20, 
        title = "GO Molecular Functions",
        font.size = 8)

dotplot(go_enrich,
        showCategory = 20, 
        title = "GO Molecular Functions",
        font.size = 8)


```

```{r}
library("clusterProfiler")
library("org.Hs.eg.db")


# reading in data from deseq2
df = read.csv("/work/daylab/lindsay/rna_seq/20240819_A549/A549_UV_timecourse_RNA_seq_01hr.csv", sep=",", header = TRUE)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$GENEID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# Exctract significant results (padj < 0.05)
sig_genes_df = subset(df, padj < 0.05)

# From significant results, we want to filter on log2fold change
genes <- sig_genes_df$log2FoldChange

# Name the vector
names(genes) <- sig_genes_df$GENEID

# omit NA values
genes <- na.omit(genes)

# filter on min log2fold change (log2FoldChange > 2)
genes <- names(genes)[abs(genes) > 1]


go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)



barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)

dotplot(go_enrich,
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)


```







