library(tximport)
library(readr)
library(tximportData)
library(DESeq2)
library(tidyverse)
library(genefilter)
library(pheatmap)

# Directory and files setup
dir <- "results/transcriptomic_data"
samples <- read.table(file.path(dir, "EBPR-Control-samples.txt"))
colnames(samples) <- c('r1', 'r2', 'sample')
rownames(samples) <- samples$sample
# split conditions between aerobic and anaerobic for calculating fold changes
samples$condition <- c("anaerobic", "anaerobic", "anaerobic", "aerobic", "aerobic", "aerobic", "aerobic")

# Load counts for each genome bin separately for DESeq calculations of fold changes to be separate
# 2767802316 IIA load
IIA_files <- file.path(dir, samples$sample, "2767802316.abundance.tsv")
names(IIA_files) <- samples$sample
IIA.txi.kallisto <- tximport(IIA_files, type="kallisto", txOut=TRUE)

# 2767802455 IIC load
IIC_files <- file.path(dir, samples$sample, "2767802455.abundance.tsv")
names(IIC_files) <- samples$sample
IIC.txi.kallisto <- tximport(IIC_files, type="kallisto", txOut=TRUE)

# Reactor4-bin.4 IIF load
IIF_files <- file.path(dir, samples$sample, "Reactor4.bin4.abundance.tsv")
names(IIF_files) <- samples$sample
IIF.txi.kallisto <- tximport(IIF_files, type="kallisto", txOut = TRUE)

# spades-bin.32 IA load
IA_files <- file.path(dir, samples$sample, "spades-bin.32.abundance.tsv")
names(IA_files) <- samples$sample
IA.txi.kallisto <- tximport(IA_files, type="kallisto", txOut = TRUE)

# DESeq Objects - by individual bin 
# For each genome, load from the tximport, filter out genes with low counts (any gene that has less than 10 counts in 3 samples)
# Generate results in log2fold


# 2767802316 IIA object
IIA.dds <- DESeqDataSetFromTximport(IIA.txi.kallisto, colData=samples, design = ~ condition)
IIA.idx <- rowSums( counts(IIA.dds) >=10 ) >=3
IIA.dds <- IIA.dds[IIA.idx, ]
IIA.dds <- DESeq(IIA.dds)
IIA.res <- results(IIA.dds)
IIA.df <- as.data.frame(IIA.res)

# 2767802455 IIC object
IIC.dds <- DESeqDataSetFromTximport(IIC.txi.kallisto, colData=samples, design = ~ condition)
IIC.idx <- rowSums( counts(IIC.dds) >=10 ) >=3
IIC.dds <- IIC.dds[IIC.idx, ]
IIC.dds <- DESeq(IIC.dds)
IIC.res <- results(IIC.dds)
IIC.df <- as.data.frame(IIC.res)

# Reactor4-bin.4 IIF object
IIF.ddsTxi <- DESeqDataSetFromTximport(IIF.txi.kallisto, colData=samples, design = ~ condition)
IIF.idx <- rowSums( counts(IIF.dds) >=10 ) >=3
IIF.dds <- IIF.dds[IIF.idx, ]
IIF.dds <- DESeq(IIF.ddsTxi)
IIF.res <- results(IIF.dds)
IIF.df <- as.data.frame(IIF.res)

# spades-bin.32 IA object
IA.ddsTxi <- DESeqDataSetFromTximport(IA.txi.kallisto, colData=samples, design = ~ condition)
IA.idx <- rowSums( counts(IA.dds) >=10 ) >=3
IA.dds <- IA.dds[IA.idx, ]
IA.dds <- DESeq(IA.ddsTxi)
IA.res <- results(IA.dds)

# --------------------------------------------

# different fold change cutoffs 
# IIA fold changes
IIA_fold <- IIA.df %>% filter(log2FoldChange > 1.5 | log2FoldChange < -1.5)
IIA_double_fold <- IIA.df %>% filter(log2FoldChange > 2 | log2FoldChange < -2)
IIA_triple_fold <- IIA.df %>% filter(log2FoldChange > 3 | log2FoldChange < -3)

# IIC fold changes
IIC_fold <- IIC.df %>% filter(log2FoldChange > 1.5 | log2FoldChange < -1.5)
IIC_double_fold <- IIC.df %>% filter(log2FoldChange > 2 | log2FoldChange < -2)
IIC_triple_fold <- IIC.df %>% filter(log2FoldChange > 3 | log2FoldChange < -3)

# IA fold changes
IA_fold <- IA.df %>% filter(log2FoldChange > 1.5 | log2FoldChange < -1.5)
IA_double_fold <- IA.df %>% filter(log2FoldChange > 2 | log2FoldChange < -2)
IA_triple_fold <- IA.df %>% filter(log2FoldChange > 3 | log2FoldChange < -3)

# IIF fold changes
IIF_fold <- IIF.df %>% filter(log2FoldChange > 1.5 | log2FoldChange < -1.5)
IIF_double_fold <- IIF.df %>% filter(log2FoldChange > 2 | log2FoldChange < -2)
IIF_triple_fold <- IIF.df %>% filter(log2FoldChange > 3 | log2FoldChange < -3)

# -------------------------

# testing heatmap functions
# get log2 fold transformed data with the `rlogTransformation` function, and turn into a matrix by assigning an assay. Can define by the top most variable genes, or a list of genes that I determine later based on COGs
rld <- rlogTransformation( IIC.dds)
topVarGenes <- head(order(-rowVars(assay(rld))), 20)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(IA.dds)[,c("r1", "r2", "sample", "condition")])
pheatmap(mat, annotation_col = df, drop_levels = TRUE)