library(tximport)
library(readr)
library(tximportData)
library(DESeq2)
library(tidyverse)
library(genefilter)
library(superheat)
library(pheatmap)
library(RColorBrewer)
library(viridis)

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
IA.df <- as.data.frame(IA.res)

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

# colors
accessory_colors <- list(group = brewer.pal(3, "Set1"))

## IIC top 20 DE genes 
IIC_metadata <- as.data.frame(colData(IIC.dds)[,c("r1", "r2", "sample", "condition")])
IIC_rld <- rlogTransformation(IIC.dds)
IIC_topVarGenes <- head(order(-rowVars(assay(IIC_rld))), 20)
IIC_top_mat <- assay(IIC_rld)[IIC_topVarGenes, ]
IIC_top_mat <- IIC_top_mat - rowMeans(IIC_top_mat)
IIC_top_heatmap <- pheatmap(IIC_top_mat, drop_levels = TRUE, cluster_cols=FALSE)

## IIC core
IIC_core <- read.csv("results/orthoTables/UW6-core-filtered.csv")
IIC_core$locus_tag <- sub("^", "2767802455_", IIC_core$locus_tag)
IIC.df <- rownames_to_column(IIC.df, var = "locus_tag") 
IIC_confirmed_core_tags <- left_join(IIC.df, IIC_core) %>% filter(!is.na(orthogroup)) %>% select(-all_or_exp_core)
IIC_core_tags <- IIC_confirmed_core_tags[['locus_tag']]
IIC_core_subset <- IIC.dds[IIC_core_tags, ]
IIC_core_rld <- rlogTransformation(IIC_core_subset)
mat <- assay(IIC_core_rld)
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(IIC.dds)[,c("r1", "r2", "sample", "condition")])
pheatmap(mat, annotation_col=df, drop_levels=TRUE, cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE)

## IIC accessory
IIC_accessory <- read.csv("results/orthoTables/UW6-accessory.csv")
IIC_accessory$locus_tag <- sub("^", "2767802455_", IIC_accessory$locus_tag)
IIC_confirmed_accessory_tags <- left_join(IIC.df, IIC_accessory) %>% filter(!is.na(orthogroup)) %>% select(-in_or_ex)
IIC_acc_tags <- IIC_confirmed_accessory_tags[['locus_tag']]
IIC_acc_subset <- IIC.dds[IIC_acc_tags, ]
IIC_acc_rld <- rlogTransformation(IIC_acc_subset)
IIC_acc_mat <- assay(IIC_acc_rld)
IIC_acc_mat <- IIC_acc_mat - rowMeans(IIC_acc_mat)

## IIC accessory plot
IIC_acc_metadata <- read.csv("results/annotations/UW6-accessory-categories.csv") %>% select(locus_tag, Category)
IIC_acc_metadata$locus_tag <- sub("^", "2767802455_", IIC_acc_metadata$locus_tag)
IIC_acc_metadata <- column_to_rownames(var="locus_tag", IIC_acc_metadata)
colors <-colorRampPalette(rev(brewer.pal(n=9,name="PuOr")))(255)
IIC_acc_mat_ordered <- IIC_acc_mat[rownames(IIC_acc_metadata), ]
IIC_accessory <- pheatmap(IIC_acc_mat, drop_levels=TRUE, show_rownames=FALSE, cluster_cols=FALSE, color=colors, annotation_row=IIC_acc_metadata, treeheight_col = 0)

## IA top 20 DE genes
IA_metadata <- as.data.frame(colData(IA.dds)[,c("r1", "r2", "sample", "condition")])
IA_rld <- rlogTransformation(IA.dds)
IA_topVarGenes <- head(order(-rowVars(assay(IA_rld))), 20)
IA_top_mat <- assay(IA_rld)[IA_topVarGenes, ]
IA_top_mat <- IA_top_mat - rowMeans(IA_top_mat)
IA_top_heatmap <- pheatmap(IA_top_mat, drop_levels = TRUE, cluster_cols=FALSE, treeheight_col=0)


## IA core
IA_core <- read.csv("results/orthoTables/UW4-core.csv")
IA_core$locus_tag <- sub("^", "spades-bin.32_", IA_core$locus_tag)
IA.df <- rownames_to_column(IA.df, var = "locus_tag") 
IA_confirmed_core_tags <- left_join(IA.df, IA_core) %>% filter(!is.na(orthogroup)) %>% select(-all_or_exp_core)
IA_core_tags <- IA_confirmed_core_tags[['locus_tag']]
IA_core_subset <- IA.dds[IA_core_tags, ]
IA_core_rld <- rlogTransformation(IA_core_subset)
IA_core_mat <- assay(IA_core_rld)
IA_core_mat <- IA_core_mat - rowMeans(IA_core_mat)
IA_df <- as.data.frame(colData(IA.dds)[,c("r1", "r2", "sample", "condition")])
pheatmap(IA_core_mat, annotation_col=IA_df, drop_levels=TRUE, cluster_rows=FALSE, cluster_cols=FALSE)

# IA accessory
IA_acc <- read.csv("results/orthoTables/UW4-accessory.csv")
IA_acc$locus_tag <- sub("^", "spades-bin.32_", IA_acc$locus_tag)
IA_confirmed_acc_tags <- left_join(IA.df, IA_acc) %>% filter(!is.na(orthogroup)) %>% select(-in_or_ex)
IA_acc_tags <- IA_confirmed_acc_tags[['locus_tag']]
IA_acc_subset <- IA.dds[IA_acc_tags, ]
IA_acc_rld <- rlogTransformation(IA_acc_subset)
IA_acc_mat <- assay(IA_acc_rld)
IA_acc_mat <- IA_acc_mat - rowMeans(IA_acc_mat)

# IA accessory plot
# IA row function categories metadata
IA_acc_categories_metadata <- read.csv("results/annotations/UW4-accessory-categories.csv") %>% select(locus_tag, Category)
IA_acc_categories_metadata$locus_tag <- sub("^", "spades-bin.32_", IA_acc_categories_metadata$locus_tag)
IA_acc_categories_metadata <- column_to_rownames(IA_acc_categories_metadata, var="locus_tag")
colors <-colorRampPalette(rev(brewer.pal(n=9,name="PuOr")))(255)
IA_accessory <- pheatmap(IA_acc_mat, drop_levels = TRUE, cluster_cols = FALSE, cluster_rows = TRUE, color=colors, annotation_row = IA_acc_categories_metadata, treeheight_col = 0,show_rownames=FALSE)

################# 
# save figures

ggsave("figures/IIC_top_20DE_heatmap.png", IIC_top_heatmap, height=20, width=20, units=c("cm"))
ggsave("figures/IA_top_20DE_heatmap.png", IA_top_heatmap, height=20, width=20, units=c("cm"))
ggsave("figures/IA_accessory_heatmap.png", IA_accessory, height=20, width=20, units=c("cm"))
ggsave("figures/IIC_accessory_heatmap.png", IIC_accessory, height=20, width=20, units=c("cm"))
