library(tximport)
library(tximportData)
library(readr)
library(tibble)
library(dplyr)

# Control experiments
dir <- "/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/R3R4/results/transcriptomic_data"
samples <- read.table(file.path(dir, 'EBPR-Control-samples.txt'), header=FALSE)
colnames(samples) <- c('r1', 'r2', 'sample')
files <- file.path(dir, samples$sample, "abundance.h5")
names(files) <- paste0("sample", 1:7)
txi.kallisto <- tximport(files, type="kallisto", txOut=TRUE)
counts <- as.data.frame(txi.kallisto)
finalcounts <- rownames_to_column(counts, var='ID')
controlRawCounts <- finalcounts[, c(1:8)]
colnames(controlRawCounts) <- c('Locus_tag', 'EBPR-Control-1045', 'EBPR-Control-1116', 'EBPR-Control-1155', 'EBPR-Control-1240', 'EBPR-Control-1315', 'EBPR-Control-1355', 'EBPR-Control-1455')
controlTPMcounts <- finalcounts[, c(1, 9:15)]
colnames(controlTPMcounts) <- c('Locus_tag', 'EBPR-Control-1045', 'EBPR-Control-1116', 'EBPR-Control-1155', 'EBPR-Control-1240', 'EBPR-Control-1315', 'EBPR-Control-1355', 'EBPR-Control-1455')
write_delim(controlRawCounts, '/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/R3R4/results/raw_tables/R3R4_EBPRcontrol_raw_counts.tsv', delim='\t')
write_delim(controlTPMcounts, '/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/R3R4/results/tpm_tables/R3R4_EBPRcontrol')

# Early oxygen



# acetate in aerobic phase