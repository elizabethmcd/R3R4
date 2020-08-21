library(tximport)
library(tximportData)
library(tidyverse)
library(reshape2)

# Control experiments
dir <- "/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/R3R4/results/transcriptomic_data"
samples <- read.table(file.path(dir, 'EBPR-Control-samples.txt'), header=FALSE)
colnames(samples) <- c('r1', 'r2', 'sample')
files <- file.path(dir, samples$sample, "abundance.h5")
names(files) <- paste0("sample", 1:7)
txi.kallisto <- tximport(files, type="kallisto", txOut=TRUE)
counts <- as.data.frame(txi.kallisto)
finalcounts <- rownames_to_column(counts, var='ID')
controlRawCounts <- finalcounts[, c(1,9:15)]
colnames(controlRawCounts) <- c('Locus_tag', 'EBPR-Control-1045', 'EBPR-Control-1116', 'EBPR-Control-1155', 'EBPR-Control-1240', 'EBPR-Control-1315', 'EBPR-Control-1355', 'EBPR-Control-1455')
controlTPMcounts <- finalcounts[, c(1:8)]
colnames(controlTPMcounts) <- c('Locus_tag', 'EBPR-Control-1045', 'EBPR-Control-1116', 'EBPR-Control-1155', 'EBPR-Control-1240', 'EBPR-Control-1315', 'EBPR-Control-1355', 'EBPR-Control-1455')
write_delim(controlRawCounts, '/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/R3R4/results/raw_tables/R3R4_EBPRcontrol_raw_counts.tsv', delim='\t')
write_delim(controlTPMcounts, '/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/R3R4/results/tpm_tables/R3R4_EBPRcontrol_tpm_counts.tsv', delim='\t')

# sum Raw counts
controlRawCounts$bin <- controlRawCounts$Locus_tag
controlRawCounts$bin <- gsub("_.*","",controlRawCounts$bin)
sumRaw <- aggregate(controlRawCounts[2:8], list(controlRawCounts$bin), sum)
sumRaw$Total <- rowSums(sumRaw[2:8])
write.csv(sumRaw, "metadata/totalCounts_AcClades.csv", quote=FALSE, row.names = FALSE)

# Aggregate and average counts across anaerobic and aerobic cycles

controlTPMcounts$bin <- controlTPMcounts$Locus_tag
controlTPMcounts$bin <- gsub("_.*","",controlTPMcounts$bin)
sumCounts = aggregate(controlTPMcounts[2:8], list(controlTPMcounts$bin), sum)
sumCounts$Total = rowSums(sumCounts[2:8])
sumCounts$Anaerobic = (sumCounts$`EBPR-Control-1045` + sumCounts$`EBPR-Control-1116` + sumCounts$`EBPR-Control-1155`) / 3
sumCounts$Aerobic = (sumCounts$`EBPR-Control-1240` + sumCounts$`EBPR-Control-1315` + sumCounts$`EBPR-Control-1355` + sumCounts$`EBPR-Control-1455`) / 4
cycles <- sumCounts %>% select(Group.1, Aerobic, Anaerobic)
colnames(cycles) <- c("Bin", "A", "N")
cycles.m <- melt(cycles, id.vars="Bin", measure.vars=c("A", "N"))

R3R4_expression <- ggplot(cycles.m, aes(x=reorder(Bin,value), y=value, fill=variable)) + geom_col(width=0.5, position="dodge") + scale_y_log10(limits=c(1,1e7), expand=c(0,0)) + coord_flip() + scale_fill_manual(values=c("grey70", "grey30")) + theme_classic()

ggsave("figures/R3R4_clades_expression_averages.png", R3R4_expression, width=15, height=10, units=c("cm"))

#############################
# N cycling genes expression profiles
nitrogen <- read.csv("results/annotations/nitrogen_cycling_annotations.csv")
nitrogen_profiles <- left_join(nitrogen, controlTPMcounts)
write.csv(nitrogen_profiles, "results/tpm_tables/nitrogen_profiles_clades.csv", quote=FALSE, row.names = FALSE)
