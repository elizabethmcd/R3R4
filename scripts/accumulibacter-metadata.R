library(tidyverse)

checkm <- read.delim("metadata/all_acc_checkm_stats.tsv", sep="\t")
acc <- read.csv("metadata/accumulibacter-genomes-refs-information.csv")

colnames(checkm) <- c("genome", "phylogeny", "completion", "contamination", "size_bp", "contigs", "gc", "x")

ac_table <- left_join(checkm, acc)
write.csv(ac_table, "metadata/all-accumulibacter-metadata.csv", quote=FALSE, row.names = FALSE)

high_acc <- checkm %>% filter(completion > 90 & contamination < 5.5)
write.csv(high_acc, "metadata/high_quality_accumulibacter_metadata.csv", quote=FALSE, row.names = FALSE)
