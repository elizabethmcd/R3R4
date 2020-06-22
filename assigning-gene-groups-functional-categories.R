library(tidyverse)

# KO Terms with Modules and Descriptions for Grouping R3R4 Transcriptomics Functions

# KEGG Categories & Modules Pulled Down from the API
ko_modules <- read_delim("metadata/kegg_data/kegg_modules_2019_07_23.tsv", delim="\t")
module_descriptions <- read_delim("metadata/kegg_data/kegg_module_attribute_list.tsv", delim="\t")
colnames(module_descriptions)[1] <- c("Module")

kegg_metadata <- left_join(ko_modules, module_descriptions)
colnames(kegg_metadata)[2] <- c("KO")

# R3R4 Gene Lists for Core and Accessory Pathways for UW4 & UW6 Genomes

# UW6
uw6_core <- read.csv("results/orthoTables/UW6-core.csv")
uw6_core_annotations <- left_join(uw6_core, kegg_metadata)
write.csv(uw6_core_annotations, "results/annotations/UW6-core-categories.csv", quote=FALSE, row.names = FALSE)
uw6_acc <- read.csv("results/orthoTables/UW6-accessory.csv")
uw6_acc_annotations <- left_join(uw6_acc, kegg_metadata)
write.csv(uw6_acc_annotations, "results/annotations/UW6-accessory-categories.csv", quote=FALSE, row.names = FALSE)

# UW4
uw4_core <- read.csv("results/orthoTables/UW4-core.csv")
uw4_core_annotations <- left_join(uw4_core, kegg_metadata)
write.csv(uw4_core_annotations, "results/annotations/UW4-core-categories.csv", quote=FALSE, row.names = FALSE)
uw4_acc <- read.csv("results/orthoTables/UW4-accessory.csv")
uw4_acc_annotations <- left_join(uw4_acc, kegg_metadata)
write.csv(uw4_acc_annotations, "results/annotations/UW4-accessory-categories.csv", quote=FALSE, row.names = FALSE)
