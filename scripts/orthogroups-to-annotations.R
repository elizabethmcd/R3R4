library(tidyverse)

# Merging annotations with protein content in ortholog groups

# annotation files
uw4_annotations <- read_delim("ref_MAGs/spades-bin.32/spades-bin.32.tsv", delim="\t") %>% filter(ftype == "CDS")
uw6_annotations <- read_delim("ref_MAGs/2767802455/2767802455.tsv", delim="\t") %>% filter(ftype == 'CDS')

# Core content


# Accessory content for each clade

# UW4 IA 
uw4_only <- read.csv("results/orthoTables/IA_UW4_accessory_orthogroups.csv", stringsAsFactors = FALSE)
colnames(uw4_only) <- c("Orthogroup", "locus_tag")
uw4_only_merged <- left_join(uw4_only, uw4_annotations) %>% select(Orthogroup, locus_tag, product)

# UW6 IIC
uw6_only <- read.csv("results/orthoTables/IIC_UW6_accessory_orthogroups.csv", stringsAsFactors = FALSE)
colnames(uw6_only) <- c("Orthogroup", "locus_tag")
uw6_only_merged <- left_join(uw6_only, uw6_annotations) %>% select(Orthogroup, locus_tag, product)
str(uw6_only)
