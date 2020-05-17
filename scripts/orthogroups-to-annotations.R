library(tidyverse)

# Merging annotations with locus tags in ortholog groups

# annotation files from prokka
uw4_prokka <- read_delim("ref_MAGs/spades-bin.32/spades-bin.32.tsv", delim="\t") %>% filter(ftype == "CDS")
uw6_prokka <- read_delim("ref_MAGs/2767802455/2767802455.tsv", delim="\t") %>% filter(ftype == 'CDS')

# KEGG annotations from KofamKOALA
uw6_kegg <- read_delim("results/orthoTables/2767802455-kofam-annotations-sig-formatted.txt", delim="\t", col_names=c("locus_tag", "KO", "Kannotation"))
uw4_kegg <- read_delim("results/orthoTables/spades-bin.32-kofam-annotations-sig-formatted.txt", delim="\t", col_names=c("locus_tag", "KO", "Kannotation"))

# Core content
# Split by clade to merge with annotations

# UW4 all core
UW4_core_orthogroups <- read.csv("results/orthoTables/2020-05-16-UW4-core.csv")
colnames(UW4_core_orthogroups) <- c("orthogroup", "locus_tag")
UW4_core_og_prokka <- left_join(UW4_core_orthogroups, uw4_prokka) %>% select(orthogroup, locus_tag, product) %>% 
  filter(product != "hypothetical protein")
UW4_core_og_kegg <- left_join(UW4_core_orthogroups, uw4_kegg) %>% 
  select(orthogroup, locus_tag, KO, Kannotation) %>% 
  filter(!is.na(KO))

# UW4 experimental core
UW4_core_exp_ogs <- read.csv("results/orthoTables/2020-05-16-UW4-experimental-core.csv")
colnames(UW4_core_exp_ogs) <- c("orthogroup", "locus_tag")
UW4_core_exp_prokka <- left_join(UW4_core_exp_ogs, uw4_prokka) %>% 
  select(orthogroup, locus_tag, product) %>% 
  filter(product != "hypothetical protein")
UW4_core_exp_kegg <- left_join(UW4_core_exp_ogs, uw4_kegg) %>% 
  select(orthogroup, locus_tag, KO, Kannotation) %>% 
  filter(!is.na(KO))

# UW6 all core
UW6_core_ogs <- read.csv("results/orthoTables/2020-05-16-UW6-core.csv")
colnames(UW6_core_ogs) <- c("orthogroup", "locus_tag")
UW6_core_ogs_prokka <- left_join(UW6_core_ogs, uw6_prokka) %>% 
  select(orthogroup, locus_tag, product) %>% 
  filter(product != "hypothetical protein")
UW6_core_ogs_kegg <- left_join(UW6_core_ogs, uw6_kegg) %>% 
  select(orthogroup, locus_tag, KO, Kannotation) %>% 
  filter(!is.na(KO))

# UW6 experimental core 
UW6_core_exp_ogs <- read.csv("results/orthoTables/2020-05-16-UW6-experimental-core.csv")
colnames(UW6_core_exp_ogs) <- c("orthogroup", "locus_tag")
Uw6_core_exp_prokka <- left_join(UW6_core_exp_ogs, uw6_prokka) %>% 
  select(orthogroup, locus_tag, product) %>% 
  filter(product != "hypothetical protein")
UW6_core_exp_kegg <- left_join(UW6_core_exp_ogs, uw6_kegg) %>% 
  select(orthogroup, locus_tag, KO, Kannotation) %>% 
  filter(!is.na(KO))

# Accessory content for each clade

# UW4 accessory content
# Only in UW4 and nothing else
# in UW4 and not UW6
UW4_only_ogs <- read.csv("results/orthoTables/2020-05-16-onlyUW4.csv")
colnames(UW4_only_ogs) <- c("orthogroup", "locus_tag")
UW4_only_prokka <- left_join(UW4_only_ogs, uw4_prokka) %>% 
  select(orthogroup, locus_tag, product) %>% 
  filter(product != "hypothetical protein")
UW4_only_kegg <- left_join(UW4_only_ogs, uw4_kegg) %>% 
  select(orthogroup, locus_tag, KO, Kannotation) %>% 
  filter(!is.na(KO))

UW4_inclusive_ogs <- read.csv("results/orthoTables/2020-05-16-UW4-inclusive.csv")
colnames(UW4_inclusive_ogs) <- c("orthogroup", "locus_tag")
UW4_inclusive_prokka <- left_join(UW4_inclusive_ogs, uw4_prokka) %>% 
  select(orthogroup, locus_tag, product) %>% 
  filter(product != "hypothetical protein")
UW4_inclusive_kegg <- left_join(UW4_inclusive_ogs, uw4_kegg) %>% 
  select(orthogroup, locus_tag, KO, Kannotation) %>% 
  filter(!is.na(KO))

# UW6 accessory content
# Only in UW6 and nothing else
# in UW6 and not UW4
UW6_only_ogs <- read.csv("results/orthoTables/2020-05-16-onlyUW6.csv")
colnames(UW6_only_ogs) <- c("orthogroup", "locus_tag")
UW6_only_prokka <- left_join(UW6_only_ogs, uw6_prokka) %>% 
  select(orthogroup, locus_tag, product) %>% 
  filter(product != "hypothetical protein")
UW6_only_kegg <- left_join(UW6_only_ogs, uw6_kegg) %>% 
  select(orthogroup, locus_tag, KO, Kannotation) %>% 
  filter(!is.na(KO))

UW6_inclusive_ogs <- read.csv("results/orthoTables/2020-05-16-UW6-inclusive.csv")
colnames(UW6_inclusive_ogs) <- c("orthogroup", "locus_tag")
UW6_inclusive_prokka <- left_join(UW6_inclusive_ogs, uw6_prokka) %>% 
  select(orthogroup, locus_tag, product) %>% 
  filter(product != "hypothetical protein")
UW6_inclusive_kegg <- left_join(UW6_inclusive_ogs, uw6_kegg) %>% 
  select(orthogroup, locus_tag, KO, Kannotation) %>% 
  filter(!is.na(KO))
