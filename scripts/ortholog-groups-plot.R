library(ggplot2)
library(ggplot2movies)
library(ComplexUpset)
library(tidyverse)

# annotations
# IA
IA_locus_tag_groups <- read.delim("results/pangenomics/pyparanoid/v2/IA-locus-tag-groupts.txt", header=FALSE, stringsAsFactors = FALSE, sep="\t")
colnames(IA_locus_tag_groups) <- c("locus_tag", "group")
IA_locus_tag_kegg <- read.delim("results/annotations/IA-locus-tags-KOs.txt", header=FALSE, stringsAsFactors = FALSE, sep="\t")
colnames(IA_locus_tag_kegg) <- c("locus_tag", "KO")
IA_groups_kegg <- left_join(IA_locus_tag_groups, IA_locus_tag_kegg)
# IIC
IIC_locus_tag_groups <- read.delim("results/pangenomics/pyparanoid/v2/IIC-locus-tag-groups.txt", header=FALSE, stringsAsFactors = FALSE, sep="\t")
colnames(IIC_locus_tag_groups) <- c("locus_tag", "group")
IIC_locus_tag_kegg <- read.delim("results/annotations/IIC-locus-tags-KOs.txt", header=FALSE, stringsAsFactors = FALSE, sep="\t")
colnames(IIC_locus_tag_kegg) <- c("locus_tag", "KO")
IIC_groups_kegg <- left_join(IIC_locus_tag_groups, IIC_locus_tag_kegg)


# kegg KO to modules
kegg_KO_modules <- read.delim("metadata/kegg_data/kegg_modules_2019_07_23.tsv", sep="\t")
colnames(kegg_KO_modules)[2] <- c("KO")
# kegg module attributes
kegg_attributes <- read.delim("metadata/kegg_data/kegg_module_attribute_list.tsv", sep="\t")
colnames(kegg_attributes)[1] <- c("Module")
kegg_table <- left_join(kegg_KO_modules, kegg_attributes)

# merged with IA and IIC
IA_groups_kegg_modules <- left_join(IA_groups_kegg, kegg_table)
IIC_groups_kegg_modules <- left_join(IIC_groups_kegg, kegg_table)

# merge IA and IIC and then kegg descriptions
colnames(IA_groups_kegg)[1] <- c("IA_locus_tag")
colnames(IIC_groups_kegg)[1] <- c("IIC_locus_tag")

IA_dedup <- IA_groups_kegg[!duplicated(IA_groups_kegg$group), ]
IIC_dedup <- IIC_groups_kegg[!duplicated(IIC_groups_kegg$group), ]
clades_merged_ko <- left_join(IA_dedup, IIC_dedup, by="group")
colnames(clades_merged_ko) <- c("IA_locus_tag", "group", "KO_IA", "IIC_locus_tag", "KO_IIC")

matrix = read.delim("results/pangenomics/pyparanoid/v2/homolog_matrix.txt", stringsAsFactors = FALSE, sep="\t")
colnames(matrix) = c("group", "UW3_IA", "UW5_IIA", "UW6_IIC", "Dechloromonas", "UW1_IIA", "SK01_IIC", "SK02_IIC", "SK12_IIF", "BA92_IB", "BA93_IA", "BA94_IIF", "HKU1_IB", "UBA2327_IIF", "UBA2315_IIF", "UBA5574_IIC", "UBA6585_IIF", "UBA6658_IA", "UWLDO_IC", "SCELSE1_IIF", "BIN19_IIC", "CANDO1_IA", "Rhodocyclus", "UW7_IIF", "UW4_IA")
matrix_binary = data.frame(lapply(matrix[2:ncol(matrix)], function(x) ifelse(x>1, 1, x)), stringsAsFactors = FALSE) %>% 
  select(UW1_IIA, UW5_IIA, UW6_IIC, UW7_IIF, UW3_IA, UW4_IA, UWLDO_IC, Dechloromonas, Rhodocyclus ) %>% 
  filter_all(any_vars(. != 0))
uw_clades = colnames(matrix_binary)[1:9]
upset(matrix_binary, uw_clades, width_ratio=8, set_sizes=FALSE, keep_empty_groups=FALSE, min_size=45, sort_sets=FALSE, stripes=c('darkblue', 'purple4', 'purple4', 'darkolivegreen4','steelblue4', 'steelblue4', 'tan3', 'gray60', 'gray60'))
