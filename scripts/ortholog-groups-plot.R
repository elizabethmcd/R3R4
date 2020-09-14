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
kegg_KO_modules <- read.delim("metadata/kegg_data/kegg_modules_2019_07_23.tsv", sep="\t", stringsAsFactors = FALSE)
colnames(kegg_KO_modules)[2] <- c("KO")
# kegg module attributes
kegg_attributes <- read.delim("metadata/kegg_data/kegg_module_attribute_list.tsv", sep="\t", stringsAsFactors = FALSE)
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

# deduplicated kegg descriptions
IIC_dedup_kegg <- left_join(IIC_dedup, kegg_table) %>% select(-IIC_locus_tag)
IA_dedup_kegg <- left_join(IA_dedup, kegg_table) %>% select(-IA_locus_tag)
clades_kegg_combined <- rbind(IIC_dedup_kegg, IA_dedup_kegg) %>% 
  filter_all(all_vars(. != 0))
combined_dedup <- clades_kegg_combined[!duplicated(clades_kegg_combined$group), ]

matrix = read.delim("results/pangenomics/pyparanoid/v2/homolog_matrix.txt", stringsAsFactors = FALSE, sep="\t")
colnames(matrix) = c("group", "UW3_IA", "UW5_IIA", "UW6_IIC", "Dechloromonas", "UW1_IIA", "SK01_IIC", "SK02_IIC", "SK12_IIF", "BA92_IB", "BA93_IA", "BA94_IIF", "HKU1_IB", "UBA2327_IIF", "UBA2315_IIF", "UBA5574_IIC", "UBA6585_IIF", "UBA6658_IA", "UWLDO_IC", "SCELSE1_IIF", "BIN19_IIC", "CANDO1_IA", "Rhodocyclus", "UW7_IIF", "UW4_IA")

matrix_groups <- matrix %>% select(group, UW3_IA)
all_group_annotations <- left_join(matrix_groups, combined_dedup) %>% 
  select(group, General_category)
all_group_annotations$General_category[is.na(all_group_annotations$Genera_category)] <- "None"

matrix_annotations <- left_join(matrix, all_group_annotations)
matrix_subset <- matrix_annotations %>% 
  select(group, UW1_IIA, UW5_IIA, UW6_IIC, UW7_IIF, UW3_IA, UW4_IA, UWLDO_IC, Dechloromonas, Rhodocyclus ) %>% 
  column_to_rownames(var="group") %>% 
  filter_all(any_vars(. != 0)) %>% 
  rownames_to_column(var="group")

###
# annotated KEGG groups don't look great because most of the groups don't have a KEGG assignment, so leave without annotations

matrix_binary = data.frame(lapply(matrix_subset[2:ncol(matrix_subset)], function(x) ifelse(x>1, 1, x)), list(matrix_subset$group), stringsAsFactors = FALSE)
colnames(matrix_binary)[10] <- c("group")
matrix_binary_annotations <- left_join(matrix_binary, all_group_annotations)

uw_clades = colnames(matrix_binary)[1:9]
clade_upset_plot <- upset(matrix_binary, uw_clades, width_ratio=8, set_sizes=FALSE, keep_empty_groups=FALSE, min_size=45, sort_sets=FALSE, base_annotations=list(
  'Intersection size'=intersection_size(counts=FALSE)), stripes=c('purple4', 'purple4', 'darkblue', 'darkolivegreen4','steelblue4', 'steelblue4', 'tan3', 'gray60', 'gray60'))
clade_upset_plot_numbers <- upset(matrix_binary, uw_clades, width_ratio=8, set_sizes=FALSE, keep_empty_groups=FALSE, min_size=45, sort_sets=FALSE, stripes=c('purple4', 'purple4', 'darkblue', 'darkolivegreen4','steelblue4', 'steelblue4', 'tan3', 'gray60', 'gray60'))

ggsave("figures/UW-clade-upset-ortholog-plot.png", clade_upset_plot, width=25, height=15, units=c("cm"))
ggsave("figures/UW-clade-upset-raw.png", clade_upset_plot_numbers, width=25, height=15, units=c("cm"))
