library(tidyverse)
library(UpSetR)

ortholog_matrix <- read.delim("results/pangenomics/pyparanoid/v2/homolog_matrix.txt", sep="\t")
colnames(ortholog_matrix)[1] <- c("group")
colnames(ortholog_matrix)[2] <- c("UW3")
colnames(ortholog_matrix)[3] <- c("UW5")
colnames(ortholog_matrix)[4] <- c("UW6")

group_descriptions <- read.delim("results/pangenomics/pyparanoid/v2/group_descriptions.txt", header=FALSE, sep="\t") %>% select(V1, V2)

colnames(group_descriptions) <- c("group", "annotation")

ortholog_table <- ortholog_matrix %>% 
  column_to_rownames(var="group")

# ortholog groups where all genomes have single gene (single copy core)
core_all <- ortholog_table %>% 
  filter_all(all_vars(. == 1)) %>% 
  rownames_to_column(var="group")
core_annotations <- left_join(core_all, group_descriptions)
core_table <- left_join(core_groups, core_annotations) %>% 
  filter(annotation != "hypothetical protein")
write.csv(core_table, "results/pangenomics/pyparanoid/core_groups_table.csv", row.names = FALSE, quote=FALSE)

# ortholog group where all accumulibacter have single gene (single copy core accumulibacter)
core_accumulibacter <- ortholog_table %>%
  filter(Dechloromonas == 0 & Rhodocyclus == 0) %>% 
  select(-Dechloromonas) %>% 
  select(-Rhodocyclus) %>% 
  filter_all(all_vars(. == 1)) %>% 
  rownames_to_column(var="group")
core_accumulibacter_table <- left_join(core_accumulibacter, group_descriptions) %>% 
  filter(annotation != "hypothetical protein")

colnames(ortholog_table) <- c("UW3_IA", "UW5_IIA", "UW6_IIC", "Dechloromonas", "UW1_IIA", "SK01_IIC", "SK02_IIC", "SK12_IIF", "BA92_IB", "BA93_IA", "BA94_IIF", "HKU1_IB", "UBA2327_IIF", "UBA2315_IIF", "UBA5574_IIC", "UBA6585_IIF", "UBA6658_IA", "UWLDO_IC", "SCELSE1_IIF", "BIN19_IIC", "CANDO1_IA", "Rhodocyclus", "UW7_IIF", "UW4_IA")

# accessory IA strict
# not in any other clade, only in all IA genomes 
acc_IA <- ortholog_table %>% 
  filter(UW5_IIA == 0 & UW6_IIC == 0 & Dechloromonas == 0 & UW1_IIA == 0 & SK01_IIC == 0 & SK02_IIC == 0 & SK12_IIF == 0 & BA92_IB == 0 & BA94_IIF == 0 & HKU1_IB == 0 & UBA2327_IIF == 0 & UBA2315_IIF == 0 & UBA5574_IIC == 0 & UBA6585_IIF == 0 & UWLDO_IC == 0 & SCELSE1_IIF == 0 & BIN19_IIC == 0 & Rhodocyclus == 0 & UW7_IIF == 0) %>% 
  select(UW3_IA, UW4_IA, CANDO1_IA, UBA6658_IA, BA93_IA) %>% 
  filter_all(all_vars(. != 0)) %>% 
  rownames_to_column(var="group")
acc_IA_table <- left_join(acc_IA, group_descriptions) %>% 
  filter(annotation != "hypothetical protein")
write.csv(acc_IA_table, "results/pangenomics/pyparanoid/IA_accessory_strict.csv", row.names = FALSE, quote=FALSE)

# accessory IA inclusive
# not in any IIC genomes
acc_IA_inclusive <- ortholog_table %>% 
  filter(UW6_IIC == 0 & Dechloromonas == 0 & SK01_IIC == 0 & UBA5574_IIC == 0 & BIN19_IIC == 0 & Rhodocyclus == 0) %>% 
  filter(UW4_IA >= 1) %>% 
  select(UW4_IA) %>% 
  rownames_to_column(var="group")
IA_inclusive_table <- left_join(acc_IA_inclusive, group_descriptions) %>% 
  filter(annotation != "hypothetical protein")
write.csv(acc_IA_inclusive, "results/pangenomics/pyparanoid/IA_accessory_inclusive.csv", quote=FALSE, row.names = FALSE)

# accessory IIC strict
# not in any other clade, only in all IIC genomes
acc_IIC <- ortholog_table %>% 
  filter(UW3_IA == 0 & UW5_IIA == 0 & Dechloromonas == 0 & UW1_IIA == 0 & SK12_IIF == 0 & BA92_IB == 0 & BA93_IA == 0 & BA94_IIF == 0 & HKU1_IB == 0 & UBA2315_IIF == 0 & UBA2327_IIF == 0 & UBA6585_IIF == 0 & UBA6658_IA == 0 & UWLDO_IC == 0 & SCELSE1_IIF == 0 & CANDO1_IA == 0 & Rhodocyclus == 0 & UW7_IIF == 0 & UW4_IA == 0) %>% 
  select(UW6_IIC, SK01_IIC, SK02_IIC, UBA5574_IIC, BIN19_IIC) %>% 
  filter_all(all_vars(. != 0)) %>% 
  rownames_to_column(var="group")
acc_IIC_table <- left_join(acc_IIC, group_descriptions) %>% 
  filter(annotation != "hypothetical protein")
write.csv(acc_IIC, "results/pangenomics/pyparanoid/IIC_accessory_strict.csv", row.names = FALSE, quote=FALSE)

# accessory IIC inclusive
# not in any IA genomes
acc_IIC_inclusive <- ortholog_table %>% 
  filter(UW3_IA == 0 & BA93_IA == 0 & UBA6658_IA == 0 & CANDO1_IA == 0 & UW4_IA == 0 & Dechloromonas == 0 & Rhodocyclus == 0) %>% 
  filter(UW6_IIC >= 1) %>% 
  select(UW6_IIC) %>% 
  rownames_to_column(var="group")
IIC_inclusive_table <- left_join(acc_IIC_inclusive, group_descriptions) %>% 
  filter(annotation != "hypothetical protein")
write.csv(acc_IIC_inclusive, "results/pangenomics/pyparanoid/IIC_accessory_inclusive.csv", quote=FALSE, row.names = FALSE)

ortho_table <- rownames_to_column(ortholog_table, var="group")

# something goes wrong below here
ortho_wide_table <- ortho_table %>% pivot_longer(cols=-group, names_to="genome", values_to="count") %>% pivot_wider(names_from = group, values_from = count)
table <- ortho_wide_table
table$clade <- c("IA", "IIA", "IIC", "outgroup", "IIA", "IIC", "IIC", "IIF", "IB", "IA", "IIF", "IB", "IIF", "IIF", "IIC", "IIF", "IA", "IC", "IIF", "IIC", "IA", "outgroup", "IIF", "IA")

# upsetR plots
# has to be binary, so change anything with more than 1 gene for a genome in a group to 1
ortholog_binary <- data.frame(lapply(ortholog_table[1:ncol(ortholog_table)], function(x) ifelse(x>1, 1, x)), stringsAsFactors = FALSE)
ortholog_uw <- ortholog_binary %>%
  select(UW1_IIA, UW3_IA, UW4_IA, UW5_IIA, UW6_IIC, UW7_IIF)
# plot of UW genomes
upset(ortholog_binary, sets = c("UW1_IIA", "UW3_IA", "UW4_IA", "UW5_IIA", "UW6_IIC", "UW7_IIF", "Dechloromonas", "Rhodocyclus"), order.by="freq")

# aggregate by clade/outgroup
exclude <- ncol(table) - 1
test <- aggregate(table[2:exclude], list(table$clade), mean) 

aggregated_clade_table <- pivot_longer(test, cols=-Group.1, names_to="group", values_to="count") %>% pivot_wider(names_from=Group.1, values_from=count)
aggregated_binary_clade_table <- data.frame(lapply(aggregated_clade_table[2:ncol(aggregated_clade_table)], function(x) ifelse(x>1, 1, x)), stringsAsFactors = FALSE)
clade_table <- data.frame(lapply(aggregated_binary_clade_table[1:ncol(aggregated_binary_clade_table)], function(x) ifelse(x<1, 0, x)), stringsAsFactors = FALSE)

upset(clade_table, sets=c("outgroup", "IA", "IB", "IC", "IIA", "IIC", "IIF"), order.by="freq", empty.intersections="on")

IIA_test <- ortholog_table %>% 
  filter(UW6_IIC == 0 & Dechloromonas == 0 & SK01_IIC == 0 & SK02_IIC == 0 & SK12_IIF == 0 & BA92_IB == 0 & BA94_IIF == 0 & HKU1_IB == 0 & UBA2327_IIF == 0 & UBA2315_IIF == 0 & UBA5574_IIC == 0 & UBA6585_IIF == 0 & UWLDO_IC == 0 & SCELSE1_IIF == 0 & BIN19_IIC == 0 & Rhodocyclus == 0 & UW7_IIF == 0 & UW3_IA == 0 & UW4_IA == 0 & CANDO1_IA == 0 & UBA6658_IA == 0 & BA93_IA == 0) %>% 
  select(UW1_IIA, UW5_IIA) %>% 
  filter_all(all_vars(. >= 1))

aggregated_clade_table %>% 
  filter(IA < 1 & IB < 1 & IC < 1 & IIC < 1 & IIF < 1 & outgroup < 1) %>% 
  select(IIA) %>% 
  filter_all(all_vars(. !=0))
     
clade_table %>% 
  filter(IA == 0 & IB == 0 & IC == 0 & IIC == 0 & IIF == 0 & outgroup == 0) %>% 
  select(IIA) %>% 
  filter_all(all_vars(. !=0))

########
# group annotations
IA_groups <- read.delim("results/pangenomics/pyparanoid/core_lists/IA_groups.txt", sep="\t", header=FALSE)
IIC_groups <- read.delim("results/pangenomics/pyparanoid/core_lists/IIC_groups.txt", sep="\t", header=FALSE)
colnames(IA_groups) <- c("IA_locus_tag", "group")
colnames(IIC_groups) <- c("IIC_locus_tag", "group")
clade_groups <- left_join(IA_groups, IIC_groups)
clade_group_annotations <- left_join(clade_groups, group_descriptions) %>% select(group, annotation, IA_locus_tag, IIC_locus_tag)
