library(tidyverse)

# Analyzing orthogroups output from OrthoFinder

# tab delimited file of orthogroups showing proteins for every genome that fell in a certain orthogroup
orthogroups <- read_delim("results/pangenomics/OrthoFinder/Results_Apr14/Orthogroups/Orthogroups.tsv", delim="\t")
colnames(orthogroups) <- c("Orthogroup", "Dechloromonas", "UW3", "UW5", "UW6", "UW1", "Thauera", "UWLDO", "UW7", "UW4")

# Accumulibacter groups - not present in outgroups
accum <- orthogroups %>% 
  filter(is.na(Dechloromonas), is.na(Thauera)) %>% 
  select(Orthogroup, UW1, UW3, UW4, UW5, UWLDO, UW6, UW7)

# present in all the clades of accumulibacter
core_accum <- accum %>% 
  filter_all(all_vars(!is.na(.)))

# select UW4 and UW6 for core
core_exp <- core_accum %>% 
  select(Orthogroup, UW4, UW6)

UW4_core <- core_exp %>% 
  select(Orthogroup, UW4) %>%
  separate_rows(UW4, sep=", ")

UW6_core <- core_exp %>% 
  select(Orthogroup, UW6) %>% 
  separate_rows(UW6, sep=", ")

# in UW4 and UW6 but doesn't necessarily have to be in all accumulibacter - "environmental-specific" core genes
experimental_core <- accum %>% 
  filter(!is.na(UW4), !is.na(UW6)) %>% 
  select(Orthogroup, UW4, UW6)

UW4_experimental_core <- experimental_core %>% 
  select(Orthogroup, UW4) %>% 
  separate_rows(UW4, sep=", ")

UW6_experimental_core <- experimental_core %>% 
  select(Orthogroup, UW6) %>% 
  separate_rows(UW6, sep=", ")

# only in UW4 - just for R3R4 genes purposes, might have to take quality into account
only_UW4 <- accum %>% 
  filter(is.na(UW1), is.na(UW5), is.na(UWLDO), is.na(UW3), is.na(UW7), is.na(UW6)) %>% 
  select(Orthogroup, UW4) %>% 
  separate_rows(UW4, sep=", ")

# only IIC from R3R4
only_UW6 <- accum %>% 
  filter(is.na(UW1), is.na(UW5), is.na(UWLDO), is.na(UW3), is.na(UW4), is.na(UW7)) %>% 
  select(Orthogroup, UW6) %>% 
  filter_all(all_vars(!is.na(.))) %>% 
  separate_rows(UW6, sep=", ")

# Accessory genes inclusive of all other clades except the reciprocal one in the experiment
# For UW6, get all genes that are in other clades but are not in UW3 and UW4, and also the outgroups
# So they either can or cannot be in other clades, don't have to be in all the other ones
uw6_inclusive <- orthogroups %>% 
  filter(is.na(Dechloromonas), is.na(Thauera), is.na(UW3), is.na(UW4)) %>% 
  filter(!is.na(UW6)) %>% 
  select(Orthogroup, UW6) %>% 
  separate_rows(UW6, sep=", ")

colnames(uw6_inclusive) <- c("Orthogroup", "locus_tag")

# For UW4, get all genes that are in other clades but are missing in IIC
uw4_inclusive <- orthogroups %>% 
  filter(is.na(Dechloromonas), is.na(Thauera), is.na(UW6)) %>% 
  filter(!is.na(UW4)) %>% 
  select(Orthogroup, UW4) %>% 
  separate_rows(UW4, sep=", ")
colnames(uw4_inclusive) <- c("Orthogroup", "locus_tag")

# write out tables
write.csv(core_exp, file="results/orthoTables/2020-05-16-core-R3R4-orthogroups.csv", quote=FALSE, row.names = FALSE)
write.csv(UW4_core, file="results/orthoTables/2020-05-16-UW4-core.csv", quote=FALSE, row.names = FALSE)
write.csv(UW6_core, file="results/orthoTables/2020-05-16-UW6-core.csv", quote=FALSE, row.names = FALSE)
write.csv(only_UW4, file="results/orthoTables/2020-05-16-onlyUW4.csv", quote=FALSE, row.names = FALSE)
write.csv(only_UW6, file="results/orthoTables/2020-05-16-onlyUW6.csv", quote=FALSE, row.names = FALSE)
write.csv(uw4_inclusive, file="results/orthoTables/2020-05-16-UW4-inclusive.csv", quote=FALSE, row.names = FALSE)
write.csv(uw6_inclusive, file="results/orthoTables/2020-05-16-UW6-inclusive.csv", quote=FALSE, row.names = FALSE)
write.csv(UW4_experimental_core, file="results/orthoTables/2020-05-16-UW4-experimental-core.csv", quote=FALSE, row.names=FALSE)
write.csv(UW6_experimental_core, file="results/orthoTables/2020-05-16-UW6-experimental-core.csv", quote=FALSE, row.names=FALSE)
