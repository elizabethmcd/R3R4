library(tidyverse)

# Analyzing orthogroups output from OrthoFinder

# tab delimited file of orthogroups showing proteins for every genome that fell in a certain orthogroup
orthogroups <- read_delim("pangenomics/OrthoFinder/Results_Apr14/Orthogroups/Orthogroups.tsv", delim="\t")
colnames(orthogroups) <- c("Orthogroup", "Dechloromonas", "UW3", "UW5", "UW6", "UW1", "Thauera", "UWLDO", "UW7", "UW4")

# core groups across all species
core <- orthogroups[complete.cases(orthogroups), ]

# Accumulibacter groups - not present in outgroups
accum <- orthogroups %>% 
  filter(is.na(Dechloromonas), is.na(Thauera)) %>% 
  select(Orthogroup, UW1, UW3, UW4, UW5, UWLDO, UW6, UW7)

# present in all the clades of accumulibacter
core_accum <- accum %>% 
  filter_all(all_vars(!is.na(.)))

# present in IIC and IA clades for transcriptomics
accum_clades <- accum %>% 
  select(Orthogroup, UW4, UW6) %>% 
  filter_all(all_vars(!is.na(.)))

# only in IA clades - includes UW3 reference from R1R2 reactors
only_IA <- accum %>% 
  filter(is.na(UW1), is.na(UW5), is.na(UWLDO), is.na(UW6), is.na(UW7)) %>% 
  select(Orthogroup, UW3, UW4) %>% 
  filter_all(all_vars(!is.na(.)))

# only in UW4 - just for R3R4 genes purposes, might have to take quality into account
only_UW4 <- accum %>% 
  filter(is.na(UW1), is.na(UW5), is.na(UWLDO), is.na(UW3), is.na(UW7), is.na(UW6)) %>% 
  select(Orthogroup, UW4)

# only IIC from R3R4
only_IIC <- accum %>% 
  filter(is.na(UW1), is.na(UW5), is.na(UWLDO), is.na(UW3), is.na(UW4), is.na(UW7)) %>% 
  select(Orthogroup, UW6) %>% 
  filter_all(all_vars(!is.na(.)))

# separate the rows in each dataframe so each protein of an orthogroup is its own row to create a nice "list" in the column
core_accum_split <- core_accum %>% 
  separate_rows(UW1, sep=",") %>%
  separate_rows(UW3, sep=",") %>% 
  separate_rows(UW4, sep=",") %>% 
  separate_rows(UW5, sep=",") %>% 
  separate_rows(UWLDO, sep=",") %>% 
  separate_rows(UW6, sep=",") %>% 
  separate_rows(UW7, sep=",")

accum_clades_split <- accum_clades %>% 
  separate_rows(UW6, sep=",") %>%
  separate_rows(UW4, sep=",")

only_UW4_split <- only_UW4 %>%
  separate_rows(UW4, sep=",")

only_IIC_split <- only_IIC %>% 
  separate_rows(UW6, sep=",")

# cbind the clades and core dataframes so the columns contain unique sets of proteins within an orthogroup

