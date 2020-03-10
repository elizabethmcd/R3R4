library(tidyverse)

# coverage and mapping results of EBPR R1R2 2013 metagenomic timepoints
covg <- read.table('/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/R1R2/EBPR-MAGs/results/2013_binning/all_covg_results.txt', header = FALSE)
colnames(covg) <- c('ref', 'meta', 'avgCov')

mapping <- read.table('/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/R1R2/EBPR-MAGs/results/2013_binning/all_mapping_covg_results.txt', header=FALSE)
colnames(mapping) <- c('ref', 'meta', 'totalReads', 'mappedReads', 'avgCov')
lowqual<- c('spades-bin.2','spades-bin.36','megahit-bin.75','spades-bin.68','3300026302-bin.54','spades-bin.77', 'spades-bin.61','spades-bin.103','spades-bin.33','megahit-bin.15','spades-bin.99')

# taxonomy


# calculating relative abundance from covg
may13 <- covg %>% filter(meta == '2013-05-13-EBPR')
may23 <- covg %>% filter(meta == '2013-05-23-EBPR')
may28 <- covg %>% filter(meta == '2013-05-28-EBPR')
may13$relativeAbundance <- may13$avgCov / sum(may13$avgCov)
may23$relativeAbundance <- may23$avgCov / sum(may23$avgCov)
may28$relativeAbundance <- may28$avgCov / sum(may28$avgCov)

# avg relative abundance across 3 samples
allAbundances <- cbind.data.frame(may13$ref, may13$relativeAbundance, may23$relativeAbundance, may28$relativeAbundance)
colnames(allAbundances) <- c('ref', 'may13RA', 'may23RA', 'may28RA')
avgAbundances <- cbind.data.frame(allAbundances$ref, rowSums(allAbundances[,2:4]) / 3)
