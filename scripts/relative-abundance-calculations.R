library(tidyverse)
library(formattable)
library(htmltools)
library(webshot)

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

###### R3R4 genome information
acc_metadata <- read.csv("metadata/high_quality_accumulibacter_metadata.csv", stringsAsFactors = FALSE) %>% 
  select(genome, completion, contamination, size_bp, contigs, gc) %>% filter(genome == "2767802316" | genome == "2767802455" | genome == "Reactor4-bin.4" | genome == "spades-bin.32")
colnames(acc_metadata) <- c("Genome", "Completion", "Contamination", "Size_Mbp", "Contigs", "GC")
acc_metadata$Size_Mbp <- acc_metadata$Size_Mbp / 1000000
acc_metadata$Size_Mbp <- round(acc_metadata$Size_Mbp, digits=2)
acc_metadata[1,1] <- c("IIA-UW5")
acc_metadata[2,1] <- c("IIC-UW6")
acc_metadata[3,1] <- c("IIF-UW7")
acc_metadata[4,1] <- c("IA-UW4")
R1R2_genomes <- formattable(acc_metadata)

export_formattable <- function(f, file, width = "60%", height = NULL,
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}
export_formattable(R1R2_genomes, "figures/R1R2-genomes-metadata-table.pdf")
