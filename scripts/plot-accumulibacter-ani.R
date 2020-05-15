library(tidyverse)
library(reshape2)
library(viridis)

# input from fastani
ani <- read.table("results/ppk1_tree/accum-ani.txt") %>% select(V1, V2, V3)
colnames(ani) <- c("bin1", "bin2", "ANI")

# accumulibacter clade information
accum <- read.csv("metadata/accumulibacter-genomes-refs-information.csv", stringsAsFactors = FALSE)
clade_order <- c("2687453699.fna",
                 "GCA_000585075.1.fna",
                 "GCA_002455435.1.fna",
                 "GCA_009467855.1.fna",
                 "spades-bin.32.fna",
                 "GCA_000585055.1.fna",
                 "GCA_000987445.1.fna",
                 "GCA_002352265.1.fna",
                 "GCA_003332265.1.fna",
                 "2767802316.fna",
                 "GCA_000024165.1.fna",
                 "2767802455.fna",
                 "GCA_000584955.2.fna",
                 "GCA_000584975.2.fna",
                 "GCA_000585035.2.fna",
                 "GCA_000987395.1.fna",
                 "GCA_002425405.1.fna",
                 "GCA_005889575.1.fna",
                 "GCA_000584995.1.fna",
                 "GCA_000585015.1.fna",
                 "GCA_000585095.1.fna",
                 "GCA_002345025.1.fna",
                 "GCA_002345285.1.fna",
                 "GCA_002433845.1.fna",
                 "GCA_003487685.1.fna",
                 "GCA_003535635.1.fna",
                 "GCA_003542235.1.fna",
                 "GCA_005524045.1.fna",
                 "Reactor4-bin.4.fna",
                 "GCA_001897745.1.fna",
                 "GCA_002304785.1.fna",
                 "GCA_900089955.1.fna",
                 "GCA_003538495.1.fna")

# get rid of genome with no detectable ppk1 sequence and couldn't confirm the clade identity
keep <- ani %>% filter(bin1 != "GCA_009467885.1.fna" & bin2 != "GCA_009467885.1.fna")

keep$bin1 <- factor(keep$bin1, levels=c(clade_order))
keep$bin2 <- factor(keep$bin2, levels=c(clade_order))

ani_heatmap <- keep %>% ggplot(aes(x=bin1, y=bin2, fill=ANI)) + geom_raster() + scale_fill_viridis() + theme(axis.text.x= element_text(angle=85, hjust=1))

ggsave(filename="figures/accumulibacter-ani.png", plot=ani_heatmap, width=10, height=7, units=c('cm'))

