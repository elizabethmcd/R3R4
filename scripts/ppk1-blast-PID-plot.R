library(tidyverse)
library(viridis)

# Pairwise BLAST results for ppk1 sequences from reference genomes

blast <- read.delim("results/ppk1_tree/orfs/indv_ppk1_seqs/accumulibacter-ppk1-blast-PID.tsv", header=FALSE, sep="\t")

blast_pid <- blast %>% select(V1, V2, V3)
colnames(blast_pid) <- c("ref1", "ref2", "PID")

names <- read.delim("results/ppk1_tree/orfs/indv_ppk1_seqs/ppk1-refs-list.txt", header=FALSE, sep="\t")
colnames(names) <- c("locus", "genome")
metadata <- read.csv("metadata/accumulibacter-genomes-refs-information.csv")

locus_order <- c('UW4_BBLLADNB_00506',
                 'GCA_000585075.1_LMKBKINC_00287',
                 'GCA_002455435.1_LJMNPMKE_00197',
                 'GCA_009467855.1_FFHDGLOA_00684',
                 '2687453699_PEBOPBNB_00179',
                 'GCA_000585055.1_JMBGGNHN_01742',
                 'GCA_002352265.1_FOOMNIKE_04153',
                 'GCA_000987445.1_DIPJCKAO_01419',
                 'GCA_003332265.1_BFLDHLJO_03541',
                 'UW5_GNJGOIBN_01044',
                 'GCA_000024165.1_DGOCCOHL_01070',
                 'GCA_000585035.2_BIOCIAGG_04430',
                 'GCA_000584955.2_BALNIADD_02612',
                 'GCA_002425405.1_NLEBCEFB_01950',
                 'GCA_005889575.1_AEPCJOJA_01401',
                 'GCA_000987395.1_IJAMOEGB_03109',
                 'GCA_000584975.2_FLOKBKCF_00068',
                 'UW6_NHBDKAOG_02407',
                 'GCA_000585095.1_DDBAAEFF_02735',
                 'GCA_005524045.1_CHBPCNLK_00370',
                 'UW7_CEKMHFBL_01627',
                 'GCA_002433845.1_PENGGAAI_03274',
                 'GCA_000585015.1_JPAJMACA_00169',
                 'GCA_002345285.1_NKDODKIB_03753',
                 'GCA_003542235.1_LHMBCMGL_01997',
                 'GCA_000584995.1_NPHLOHBO_00818',
                 'GCA_003535635.1_FMDBLDOA_00549',
                 'GCA_003487685.1_BPGNMNNK_01401',
                 'GCA_002345025.1_KFEAILBB_02335',
                 'GCA_001897745.1_LNFOKPNO_00497',
                 'GCA_002304785.1_GBHPAPPJ_02574',
                 'GCA_003538495.1_COGBHOPJ_04777',
                 'GCA_900089955.1_LKKGGHDO_04133')

blast_pid$ref1 <- factor(blast_pid$ref1, levels=c(locus_order))
blast_pid$ref2 <- factor(blast_pid$ref2, levels=c(locus_order))

ppk_heatmap <- blast_pid %>% ggplot(aes(x=ref1, y=ref2, fill=PID)) + geom_raster() + scale_fill_viridis(option="magma") + theme(axis.text.x= element_text(angle=85, hjust=1))

ggsave("figures/ppk_blast_PID_heatmap.png", ppk_heatmap, width=20, height=15, units=c("cm"))
