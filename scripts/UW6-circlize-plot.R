library(circlize)
library(tidyverse)

# files
gc_info <- read.table("results/genome_track/gc_calculations_oriented.txt")
colnames(gc_info) <- c("gc_content", "gc_skew", "gc_culm")
cds_pos <- read.table("results/genome_track/UW6-cds-positive-modf.bed", header=TRUE)
cds_neg <- read.table("results/genome_track/UW6-cds-negative-modf.bed", header=TRUE)
trna <- read.table("results/genome_track/UW6-tRNA-modf.bed", header=TRUE)
rrna <- read.table("results/genome_track/UW6-rRNA-modf.bed", header=TRUE)
cytoband.df <- read.table("results/genome_track/UW6_cytoband_file.txt")

# check cds
cds_pos$size <- cds_pos$end - cds_pos$start
cds_post_modf <- cds_pos %>% filter(size > 0) %>% select(chr, start, end, cds_pos)

# illumina coverage files
covg <- read.table("results/genome_track/UW6.regions.modified.bed", sep="\t")

# bed with gc and coverage together
final_bed <- covg
colnames(final_bed) <- c("chr", "start", "end", "covg")
final_bed$gc_content <- gc_info$gc_content
final_bed$gc_scew <- gc_info$gc_skew
final_bed$gc_culm <- gc_info$gc_culm

# plotting initial track
circos.clear()
circos.par(start.degree=75, gap.after=30, cell.padding=c(0,0,0,0))
circos.initializeWithIdeogram(cytoband.df, plotType=NULL)
circos.track(ylim=c(0,1), panel.fun=function(x,y){
}, track.height=0.05, bg.col="#FF6E3C")

# maximum ranges to display
maximum_illumina <- mean(final_bed$covg) * 1.2
minimum_illumina <- min(final_bed$covg)

# adjust max for plotting purposes
final_bed$covg[which(final_bed$covg > maximum_illumina)] <- maximum_illumina

# coverage adjusted track
circos.genomicTrackPlotRegion(final_bed, numeric.column=c("covg") , ylim =c(minimum_illumina, maximum_illumina), panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, numeric.column="covg", baseline = 0, col = "#404272", border=NA, area=TRUE, ...)
}, track.height = 0.2, bg.border=NA)

# gc_content track
circos.genomicTrackPlotRegion(final_bed, numeric.column="gc_content" , ylim = c(min(final_bed$gc_content),max(final_bed$gc_content)), panel.fun = function(region, value, ...) {
  
  circos.genomicLines(region, value, numeric.column="gc_content", baseline = 0, col = rgb(0,0,0,0.7), border = NA, ...)
  
},track.height = 0.1)

final_bed$gc_skew <- gc_info$gc_skew / max(gc_info$gc_skew)
final_bed$gc_culm <- gc_info$gc_culm / max(abs(gc_info$gc_culm))

# plot GC skew and culmulative skew on the same track.
circos.genomicTrackPlotRegion(final_bed, numeric.column="gc_culm", ylim = c(-1, 1), panel.fun = function(region, value, ...) {
  
  circos.genomicLines(region, value, numeric.column="gc_culm", baseline = 0, col = rgb(0,0,0,0.5), border = NA, area=TRUE,...)
  
  circos.genomicLines(region, value, numeric.column="gc_skew", baseline = 0, col = rgb(0,0,0,1), border = NA,...)
  
}, bg.border = NA, track.height = 0.1)

# plot positive coding sequence track
circos.genomicTrackPlotRegion(cds_post_modf, numeric.column="cds_pos", ylim = c(0, 1), panel.fun = function(region, value, ...) {
  
  # in every region, print a vertical line at the start region.
  for (k in seq_len(nrow(region))){
    # plot vertical lines for each START of the ORF from bottom to top.
    circos.lines(rep(region[k, 1], 2), c(0, 1), lwd = 0.5, straight = TRUE, col = rgb(0,0,0,0.5))
  }
  
}, bg.border = NA, track.height = 0.05)

# plot negative coding sequence track
circos.genomicTrackPlotRegion(cds_neg, numeric.column="cds_neg", ylim = c(0, 1), panel.fun = function(region, value, ...) {
  
  # in every region, print a vertical line at the start region.
  for (k in seq_len(nrow(region))){
    # plot vertical lines for each START of the orf from bottom to top.
    circos.lines(rep(region[k, 1], 2), c(0, 1), lwd = 0.5, straight = TRUE, col = rgb(0,0,0,0.5))
  }
  
}, bg.border = NA, track.height = 0.05)

# plot tRNA and rRNA track
circos.genomicTrackPlotRegion(list(trna, rrna), ylim = c(0, 1), panel.fun = function(region, value, ...) {
  
  # this will iterate through the list of bed dataframes (trna and rrna)
  i = getI(...)
  
  # plot tRNAs
  if (i == 1){
    
    for (k in seq_len(nrow(region))){
      # plot vertical lines for each START of the orf from bottom to top.
      circos.lines(rep(region[k, 1], 2), c(0, 1), lwd = 2, straight = TRUE, col = rgb(0,0,0,0.2))
    }
    
  } else {
    # plot rRNAs as red
    for (k in seq_len(nrow(region))){
      # plot vertical lines for each START of the orf from bottom to top.
      circos.lines(rep(region[k, 1], 2), c(0, 1), lwd = 2, straight = TRUE, col = rgb(1,0,0,1))
    }
  }
}, bg.border = NA, track.height = 0.05)
