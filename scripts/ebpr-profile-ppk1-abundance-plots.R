# R3R4 time series P/acetate data and clade abundance

library(tidyverse)

metadata <- read.csv("metadata/R3-chemical-metadata.csv")

# values normalized by VSS
metadata$VSS_TSS <- metadata$VSS / metadata$TSS
metadata$P_VSS <- metadata$P / metadata$VSS
metadata$ace_VSS <- metadata$Acetate / metadata$VSS

# Timepoints within the start of the cycle
cycle <- metadata %>% 
  filter(cycle_minutes >= 0)

# Get rid of negative P values and reformat dates
cycle_cleaned <- cycle %>% 
  filter(P >= 0)

# Get cycles that have more than 10 data entries (get rid of data entries that have small amount of P tests/acetate taken for that exact minute of the cycle)
good_cycle_minutes <- cycle_cleaned %>% 
  group_by(cycle_minutes) %>% 
  summarize(count=n()) %>% 
  filter(count > 10) %>% 
  pull(cycle_minutes)

# Group by cycle minute to get average, standard deviation
P_stats <- cycle_cleaned %>%
  filter(cycle_minutes %in% good_cycle_minutes) %>% 
  group_by(cycle_minutes) %>% 
  summarize(mean_P = mean(P_VSS), max_P = max(P_VSS), min_P = min(P_VSS), med_P = median(P_VSS))

# convert to dates to select days where acetate tests were performed (sum greater than 0)
cycle_cleaned$Date <- as.Date(cycle_cleaned$Date)
acetate_days <- aggregate(cycle_cleaned$ace_VSS, by=list(cycle_cleaned$Date), sum) %>% filter(x > 0) %>% pull(Group.1)

ace_stats <- cycle_cleaned %>% 
  filter(Date %in% acetate_days) %>%
  filter(cycle_minutes %in% good_cycle_minutes) %>%
  group_by(cycle_minutes) %>% 
  summarize(mean_ace = mean(ace_VSS), max_ace = max(ace_VSS), min_ace = min(ace_VSS), med_ace = median(ace_VSS))

# separate P and Acetate plots for all sampling days across reactor operation

P_stats %>% ggplot(aes(x=cycle_minutes, y=med_P)) + geom_line() + geom_point() + geom_errorbar(aes(ymin=min_P, ymax=max_P), width=10, position=position_dodge(1))

ace_stats %>% ggplot(aes(x=cycle_minutes, y=mean_ace)) + geom_point() + geom_line() + geom_errorbar(aes(ymin=min_ace, ymax=max_ace), width=10, position=position_dodge(1))

# Data for just the sampling day in question

normal <- read.csv("metadata/R3R4-normal-cycle-chemical-data.csv")

# Normalized by VSS
normal$VSS_TSS <- normal$VSS / normal$TSS
normal$P_VSS <- normal$P / normal$VSS
normal$ace_VSS <- normal$Acetate / normal$VSS

profile <- normal %>% ggplot(aes(x=Minutes)) + geom_line(aes(y=P_VSS), colour="#2BAEB3", size=2.5) + geom_line(aes(y=ace_VSS), colour="#FF6E3C", size=2.5) + scale_x_continuous(breaks=seq(0,280, by=50), expand=c(0,0)) + scale_y_continuous(breaks=seq(0, 0.20, by=0.04), expand=c(0,0)) + geom_vline(xintercept=100, linetype="dotted", size=1) + theme_classic()

profile_transparent <- profile +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )

# qPCR data

clades <- read.csv("metadata/R3R4-accumulibacter-clade-abundance-normal-cycle.csv")

clade_order <- c("IA", "IIC", "IIA")

abundance <- clades %>% ggplot(aes(x=Clade, y=Mean, fill = as.factor(Clade))) + geom_bar(stat="identity", colour="black", size=1) + scale_fill_manual(values=c("#2BAEB3", "#404272", "#FF6E3C")) + geom_errorbar(aes(ymin=(Mean - Dev), ymax=(Mean + Dev)), width=.3, position=position_dodge(1)) + scale_y_log10(limits=c(1,1e6), expand=c(0,0), breaks = scales::log_breaks(n=7) ) + scale_x_discrete(limits = clade_order) + theme_classic()

# save plots

ggsave(plot=profile, filename="figures/R3R4-EBPR-cycle-profile.png", width=10, height=5, units=c("in"))

ggsave(plot=abundance, filename="figures/R3R4-accumulibacter-ppk1-abundance.png", width=10, height=10, units=c("in"))      

ggsave(plot= profile_transparent, filename="figures/R3R4-EBPR-cycle-profile-transparent.png", width=10, height=5, units=c("in"), bg="transparent")
