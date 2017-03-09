# set working directory
setwd("Desktop/Protocols/R programming/Kornfeld Koding Klub/Nils qPCR analysis pipeline/")
# load libraries
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)
library(grid)
library(ggplot2)
library(gridExtra)

# read CVS file
qPCR.analysis <- read_delim("2017-01-20 1-Hep 48334 LNA TFX Metabolic Genes IV_analysis/Results-Tabelle 1.csv", delim=';')

# or qPCR.analysis2 <- read.csv2("2017-01-20 1-Hep 48334 LNA TFX Metabolic Genes IV_analysis/Results-Tabelle 1.csv")

# select relevant values
qPCR.analysis.relevant <- select(qPCR.analysis, `Sample Name`, `Target Name`, `CT`)

# define vectors for targets and samples
targets <- c("Fbp1", "G6PC", "Pck1", "Slc2a1", "Slc2a2", "Slc2a4", "Acaca", "Fabp4", "Srebf1", "Acox1", "Cpt1a", "Hadhb", "IRS1", "IRS2", "Foxo1", "HPRT")
samples <- c("scr 50-1", "scr 50-2", "scr 50-3", "scr 100-1", "scr 100-2", "scr 100-3", "LNA 50-1", "LNA 50-2", "LNA 50-3", "LNA 100-1", "LNA 100-2", "LNA 100-3", "scr 50-1 + Ins", "scr 50-2 + Ins", "scr 50-3 + Ins", "scr 100-1 + Ins", "scr 100-2 + Ins", "scr 100-3 + Ins", "LNA 50-1 + Ins", "LNA 50-2 + Ins", "LNA 50-3 + Ins", "LNA 100-1 + Ins", "LNA 100-2 + Ins", "LNA 100-3 + Ins")

# rename targets and samples
qPCR.analysis.relevant$`Sample Name` <- rep(samples)
qPCR.analysis.relevant$`Target Name` <- rep(targets, each=24)

# spread by Target Name and CT
# qPCR.analysis.relevant2 <- spread(qPCR.analysis.relevant, `Target Name`, `CT`)

# data arithmetics

qPCR.analysis.relevant$deltaCT = 0

# qPCR.analysis.math <- mutate(qPCR.analysis.relevant2, deltaCT = `HPRT` - `Fbp1`)

for (i in 1:length(samples)) {
  current_sample <- samples[i]
  for (j in 1:length(targets)) {
    current_target <- targets[j]
    current_ct <- qPCR.analysis.relevant[qPCR.analysis.relevant$`Sample Name` == current_sample & qPCR.analysis.relevant$`Target Name` == current_target, 'CT' ]
    housekeeping_ct <- qPCR.analysis.relevant[qPCR.analysis.relevant$`Sample Name` == current_sample & qPCR.analysis.relevant$`Target Name` == 'HPRT', 'CT' ]
    qPCR.analysis.relevant[qPCR.analysis.relevant$`Sample Name` == current_sample & qPCR.analysis.relevant$`Target Name` == current_target, 'deltaCT' ] = housekeeping_ct - current_ct
  }
}
  
qPCR.analysis.relevant$pow_deltaCT = 0
  
# qPCR.analysis.relevant2 <- mutate(qPCR.analysis.relevant, pow_deltaCT = 2^deltaCT)

for (k in 1:length(samples)) {
  current_sample2 <- samples[k]
  for (l in 1:length(targets)) {
    current_target2 <- targets [l]
    deltaCTvalue <- qPCR.analysis.relevant[qPCR.analysis.relevant$`Sample Name` == current_sample2 & qPCR.analysis.relevant$`Target Name` == current_target2, 'deltaCT']
    qPCR.analysis.relevant[qPCR.analysis.relevant$`Sample Name` == current_sample2 & qPCR.analysis.relevant$`Target Name` == current_target2, 'pow_deltaCT'] = 2^deltaCTvalue
  }
}

qPCR.analysis.relevant$fold_change = 0

# assign control groups
# control_group1 <- c("scr 50-1", "scr 50-2", "scr 50-3") 
# control_group2 <- c("scr 100-1", "scr 100-2", "scr 100-3")
  
# divide data frame
samples_50nM <- c("scr 50-1", "scr 50-2", "scr 50-3", "LNA 50-1", "LNA 50-2", "LNA 50-3", "scr 50-1 + Ins", "scr 50-2 + Ins", "scr 50-3 + Ins", "LNA 50-1 + Ins", "LNA 50-2 + Ins", "LNA 50-3 + Ins")
samples_100nM <- c("scr 100-1", "scr 100-2", "scr 100-3", "LNA 100-1", "LNA 100-2", "LNA 100-3", "scr 100-1 + Ins", "scr 100-2 + Ins", "scr 100-3 + Ins", "LNA 100-1 + Ins", "LNA 100-2 + Ins", "LNA 100-3 + Ins")

qPCR.analysis.relevant$Molarity = 0

qPCR.analysis.relevant <- mutate(qPCR.analysis.relevant, Molarity = ifelse(qPCR.analysis.relevant$`Sample Name` %in% c("scr 50-1", "scr 50-2", "scr 50-3", "LNA 50-1", "LNA 50-2", "LNA 50-3", "scr 50-1 + Ins", "scr 50-2 + Ins", "scr 50-3 + Ins", "LNA 50-1 + Ins", "LNA 50-2 + Ins", "LNA 50-3 + Ins"), '50', '100'))

qPCR.analysis.relevant.50nM <- filter(qPCR.analysis.relevant, qPCR.analysis.relevant$Molarity == '50')
qPCR.analysis.relevant.100nM <- filter(qPCR.analysis.relevant, qPCR.analysis.relevant$Molarity == '100')

for (m in 1:length(samples_50nM)) {
  current_sample3 <- samples_50nM[m]
  for (n in 1:length(targets)) {
    current_target3 <- targets[n]
    pow_deltaCTvalue2 <- qPCR.analysis.relevant.50nM[qPCR.analysis.relevant.50nM$`Sample Name` == current_sample3 & qPCR.analysis.relevant.50nM$`Target Name` == current_target3, 'pow_deltaCT']
    pow_deltaCT_control1 <- qPCR.analysis.relevant.50nM[qPCR.analysis.relevant.50nM$`Sample Name` == 'scr 50-1' & qPCR.analysis.relevant.50nM$`Target Name` == current_target3, 'pow_deltaCT']
    pow_deltaCT_control2 <- qPCR.analysis.relevant.50nM[qPCR.analysis.relevant.50nM$`Sample Name` == 'scr 50-2' & qPCR.analysis.relevant.50nM$`Target Name` == current_target3, 'pow_deltaCT']
    pow_deltaCT_control3 <- qPCR.analysis.relevant.50nM[qPCR.analysis.relevant.50nM$`Sample Name` == 'scr 50-3' & qPCR.analysis.relevant.50nM$`Target Name` == current_target3, 'pow_deltaCT']
    pow_deltaCT_controls <- c(pow_deltaCT_control1$pow_deltaCT, pow_deltaCT_control2$pow_deltaCT, pow_deltaCT_control3$pow_deltaCT)
    pow_deltaCT_controls_mean <- mean(pow_deltaCT_controls)
    qPCR.analysis.relevant.50nM[qPCR.analysis.relevant.50nM$`Sample Name` == current_sample3 & qPCR.analysis.relevant.50nM$`Target Name` == current_target3, 'fold_change'] = pow_deltaCTvalue2/pow_deltaCT_controls_mean
  }
}

for (o in 1:length(samples_100nM)) {
  current_sample4 <- samples_100nM[o]
  for (p in 1:length(targets)) {
    current_target4 <- targets[p]
    pow_deltaCTvalue3 <- qPCR.analysis.relevant.100nM[qPCR.analysis.relevant.100nM$`Sample Name` == current_sample4 & qPCR.analysis.relevant.100nM$`Target Name` == current_target4, 'pow_deltaCT']
    pow_deltaCT_control4 <- qPCR.analysis.relevant.100nM[qPCR.analysis.relevant.100nM$`Sample Name` == 'scr 100-1' & qPCR.analysis.relevant.100nM$`Target Name` == current_target4, 'pow_deltaCT']
    pow_deltaCT_control5 <- qPCR.analysis.relevant.100nM[qPCR.analysis.relevant.100nM$`Sample Name` == 'scr 100-2' & qPCR.analysis.relevant.100nM$`Target Name` == current_target4, 'pow_deltaCT']
    pow_deltaCT_control6 <- qPCR.analysis.relevant.100nM[qPCR.analysis.relevant.100nM$`Sample Name` == 'scr 100-3' & qPCR.analysis.relevant.100nM$`Target Name` == current_target4, 'pow_deltaCT']
    pow_deltaCT_controls2 <- c(pow_deltaCT_control4$pow_deltaCT, pow_deltaCT_control5$pow_deltaCT, pow_deltaCT_control6$pow_deltaCT)
    pow_deltaCT_controls_mean2 <- mean(pow_deltaCT_controls2)
    qPCR.analysis.relevant.100nM[qPCR.analysis.relevant.100nM$`Sample Name` == current_sample4 & qPCR.analysis.relevant.100nM$`Target Name` == current_target4, 'fold_change'] = pow_deltaCTvalue3/pow_deltaCT_controls_mean2
  }
}

# Select relevant information to perform statistics and plotting
qPCR.analysis.relevant.50nM.plot <- select(qPCR.analysis.relevant.50nM, `Sample Name`, `Target Name`, `fold_change`)
qPCR.analysis.relevant.100nM.plot <- select(qPCR.analysis.relevant.100nM, `Sample Name`, `Target Name`, `fold_change`)

# Statistics

# for loop überlegen ...

summary.50 <- qPCR.analysis.relevant.50nM.plot %>% group_by(`Sample Name`) %>%
  summarise(
    count = n(),
    mean = mean(fold_change, na.rm=TRUE),
    sd = sd(fold_change, na.rm = TRUE)
  )
show(summary.50)

summary.100 <- qPCR.analysis.relevant.100nM.plot %>% group_by(`Sample Name`) %>%
  summarise(
    count = n(),
    mean = mean(fold_change, na.rm=TRUE),
    sd = sd(fold_change, na.rm = TRUE)
  )
show(summary.100)
