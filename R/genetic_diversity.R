#### trait regression
rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('phylolm')
library('ggplot2')
library("ape")
library('latex2exp')
library('ggpubr')

df.species <- read.table("data/demography/weibull_results_clean_species.csv", 
                         header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
rownames(df.species) <- df.species$Species
div <-  read.table("data/breseq/genetic_diversity.txt",sep = "\t",header = TRUE)
div$Pi.log10 <- log10(div$Pi)
div$Theta.log10 <- log10(div$Theta)
div$TD.log10 <- log10(div$Tajimas_D)

# get mean of log 10 transformed observations
div.mean <- aggregate(div[, c('Pi.log10', 'Theta.log10', 'TD.log10', 'dN_dS_total', 'dN_dS_fixed')], list(div$Species), mean)
colnames(div.mean)[1] <- "Species"
colnames(div.mean)[2] <- "Pi.log10"
colnames(div.mean)[3] <- "Theta.log10"
colnames(div.mean)[4] <- "TD.log10"
colnames(div.mean)[5] <- "dN_dS_total"
colnames(div.mean)[6] <- "dN_dS_fixed"

df.species <- merge(div.mean, df.species,by="Species")
plot(df.species$dN_dS_fixed, df.species$alpha)

# make plots



# plot diversity figures in supplement

