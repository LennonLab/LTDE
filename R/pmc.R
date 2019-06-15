rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

# unblock to install package
#library("devtools")
#install_github("cboettig/pmc")

library("pmc")
library("ape")
library("ggplot2")
library("reshape")
library("latex2exp")
library('ggpubr')

library("tidyr")
library("dplyr")

df <- read.table("data/demography/weibull_results_clean_species.csv", 
                 header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
df<-df[!(df$Species=="KBS0727"),]
rownames(df) <- df$Species

# Load ML tree
ml.tree <- read.tree("data/tree/RAxML_bipartitionsBranchLabels.ltde")
# Define the outgroup
outgroup <- match("NC_005042.1.353331-354795", ml.tree$tip.label)
# Create a rooted tree {ape} using the outgroup
ml.rooted <- root(ml.tree, outgroup, resolve.root = TRUE)
# Keep rooted but drop outgroup branch
ml.rooted <- drop.tip(ml.rooted, c("NC_005042.1.353331-354795"))
is.ultrametric(ml.rooted)
ml.rooted.um  <- chronos(ml.rooted)
is.ultrametric(ml.rooted.um)

beta.log10 <- df$beta.log10
names(beta.log10) <- df$Species
alpha <- df$alpha
names(alpha) <- df$Species
mttf.log10 <- df$mttf.log10
names(mttf.log10) <- df$Species


iter <- 1000
BM.OU.beta.log10 <- pmc(ml.rooted.um, beta.log10, "BM", "OU", nboot = iter)
BM.OU.alpha <- pmc(ml.rooted.um, alpha, "BM", "OU", nboot = iter)
BM.OU.mttf.log10 <- pmc(ml.rooted.um, mttf.log10, "BM", "OU", nboot = iter)

BM.GBM.beta.log10 <- pmc(ml.rooted.um, beta.log10, "BM", "trend", nboot = iter)
BM.GBM.alpha <- pmc(ml.rooted.um, alpha, "BM", "trend", nboot = iter)
BM.GBM.mttf.log10 <- pmc(ml.rooted.um, mttf.log10, "BM", "trend", nboot = iter)

# run after this
GBM.OU.beta.log10 <- pmc(ml.rooted.um, beta.log10, "trend", "OU", nboot = iter)
GBM.OU.alpha <- pmc(ml.rooted.um, alpha, "trend", "OU", nboot = iter)
GBM.OU.mttf.log10 <- pmc(ml.rooted.um, mttf.log10, "trend", "OU", nboot = iter)

p_value.BM.OU.beta.log10 <- length(BM.OU.beta.log10$null[BM.OU.beta.log10$null > BM.OU.beta.log10$lr] ) / iter
p_value.BM.GBM.beta.log10 <- length(BM.GBM.beta.log10$null[BM.GBM.beta.log10$null > BM.GBM.beta.log10$lr] ) / iter
p_value.GBM.OU.beta.log10 <- length(GBM.OU.beta.log10$null[GBM.OU.beta.log10$null > GBM.OU.beta.log10$lr] ) / iter

p_value.BM.OU.alpha <- length(BM.OU.alpha$null[BM.OU.alpha$null > BM.OU.alpha$lr] ) / iter
p_value.BM.GBM.alpha <- length(BM.GBM.alpha$null[BM.GBM.alpha$null > BM.GBM.alpha$lr] ) / iter
p_value.GBM.OU.alpha <- length(GBM.OU.alpha$null[GBM.OU.alpha$null > GBM.OU.alpha$lr] ) / iter

p_value.BM.OU.mttf.log10 <- length(BM.OU.mttf.log10$null[BM.OU.mttf.log10$null > BM.OU.mttf.log10$lr] ) / iter
p_value.BM.GBM.mttf.log10 <- length(BM.GBM.mttf.log10$null[BM.GBM.mttf.log10$null > BM.GBM.mttf.log10$lr] ) / iter
p_value.GBM.OU.mttf.log10 <- length(GBM.OU.mttf.log10$null[GBM.OU.mttf.log10$null > GBM.OU.mttf.log10$lr] ) / iter



df.summary <- cbind(parameter=c("beta.log10","beta.log10","beta.log10","alpha","alpha","alpha", "mttf.log10", "mttf.log10", "mttf.log10"), 
                    test=c("BM_OU","BM_GBM", "GBM_OU", "BM_OU","BM_GBM", "GBM_OU", "BM_OU","BM_GBM", "GBM_OU"), 
                    llr=c(BM.OU.beta.log10$lr, BM.GBM.beta.log10$lr, GBM.OU.beta.log10$lr, BM.OU.alpha$lr, BM.GBM.alpha$lr, GBM.OU.alpha$lr, BM.OU.mttf.log10$lr, BM.GBM.mttf.log10$lr, GBM.OU.mttf.log10$lr),
                    p_value=c(p_value.BM.OU.beta.log10, p_value.BM.GBM.beta.log10, p_value.GBM.OU.beta.log10, p_value.BM.OU.alpha, p_value.BM.GBM.alpha, p_value.GBM.OU.alpha, p_value.BM.OU.mttf.log10, p_value.BM.GBM.mttf.log10, p_value.GBM.OU.mttf.log10),
                    sigsq_BM=c(BM.OU.beta.log10$A$opt$sigsq, BM.GBM.beta.log10$A$opt$sigsq, NA, BM.OU.alpha$A$opt$sigsq, BM.GBM.alpha$A$opt$sigsq, NA, BM.OU.mttf.log10$A$opt$sigsq, BM.GBM.mttf.log10$A$opt$sigsq, NA),
                    sigsq_OU=c(BM.OU.beta.log10$B$opt$sigsq, NA, GBM.OU.beta.log10$B$opt$sigsq, BM.OU.alpha$B$opt$sigsq, NA,  GBM.OU.alpha$B$opt$sigsq, BM.OU.mttf.log10$B$opt$sigsq, NA,  GBM.OU.mttf.log10$B$opt$sigsq),
                    alpha_OU=c(BM.OU.beta.log10$B$opt$alpha, NA, GBM.OU.beta.log10$B$opt$alpha, BM.OU.alpha$B$opt$alpha, NA,  GBM.OU.alpha$B$opt$alpha, BM.OU.mttf.log10$B$opt$alpha, NA,  GBM.OU.mttf.log10$B$opt$alpha),
                    sigsq_GBM=c(NA, BM.GBM.beta.log10$B$opt$sigsq, GBM.OU.beta.log10$A$opt$sigsq, NA, BM.GBM.alpha$B$opt$sigsq, GBM.OU.alpha$A$opt$sigsq, NA, BM.GBM.mttf.log10$B$opt$sigsq, GBM.OU.mttf.log10$A$opt$sigsq),
                    slope_GBM=c(NA, BM.GBM.beta.log10$B$opt$slope, GBM.OU.beta.log10$A$opt$slope, NA, BM.GBM.alpha$B$opt$slope, GBM.OU.alpha$A$opt$slope, NA, BM.GBM.mttf.log10$B$opt$slope, GBM.OU.mttf.log10$A$opt$slope)
                    )

df.summary <- as.data.frame(df.summary)
df.summary$llr <- as.numeric(as.character(df.summary$llr))
df.summary$p_value <- as.numeric(as.character(df.summary$p_value))

write.csv(df.summary, file = "data/pmc/pmc_summary.csv")



# write code for table here
#######




