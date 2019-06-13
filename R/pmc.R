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


iter <- 1000
BM.OU.beta.log10 <- pmc(ml.rooted.um, beta.log10, "BM", "OU", nboot = iter)
BM.OU.alpha <- pmc(ml.rooted.um, alpha, "BM", "OU", nboot = iter)
BM.GBM.beta.log10 <- pmc(ml.rooted.um, beta.log10, "BM", "trend", nboot = iter)
BM.GBM.alpha <- pmc(ml.rooted.um, alpha, "BM", "trend", nboot = iter)
# run after this
GBM.OU.beta.log10 <- pmc(ml.rooted.um, beta.log10, "trend", "OU", nboot = iter)
GBM.OU.alpha <- pmc(ml.rooted.um, alpha, "trend", "OU", nboot = iter)


#beta.log10.ll <- do.call(c, list(BM.OU.beta.log10$null, BM.OU.beta.log10$test, BM.GBM.beta.log10$null, BM.GBM.beta.log10$test, GBM.OU.beta.log10$null, GBM.OU.beta.log10$test))
#beta.log10.test <- do.call(c, list(replicate(iter*2, "BM_OU"), replicate(iter*2, "BM_GBM"), replicate(iter*2, "GBM_OU")))
#beta.log10.model <- do.call(c, list(replicate(iter, "BM"), replicate(iter, "OU"), replicate(iter, "BM"), replicate(iter, "GBM"), replicate(iter, "OU"), replicate(iter, "GBM")))
#beta.log10.df <- data.frame(llr = mttf.ll, test= mttf.test, model = mttf.model)

#alpha.ll <- do.call(c, list(BM.OU.alpha$null, BM.OU.alpha$test, BM.GBM.alpha$null, BM.GBM.alpha$test, GBM.OU.alpha$null, GBM.OU.alpha$test))
#alpha.test <- do.call(c, list(replicate(iter*2, "BM_OU"), replicate(iter*2, "BM_GBM"), replicate(iter*2, "OU_GBM")))
#alpha.model <- do.call(c, list(replicate(iter, "BM"), replicate(iter, "OU"), replicate(iter, "BM"), replicate(iter, "GBM"), replicate(iter, "OU"), replicate(iter, "GBM")))
#alpha.df <- data.frame(llr = alpha.ll, test= alpha.test, model = alpha.model)

#write.csv(mttf.df, file = "data/pmc/pmc_beta_log10.csv")
#write.csv(alpha.df, file = "data/pmc/pmc_alpha.csv")

p_value.BM.OU.beta.log10 <- length(BM.OU.beta.log10$null[BM.OU.beta.log10$null < BM.OU.beta.log10$lr] ) / iter
p_value.BM.GBM.beta.log10 <- length(BM.GBM.beta.log10$null[BM.GBM.beta.log10$null < BM.GBM.beta.log10$lr] ) / iter
p_value.GBM.OU.beta.log10 <- length(GBM.OU.beta.log10$null[GBM.OU.beta.log10$null < GBM.OU.beta.log10$lr] ) / iter

p_value.BM.OU.alpha <- length(BM.OU.alpha$null[BM.OU.alpha$null < BM.OU.alpha$lr] ) / iter
p_value.BM.GBM.alpha <- length(BM.GBM.alpha$null[BM.GBM.alpha$null < BM.GBM.alpha$lr] ) / iter
p_value.GBM.OU.alpha <- length(GBM.OU.alpha$null[GBM.OU.alpha$null < GBM.OU.alpha$lr] ) / iter


df.summary <- cbind(parameter=c("beta.log10","beta.log10","beta.log10","alpha","alpha","alpha"), 
                    test=c("BM_OU","BM_GBM", "GBM_OU", "BM_OU","BM_GBM", "GBM_OU"), 
                    llr=c(BM.OU.beta.log10$lr, BM.GBM.beta.log10$lr, GBM.OU.beta.log10$lr, BM.OU.alpha$lr, BM.GBM.alpha$lr, GBM.OU.alpha$lr),
                    p_value=c(p_value.BM.OU.beta.log10, p_value.BM.GBM.beta.log10, p_value.GBM.OU.beta.log10, p_value.BM.OU.alpha, p_value.BM.GBM.alpha, p_value.GBM.OU.alpha),
                    sigsq_BM=c(BM.OU.beta.log10$A$opt$sigsq, BM.GBM.beta.log10$A$opt$sigsq, NA, BM.OU.alpha$A$opt$sigsq, BM.GBM.alpha$A$opt$sigsq, NA),
                    sigsq_OU=c(BM.OU.beta.log10$B$opt$sigsq, NA, GBM.OU.beta.log10$B$opt$sigsq, BM.OU.alpha$B$opt$sigsq, NA,  GBM.OU.alpha$B$opt$sigsq),
                    alpha_OU=c(BM.OU.beta.log10$B$opt$alpha, NA, GBM.OU.beta.log10$B$opt$alpha, BM.OU.alpha$B$opt$alpha, NA,  GBM.OU.alpha$B$opt$alpha),
                    sigsq_GBM=c(NA, BM.GBM.beta.log10$B$opt$sigsq, GBM.OU.beta.log10$A$opt$sigsq, NA, BM.GBM.alpha$B$opt$sigsq, GBM.OU.alpha$A$opt$sigsq),
                    slope_GBM=c(NA, BM.GBM.beta.log10$B$opt$slope, GBM.OU.beta.log10$A$opt$slope, NA, BM.GBM.alpha$B$opt$slope, GBM.OU.alpha$A$opt$slope)
                    )

df.summary <- as.data.frame(df.summary)
df.summary$llr <- as.numeric(as.character(df.summary$llr))
df.summary$p_value <- as.numeric(as.character(df.summary$p_value))

write.csv(df.summary, file = "data/pmc/pmc_summary.csv")



# write code for table here
#######




