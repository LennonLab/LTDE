rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

# unblock to install package
#library("devtools")
#install_github("cboettig/pmc")

library("pmc")
library("ape")

df <- read.table("data/demography/weibull_results_clean_species.csv", 
                 header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
df<-df[!(df$Species=="KBS0812"),]
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
ml.rooted.um.prunned <- drop.tip(ml.rooted.um, 
                                 ml.rooted.um$tip.label[na.omit(match(c('KBS0812'),
                                                                      ml.rooted.um$tip.label))])
mttf <- df$mttf
names(mttf) <- df$Species
alpha <- df$alpha
names(alpha) <- df$Species
iter <- 1000
pmc.mttf <- pmc(ml.rooted.um.prunned, mttf, "BM", "OU", nboot = iter)
pmc.alpha <- pmc(ml.rooted.um.prunned, alpha, "BM", "OU", nboot = iter)

pmc.mttf.df <- data.frame(null = pmc.mttf$null, test = pmc.mttf$test)
colnames(pmc.mttf.df) <- c("BM", "OU")
write.csv(pmc.mttf.df, file = "data/pmc/pmc_mttf.csv")

pmc.alpha.df <- data.frame(null = pmc.alpha$null, test = pmc.alpha$test)
colnames(pmc.alpha.df) <- c("BM", "OU")
write.csv(pmc.alpha.df, file = "data/pmc/pmc_alpha.csv")

df.summary <- cbind(Parameter=c("mttf","alpha"), LikelihoodRatio=c(pmc.mttf$lr,pmc.alpha$lr))
df.summary <- as.data.frame(df.summary)
write.csv(df.summary, file = "data/pmc/pmc_summary.csv")




