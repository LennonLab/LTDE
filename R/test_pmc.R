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
mttf <- log10(df$mttf)
names(mttf) <- df$Species
alpha <- df$alpha
names(alpha) <- df$Species
iter <- 1000
BM.OU.mttf <- pmc(ml.rooted.um.prunned, mttf, "BM", "OU", nboot = iter)
BM.OU.alpha <- pmc(ml.rooted.um.prunned, alpha, "BM", "OU", nboot = iter)
BM.GBM.mttf <- pmc(ml.rooted.um.prunned, mttf, "BM", "trend", nboot = iter)
BM.GBM.alpha <- pmc(ml.rooted.um.prunned, alpha, "BM", "trend", nboot = iter)
# run after this
OU.GBM.mttf <- pmc(ml.rooted.um.prunned, mttf, "OU", "trend", nboot = iter)
OU.GBM.alpha <- pmc(ml.rooted.um.prunned, alpha, "OU", "trend", nboot = iter)

mttf.ll <- do.call(c, list(BM.OU.mttf$null, BM.OU.mttf$test, BM.GBM.mttf$null, BM.GBM.mttf$test, OU.GBM.mttf$null, OU.GBM.mttf$test))
mttf.test <- do.call(c, list(replicate(iter*2, "BM_OU"), replicate(iter*2, "BM_GBM"), replicate(iter*2, "OU_GBM")))
mttf.model <- do.call(c, list(replicate(iter, "BM"), replicate(iter, "OU"), replicate(iter, "BM"), replicate(iter, "GBM"), replicate(iter, "OU"), replicate(iter, "GBM")))
mttf.df <- data.frame(llr = mttf.ll, test= mttf.test, model = mttf.model)

alpha.ll <- do.call(c, list(BM.OU.alpha$null, BM.OU.alpha$test, BM.GBM.alpha$null, BM.GBM.alpha$test, OU.GBM.alpha$null, OU.GBM.alpha$test))
alpha.test <- do.call(c, list(replicate(iter*2, "BM_OU"), replicate(iter*2, "BM_GBM"), replicate(iter*2, "OU_GBM")))
alpha.model <- do.call(c, list(replicate(iter, "BM"), replicate(iter, "OU"), replicate(iter, "BM"), replicate(iter, "GBM"), replicate(iter, "OU"), replicate(iter, "GBM")))
alpha.df <- data.frame(llr = alpha.ll, test= alpha.test, model = alpha.model)

write.csv(mttf.df, file = "data/pmc/pmc_mttf.csv")
write.csv(alpha.df, file = "data/pmc/pmc_alpha.csv")

df.summary <- cbind(parameter=c("mttf","mttf","mttf","alpha","alpha","alpha"), 
                    test=c("BM_OU","BM_GBM", "OU_GBM", "BM_OU","BM_GBM", "OU_GBM"), 
                    llr=c(BM.OU.mttf$lr, BM.GBM.mttf$lr, OU.GBM.mttf$lr, BM.OU.alpha$lr, BM.GBM.alpha$lr, OU.GBM.alpha$lr))
df.summary <- as.data.frame(df.summary)
write.csv(df.summary, file = "data/pmc/pmc_summary.csv")


# get p-value for log-likelihood 

              
g1 <- ggplot(pmc.alpha.df.reshape, aes(x=ll, fill=model)) + 
  geom_density(alpha=0.7) +
  ylab("Density") +
  xlab(TeX("Log-likelihood difference, $\\delta$")) +
  theme_bw() +
  geom_vline(xintercept=pmc.mttf$lr, linetype="dashed", color = "black") +

  
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggtitle(TeX("\\bar{T}_{d},  log_{10}")) +
  scale_fill_brewer(palette="Set1", name="Model", 
                    labels=c("Brownian motion", "Ornstein–Uhlenbeck"))
  




#ggsave(file="figs/PMC.png", g1, width=10,height=10, units='in', dpi=600)






