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
df<-df[!(df$Species=="KBS0727" | df$Species=="KBS0812"),]
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


mttf <- df$mttf.log10
names(mttf) <- df$Species
alpha <- df$alpha
names(alpha) <- df$Species

test.mttf <- pmc(ml.rooted.um.prunned, mttf, "lambda", "OU", nboot = 100)
test.alpha <- pmc(ml.rooted.um.prunned, alpha, "lambda", "OU", nboot = 100)

length(test.mttf$null[test.mttf$null > test.mttf$lr]) /100
length(test.alpha$null[test.alpha$null > test.alpha$lr]) /100


dists <- data.frame(null = test.mttf$null, test = test.mttf$test)

dists %>% 
  gather(dist, value) %>%
  ggplot(aes(value, fill = dist)) + 
  geom_density(alpha = 0.5) + 
  geom_vline(xintercept = test.mttf$lr)

length(dists$null[dists$null > test.mttf$lr]) /100



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
df.summary$llr <- as.numeric(as.character(df.summary$llr))
write.csv(df.summary, file = "data/pmc/pmc_summary.csv")




# get p-value for log-likelihood 

BM_OU.mttf <- ggplot(mttf.df[mttf.df$test == "BM_OU", c("llr", "model")], aes(x=llr, fill=model)) + 
  geom_density(alpha=0.7) +
  ylab("Density") +
  xlab(TeX("Log-likelihood difference, $\\delta$")) +
  theme_bw() +
  geom_vline(xintercept= df.summary[df.summary$parameter == "mttf" & df.summary$test == "BM_OU" ,]$llr, linetype="dashed", color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggtitle(TeX("\\bar{T}_{d},  log_{10}")) +
  scale_fill_brewer(palette="Set1", name="Model", 
                    labels=c("Brownian motion", "Ornstein–Uhlenbeck"))
  


BM_GBM.mttf <- ggplot(mttf.df[mttf.df$test == "BM_GBM", c("llr", "model")], aes(x=llr, fill=model)) + 
  geom_density(alpha=0.7) +
  ylab("Density") +
  xlab(TeX("Log-likelihood difference, $\\delta$")) +
  theme_bw() +
  geom_vline(xintercept= df.summary[df.summary$parameter == "mttf" & df.summary$test == "BM_GBM" ,]$llr, linetype="dashed", color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggtitle(TeX("\\bar{T}_{d},  log_{10}")) +
  scale_fill_brewer(palette="Set1", name="Model", 
                    labels=c("Brownian motion", "Geometric Brownian motion"))


ggsave(file="figs/test_pmc1.png", BM_OU.mttf, units='in', dpi=600)
ggsave(file="figs/test_pmc2.png", BM_GBM.mttf, units='in', dpi=600)


g <- ggarrange(BM_OU.mttf, BM_GBM.mttf,                                              # First row with scatter plot
               nrow = 1, ncol =2,
               labels = "auto") 


ggsave(file="figs/test_pmc.png", g, units='in', dpi=600)



