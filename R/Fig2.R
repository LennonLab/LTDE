#### trait regression
rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('phylolm')

df.species <- read.table("data/demography/weibull_results_clean_species.csv", 
                         header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
rownames(df.species) <- df.species$Species
# remove bacillus
df.species.no_812<-df.species[!(df.species$Species=="KBS0812"),]

df.iRep <- read.table("data/iRep_clean.txt", 
                   header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
df.iRep.species <- aggregate(df.iRep[, c('iRep')], list(df.iRep$strain), mean)
colnames(df.iRep.species) <- c("Species","iRep")
rownames(df.iRep.species) <- df.iRep.species$Species
df.iRep.merge <- merge(df.species.no_812, df.iRep.species, by=0)
rownames(df.iRep.merge) <- df.iRep.merge$Species.x

# Load ML tree
ml.tree <- read.tree("data/tree/RAxML_bipartitionsBranchLabels.ltde_seqs")
# Define the outgroup
outgroup <- match("NC_005042.1_353331-354795", ml.tree$tip.label)
# Create a rooted tree {ape} using the outgroup
ml.rooted <- root(ml.tree, outgroup, resolve.root = TRUE)
# Keep rooted but drop outgroup branch
ml.rooted <- drop.tip(ml.rooted, c("NC_005042.1_353331-354795"))
is.ultrametric(ml.rooted)
ml.rooted.um  <- chronos(ml.rooted)
is.ultrametric(ml.rooted.um)
# pmc is having trouble converting the chronos object to a phy type object 
# just save the tree and re-load it
write.tree(ml.rooted.um, file = "data/tree/test.txt")
re.ml.rooted.um <- read.tree("data/tree/test.txt")
re.ml.rooted.um.prunned <- drop.tip(re.ml.rooted.um, 
                                    re.ml.rooted.um$tip.label[na.omit(match(c('KBS0812'),
                                                                            re.ml.rooted.um$tip.label))])

re.ml.rooted.um.prunned<-drop.tip(re.ml.rooted.um.prunned,re.ml.rooted.um.prunned$tip.label[-match(df.iRep.merge$Row.names, re.ml.rooted.um.prunned$tip.label)])

# Run a phylogeny-corrected regression with no bootstrap replicates
fit.phy.mttf <- phylolm(iRep ~ mttf.mean, data = df.iRep.merge, 
                        re.ml.rooted.um.prunned, model = 'lambda', lower.bound = 0, upper.bound = 1, boot = 0)
fit.phy.beta <- phylolm(iRep ~ beta.mean, data = df.iRep.merge, 
                         re.ml.rooted.um.prunned, model = 'lambda', lower.bound = 0, upper.bound = 1, boot = 0)
fit.phy.alpha <- phylolm(iRep ~ alpha.mean, data = df.iRep.merge, 
                   re.ml.rooted.um.prunned, model = 'lambda', lower.bound = 0, upper.bound = 1, boot = 0)


phy.mttf.plot <- ggplot(data = df.iRep.merge, aes(x = log10(beta.mean), y = alpha.mean)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  xlab(TeX("$\\bar{\\beta}$,  $\\log_{10}$") ) + 
  ylab(TeX("$\\bar{\\alpha}$")) +
  stat_function(fun = function(x) fit.phy$coefficients[1] + fit.phy$coefficients[2] * x) + 
  # how to plot CIs for this?
  #stat_function(fun = function(x) CI.2.5.inter + CI.2.5.slope  * log(x) + 
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()

# make plot
# Don't need slopes, since it's non-significant
CI.2.5.inter <- fit.phy.mttf$bootconfint95[1,1]
CI.97.5.inter <- fit.phy.mttf$bootconfint95[2,1]
CI.2.5.slope <- fit.phy.mttf$bootconfint95[1,2]
CI.97.5.slope <- fit.phy.mttf$bootconfint95[2,2]
