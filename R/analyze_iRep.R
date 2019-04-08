#### trait regression
rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('phylolm')
library('ggplot2')

df.species <- read.table("data/demography/weibull_results_clean_species.csv", 
                         header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
rownames(df.species) <- df.species$Species
# remove bacillus
df.species.no_812<-df.species[!(df.species$Species=="KBS0812"),]
df.species.no_812<-df.species.no_812[!(df.species.no_812$Species=="KBS0727"),]

df.iRep <- read.table("data/iRep_clean.txt", 
                   header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
df.iRep.species <- aggregate(df.iRep[, c('iRep')], list(df.iRep$strain), mean)
colnames(df.iRep.species) <- c("Species","iRep")
rownames(df.iRep.species) <- df.iRep.species$Species
df.iRep.merge <- merge(df.species.no_812, df.iRep.species, by=0)
rownames(df.iRep.merge) <- df.iRep.merge$Species.x


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
ml.rooted.um.prunned<-drop.tip(ml.rooted.um.prunned, ml.rooted.um.prunned$tip.label[-match(df.iRep.merge$Row.names, ml.rooted.um.prunned$tip.label)])

# Run a phylogeny-corrected regression with no bootstrap replicates
fit.phy.mttf <- phylolm(iRep ~ mttf.mean, data = df.iRep.merge, 
                        ml.rooted.um.prunned, model = 'lambda', lower.bound = 0, upper.bound = 1, boot = 0)
fit.phy.beta <- phylolm(iRep ~ beta.mean, data = df.iRep.merge, 
                        ml.rooted.um.prunned, model = 'lambda', lower.bound = 0, upper.bound = 1, boot = 0)
fit.phy.alpha <- phylolm(iRep ~ alpha.mean, data = df.iRep.merge, 
                         ml.rooted.um.prunned, model = 'lambda', lower.bound = 0, upper.bound = 1, boot = 0)


phy.beta.plot <- ggplot(data = df.iRep.merge, aes(x = log10(beta.mean), y = iRep)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  xlab(TeX("$\\bar{\\beta}$,  $\\log_{10}$") ) + 
  ylab(TeX("$\\bar{\\iRep}$")) +
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw()

phy.alpha.plot <- ggplot(data = df.iRep.merge, aes(x = alpha.mean, y = iRep)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  xlab(TeX("$\\bar{\\alpha}$") ) + 
  ylab(TeX("$\\bar{\\iRep}$")) +
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw()

phy.mttf.plot <- ggplot(data = df.iRep.merge, aes(x = mttf.mean, y = iRep)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  xlab(TeX("$\\T_{death}$ (days),  $\\log_{10}$") ) + 
  ylab(TeX("$\\bar{\\iRep}$")) +
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw()


phy.beta.plot <- phy.beta.plot + theme(axis.title.x = element_text(color="black", size=14), 
                                         axis.title.y = element_text(color="black", size=14))

phy.alpha.plot <- phy.alpha.plot + theme(axis.title.x = element_text(color="black", size=14), 
                                         axis.title.y = element_text(color="black", size=14))

phy.mttf.plot <- phy.mttf.plot + theme(axis.title.x = element_text(color="black", size=14), 
                                         axis.title.y = element_text(color="black", size=14))


g <- ggarrange(phy.beta.plot, phy.alpha.plot, phy.mttf.plot,                                              # First row with scatter plot
               # Second row with box and dot plots
               ncol = 1, nrow = 3,
               labels = "auto")#, label.y = c(1, 0.5, 0.25)                                     # Labels of the scatter plot
#) 


ggsave(file="figs/Fig2.png", g,width=5,height=15, units='in', dpi=600)

