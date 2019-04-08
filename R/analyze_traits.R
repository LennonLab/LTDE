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
# and KBS0801, since there's no yield parameter
df.species.no_812<-df.species[!(df.species$Species=="KBS0812"),]
df.species.no_812<-df.species.no_812[!(df.species.no_812$Species=="KBS0801"),]

df.species.no_812<-df.species.no_812[!(df.species.no_812$Species=="KBS0727"),]

traits <-  read.table("data/traits/traits.txt", 
                      header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

traits.merge <- merge(df.species.no_812, traits, by="row.names")
rownames(traits.merge) <- traits.merge$Species
# save table
traits.merge <- traits.merge[, !(colnames(traits.merge) %in% c("Row.names"))]

#df.iRep <- read.table("data/iRep_clean.txt", 
#                      header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
#df.iRep.species <- aggregate(df.iRep[, c('iRep')], list(df.iRep$strain), mean)
#colnames(df.iRep.species) <- c("Species","iRep")
#rownames(df.iRep.species) <- df.iRep.species$Species
#df.iRep.merge <- merge(traits.merge, df.iRep.species, by=0)
#rownames(df.iRep.merge) <- df.iRep.merge$Species.x

traits.irep.merge <- traits.merge[, !(colnames(traits.merge) %in% c("Species.x", "Row.names", "Species.y"))]
write.table(traits.irep.merge, "data/traits/traits_weibull.txt", sep="\t")

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
                                 ml.rooted.um$tip.label[na.omit(match(c('KBS0812', 'KBS0801'),
                                                                      ml.rooted.um$tip.label))])
ml.rooted.um.prunned<-drop.tip(ml.rooted.um.prunned, ml.rooted.um.prunned$tip.label[-match(traits.merge$Row.names, ml.rooted.um.prunned$tip.label)])

traits.merge$mttf.mean.log10 <- log10(traits.merge$mttf.mean)

fit.trait.mttf.select <- phylostep(mttf.mean.log10 ~ A + umax + Lag, starting.formula=NULL, data= traits.merge, phy=ml.rooted.um.prunned, model = 'lambda', lower.bound = 0, upper.bound = 1)
fit.trait.alpha.select <- phylostep(alpha.mean ~ A + umax + Lag, starting.formula=NULL, data= traits.merge, phy=ml.rooted.um.prunned, model = 'lambda', lower.bound = 0, upper.bound = 1)
#fit.trait.beta.select <- phylostep(beta.mean ~ A + umax + Lag, starting.formula=NULL, data= traits.merge, phy=ml.rooted.um.prunned, model = 'lambda', lower.bound = 0, upper.bound = 1)

summary(fit.trait.mttf.select)
summary(fit.trait.alpha.select)

fit.trait.alpha.phylo <- phylolm(alpha.mean ~ A + Lag, data = traits.merge, 
                                 ml.rooted.um.prunned, model = 'lambda', lower.bound = 0, upper.bound = 1, boot = 100)
summary(fit.trait.alpha.phylo)

phy.mttf.umax.plot <- ggplot(data = traits.merge, aes(x = umax, y = mttf.mean)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("$\\bar{\\T_{death}$} (days),  $\\log_{10}$") ) +  
  xlab(TeX("$\\mu_{max}$")) +
  scale_y_log10() + 
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())#,
        #axis.text.x=element_blank(),
        #axis.title.x=element_blank())


phy.mttf.yield.plot <- ggplot(data = traits.merge, aes(x = A, y = mttf.mean)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("$\\bar{\\T_{death}$} (days),  $\\log_{10}$") ) + 
  xlab(TeX("Yield")) +
  scale_y_log10() + 
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())#,
        #axis.text.x=element_blank(),
        #axis.title.x=element_blank())



phy.mttf.lag.plot <- ggplot(data = traits.merge, aes(x = Lag, y = mttf.mean)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("$\\bar{\\T_{death}$} (days),  $\\log_{10}$") ) + 
  xlab(TeX("Lag Time")) +
  scale_y_log10() + 
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



phy.alpha.umax.plot <- ggplot(data = traits.merge, aes(x = umax, y = alpha.mean)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("$\\bar{\\alpha}$") ) + 
  xlab(TeX("$\\mu_{max}$")) +
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


phy.alpha.yield.plot <- ggplot(data = traits.merge, aes(x = A, y = alpha.mean)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("$\\bar{\\alpha}$") ) + 
  xlab(TeX("Yield")) +
  #scale_y_continuous(limits = c(0, 1)) +
  stat_function(fun = function(x) fit.trait.alpha.phylo$coefficients[1] + fit.trait.alpha.phylo$coefficients[2] * x) + 
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


phy.alpha.lag.plot <- ggplot(data = traits.merge, aes(x = Lag, y = alpha.mean)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("$\\bar{\\alpha}$") ) + 
  xlab(TeX("Lag Time")) +
  #scale_y_continuous(limits = c(0, 1)) +
  stat_function(fun = function(x) fit.trait.alpha.phylo$coefficients[1] + fit.trait.alpha.phylo$coefficients[3] * x) + 
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


g <- ggarrange(phy.mttf.umax.plot, phy.mttf.yield.plot, phy.mttf.lag.plot,
               phy.alpha.umax.plot, phy.alpha.yield.plot, phy.alpha.lag.plot,
               # First row with scatter plot
               # Second row with box and dot plots
               ncol = 3, nrow = 2,
               labels = "auto")#, label.y = c(1, 0.5, 0.25)    

ggsave(file="figs/traits.png", g, width=15,height=10, units='in', dpi=600)

