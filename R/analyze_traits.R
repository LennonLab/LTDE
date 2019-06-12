#### trait regression
rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('phylolm')
library('ggplot2')

df.species <- read.table("data/demography/weibull_results_clean_species.csv", 
                         header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
rownames(df.species) <- df.species$Species
traits <-  read.table("data/traits/traits.txt", 
                      header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

traits.merge <- merge(df.species, traits, by="row.names")
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

traits.merge$min.T.birth.log10 <- log10(1/traits.merge$umax)
#KBS0801 does not have a yield measurement
traits.merge.modelSel <- traits.merge[!(row.names(traits.merge) %in% c('KBS0801')), ]
ml.rooted.um.modelSel <- drop.tip(ml.rooted.um, c("KBS0801"))

fit.trait.mttf.select <- phylostep(mttf.log10 ~ A + min.T.birth.log10 + Lag, starting.formula=NULL, data= traits.merge.modelSel, phy=ml.rooted.um.modelSel, model = 'lambda', lower.bound = 0, upper.bound = 1)
fit.trait.alpha.select <- phylostep(alpha ~ A + min.T.birth.log10 + Lag, starting.formula=NULL, data= traits.merge.modelSel, phy=ml.rooted.um.modelSel, model = 'lambda', lower.bound = 0, upper.bound = 1)

summary(fit.trait.mttf.select)
summary(fit.trait.alpha.select)

fit.trait.alpha.phylo <- phylolm(alpha ~ A + Lag, data = traits.merge.modelSel, 
                                 ml.rooted.um.modelSel, model = 'BM', lower.bound = 0, upper.bound = 1, boot = 100)
summary(fit.trait.alpha.phylo)

phy.mttf.umax.plot <- ggplot(data = traits.merge.modelSel, aes(x = min.T.birth.log10, y = 10** mttf.log10)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Mean time to death, $\\bar{T}_{d}$ (days)") ) + 
  xlab(TeX("Min. time to birth $\\T_{b, min}$ (hours)") ) + 
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme(axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())#,
        #axis.text.x=element_blank(),
        #axis.title.x=element_blank())


phy.mttf.yield.plot <- ggplot(data = traits.merge, aes(x = A, y = 10**mttf.log10)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Mean time to death, $\\bar{T}_{d}$ (days)") ) + 
  xlab(TeX("Yield, $OD_{600}$")) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())#,
        #axis.text.x=element_blank(),
        #axis.title.x=element_blank())



phy.mttf.lag.plot <- ggplot(data = traits.merge, aes(x = Lag, y = 10**mttf.log10)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Mean time to death, $\\bar{T}_{d}$ (days)") ) + 
  xlab(TeX("Lag time (hours)")) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



phy.alpha.umax.plot <- ggplot(data = traits.merge, aes(x = umax, y = alpha)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Mean shape paramater, $\\bar{k}$")) +
  xlab(TeX("Min. time to birth $\\T_{b, min}$ (hours)") ) + 
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


phy.alpha.yield.plot <- ggplot(data = traits.merge, aes(x = A, y = alpha)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Mean shape paramater, $\\bar{k}$")) +
  xlab(TeX("Yield, $OD_{600}$")) +
  #scale_y_continuous(limits = c(0, 1)) +
  stat_function(fun = function(x) fit.trait.alpha.phylo$coefficients[1] + fit.trait.alpha.phylo$coefficients[2] * x) + 
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


phy.alpha.lag.plot <- ggplot(data = traits.merge, aes(x = Lag, y = alpha)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Mean shape paramater, $\\bar{k}$")) +
  xlab(TeX("Lag time (hours)")) +
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

