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

traits.merge <- merge(df.species.no_812, traits, by=0)
rownames(traits.merge) <- traits.merge$Species


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

fit.trait.mttf.select <- phylostep(mttf.mean ~ A + umax + Lag + umax:Lag, starting.formula=NULL, data= traits.merge, phy=ml.rooted.um.prunned, model = 'lambda', lower.bound = 0, upper.bound = 1)
fit.trait.alpha.select <- phylostep(alpha.mean ~ A + umax + Lag + umax:Lag, starting.formula=NULL, data= traits.merge, phy=ml.rooted.um.prunned, model = 'lambda', lower.bound = 0, upper.bound = 1)
fit.trait.beta.select <- phylostep(beta.mean ~ A + umax + Lag + umax:Lag, starting.formula=NULL, data= traits.merge, phy=ml.rooted.um.prunned, model = 'lambda', lower.bound = 0, upper.bound = 1)


summary(fit.trait.mttf.select)
summary(fit.trait.alpha.select)
summary(fit.trait.beta.select)

summary(lm(traits.merge$Lag ~traits.merge$umax))
plot(traits.merge$umax, traits.merge$Lag)
