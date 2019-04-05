rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('phylolm')
library('ape')

df.species <- read.table("data/demography/weibull_results_clean_species.csv", 
                         header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
rownames(df.species) <- df.species$Species
# remove bacillus
df.species.no_812<-df.species[!(df.species$Species=="KBS0812"),]

traits <-  read.table("data/traits/traits.txt", 
                      header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

traits.merge <- merge(df.species.no_812, traits, by=0)
rownames(traits.merge) <- traits.merge$Species


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

re.ml.rooted.um.prunned<-drop.tip(re.ml.rooted.um.prunned, re.ml.rooted.um.prunned$tip.label[-match(rownames(traits.merge), re.ml.rooted.um.prunned$tip.label)])


plot(log10(traits.merge$alpha.mean), traits.merge$umax)
plot(log10(traits.merge$alpha.mean), traits.merge$A)
plot(log10(traits.merge$alpha.mean), traits.merge$Lag)

fit.phy <- phylolm(A  ~ alpha.mean, data = traits.merge, 
                   re.ml.rooted.um.prunned, model = 'lambda', boot = 0, lower.bound = 0, upper.bound = 1)

summary(fit.phy)

# find traits that explain mean time to failure

