rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('phytools')
library('ape')

#source("http://www.phytools.org/utilities/v4.6/utilities.R")
#source("http://www.phytools.org/phyl.pca/v0.5/phyl.pca.R")
#require("corpcor")

df.met <- read.table("data/metab_paths/module_by_taxon.txt", 
                     header = TRUE, sep = "\t", row.names = 1)
df.met <- subset(df.met, select = -c(KBS0812, KBS0727))#, KBS0710, KBS0721))
# remove rows with all ones
df.met<- t(df.met[apply(df.met[,-1], 1, function(x) !all(x==1)),])

df.cog <- read.table("data/COGs/cog_by_genome.txt", 
                     header = TRUE, sep = "\t", row.names = 1)
df.cog <- subset(t(df.cog), select = -c(KBS0812, KBS0727))#, KBS0710, KBS0721))
# remove rows with all ones
df.cog<- t(df.cog[apply(df.cog[,-1], 1, function(x) !all(x==1)),])

# Load ML tree
ml.tree <- read.tree("data/tree/RAxML_bipartitionsBranchLabels.ltde")
# Define the outgroup
outgroup <- match("NC_005042.1.353331-354795", ml.tree$tip.label)
# Create a rooted tree {ape} using the outgroup
ml.rooted <- root(ml.tree, outgroup, resolve.root = TRUE)
# Keep rooted but drop outgroup branch
ml.rooted <- drop.tip(ml.rooted, c("NC_005042.1.353331-354795", "KBS0812"))
is.ultrametric(ml.rooted)
ml.rooted.um  <- chronos(ml.rooted)
is.ultrametric(ml.rooted.um)
ml.rooted.um.prunned <- drop.tip(ml.rooted.um, 
                                 ml.rooted.um$tip.label[na.omit(match(c('KBS0812'),
                                                                      ml.rooted.um$tip.label))])

ltde.met.ppca <- phyl.pca(ml.rooted, df.met, method = "BM")
ltde.cog.ppca <- phyl.pca(ml.rooted, df.cog, method = "BM")

plot(ltde.met.ppca$S[,1], 
     ltde.met.ppca$S[,2], ylab = "", xlab = "")

plot(ltde.cog.ppca$S[,1], 
     ltde.cog.ppca$S[,2], ylab = "", xlab = "")

#diag(ltde.ppca$Eval)[1] / sum(diag(ltde.ppca$Eval))


# run regression
df.species <- read.table("data/demography/weibull_results_clean_species.csv", 
                         header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
rownames(df.species) <- df.species$Species
# remove bacillus
df.species.no_812<-df.species[!(df.species$Species=="KBS0812"),]
df.species.no_812<-df.species.no_812[!(df.species.no_812$Species=="KBS0727"),]
rownames(df.species.no_812) <- lapply(rownames(df.species.no_812), as.character)
rownames(ltde.met.ppca$S) <- lapply(rownames(ltde.met.ppca$S), as.character)

df.PCA.met.merge <- merge(df.species.no_812, ltde.met.ppca$S,by="row.names")
df.species.no_812.2 <- cbind(df.species.no_812)
df.PCA.cog.merge <- merge(df.species.no_812.2, ltde.cog.ppca$S,by="row.names")

plot(df.PCA.met.merge$PC1, df.PCA.met.merge$alpha.mean)
summary(lm( df.PCA.met.merge$alpha.mean ~ df.PCA.met.merge$PC1))



plot(df.PCA.cog.merge$PC1, df.PCA.cog.merge$alpha.mean)
summary(lm( df.PCA.cog.merge$alpha.mean ~ df.PCA.cog.merge$PC1))


# work on figures
# 2 x 2 for PCAs
# 2 x 2 for PC regression


