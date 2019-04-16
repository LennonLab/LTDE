rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('phytools')
library('ape')
library('ggplot2')
library('phylolm')


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
# run regression
df.species <- read.table("data/traits/traits_weibull.txt", 
                         header = TRUE, sep = "\t", 
                         row.names = 1, stringsAsFactors = FALSE)
#rownames(df.species) <- df.species$Species
# remove bacillus
#df.species.no_812<-df.species[!(df.species$Species=="KBS0812"),]
#df.species.no_812<-df.species.no_812[!(df.species.no_812$Species=="KBS0727"),]
rownames(df.species) <- lapply(rownames(df.species), as.character)
rownames(ltde.met.ppca$S) <- lapply(rownames(ltde.met.ppca$S), as.character)

df.PCA.met.merge <- merge(df.species, ltde.met.ppca$S,by="row.names")
df.species.2 <- cbind(df.species)
df.PCA.cog.merge <- merge(df.species.2, ltde.cog.ppca$S,by="row.names")

select.me <- c('Species', "alpha.mean", "mttf.mean","A","umax", "Lag", 'PC1', "PC2")
df.PCA.met.merge.subset <- df.PCA.met.merge[,select.me]
df.PCA.cog.merge.subset <- df.PCA.cog.merge[,select.me]

met.explainvar1 <- round(diag(ltde.met.ppca$Eval)[1] / sum(diag(ltde.met.ppca$Eval)), 3) * 100
met.explainvar2 <- round(diag(ltde.met.ppca$Eval)[2] / sum(diag(ltde.met.ppca$Eval)), 3) * 100
cog.explainvar1 <- round(diag(ltde.cog.ppca$Eval)[1] / sum(diag(ltde.cog.ppca$Eval)), 3) * 100
cog.explainvar2 <- round(diag(ltde.cog.ppca$Eval)[2] / sum(diag(ltde.cog.ppca$Eval)), 3) * 100



met.PCA <- ggplot(data = df.PCA.met.merge.subset, aes(x = PC1, y = PC2)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  xlab(paste("KEGG PC 1 (", met.explainvar1, "%)", sep = "")) + 
  ylab(paste("KEGG PC 2 (", met.explainvar2, "%)", sep = "")) +
  scale_x_continuous(limits = c(-4, 4)) +
  scale_y_continuous(limits = c(-4, 4)) +
  geom_vline(xintercept=0, linetype="dashed", color = "black") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
 

cog.PCA <- ggplot(data = df.PCA.cog.merge.subset, aes(x = PC1, y = PC2)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  xlab(paste("COG ", "PC 1 (", cog.explainvar1, "%)", sep = "")) + 
  ylab(paste("COG ", "PC 2 (", cog.explainvar2, "%)", sep = "")) +
  scale_x_continuous(limits = c(-3, 3)) +
  scale_y_continuous(limits = c(-3, 3)) +
  geom_vline(xintercept=0, linetype="dashed", color = "black") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



g <- ggarrange(cog.PCA, met.PCA,
               # First row with scatter plot
               # Second row with box and dot plots
               ncol = 2, nrow = 1,
               labels = "auto")#, label.y = c(1, 0.5, 0.25)    


ggsave(file="figs/pPCA.png", g, width=10,height=5, units='in', dpi=600)

#fit.phy <- phylolm(alpha  ~ log10(beta), data = df.species.no_812, 
#                   ml.rooted.um.prunned, model = 'OUrandomRoot', boot = 10)


summary(lm(log10(df.PCA.cog.merge.subset$mttf.mean) ~ df.PCA.cog.merge.subset$PC1))
summary(lm(df.PCA.cog.merge.subset$umax ~ df.PCA.cog.merge.subset$PC1))
#ml.rooted.um.prunned.prunned<-drop.tip(ml.rooted.um.prunned, ml.rooted.um.prunned$tip.label[-match(df.PCA.cog.merge.subset$Species, ml.rooted.um.prunned$tip.label)])



summary(lm(log10(df.PCA.met.merge.subset$mttf.mean) ~ df.PCA.met.merge.subset$PC1))
summary(lm(df.PCA.met.merge.subset$umax ~ df.PCA.met.merge.subset$PC1))

pc1.mttf.cog.plot <- ggplot(data = df.PCA.cog.merge.subset, aes(x = PC1, y = mttf.mean)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  #scale_y_log10() + 
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  ylab(TeX("$\\bar{T}_{d}$ (days)") ) + 
  xlab(paste("COG ", "PC 1 (", cog.explainvar1, "%)", sep = "")) + 
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


pc1.umax.cog.plot <- ggplot(data = df.PCA.cog.merge.subset, aes(x = PC1, y = umax)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  geom_smooth(method=lm, se=TRUE, color="black") +
  ylab(TeX("$\\mu_{max}$")) +
  xlab(paste("COG ", "PC 1 (", cog.explainvar1, "%)", sep = "")) + 
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())





pc1.mttf.met.plot <- ggplot(data = df.PCA.met.merge.subset, aes(x = PC1, y = mttf.mean)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  ylab(TeX("$\\bar{T}_{d}$ (days)") ) + 
  xlab(paste("KEGG PC 1 (", met.explainvar1, "%)", sep = "")) + 
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


pc1.umax.met.plot <- ggplot(data = df.PCA.met.merge.subset, aes(x = PC1, y = umax)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  geom_smooth(method=lm, se=TRUE, color="black") +
  ylab(TeX("$\\mu_{max}$")) +
  xlab(paste("KEGG PC 1 (", met.explainvar1, "%)", sep = "")) + 
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


g.2 <- ggarrange(pc1.mttf.cog.plot, pc1.umax.cog.plot, pc1.mttf.met.plot, pc1.umax.met.plot,
               # First row with scatter plot
               # Second row with box and dot plots
               ncol = 2, nrow = 2,
               labels = "auto")#, label.y = c(1, 0.5, 0.25)    


ggsave(file="figs/pPCR.png", g.2, width=10,height=10, units='in', dpi=600)

