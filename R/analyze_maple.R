rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

set.seed(123456)

library('ggplot2')
library('latex2exp')
library('ggpubr')
library('vegan')

#source("http://www.phytools.org/utilities/v4.6/utilities.R")
#source("http://www.phytools.org/phyl.pca/v0.5/phyl.pca.R")
#require("corpcor")

df.met <- read.table("data/genomes/genomes_ncbi_maple.txt", 
                     header = TRUE, sep = "\t", row.names = 1)
df.met <- subset(df.met, select = -c( KBS0727))#, KBS0710, KBS0721))
# remove rows with all ones
df.met<- t(df.met[apply(df.met[,-1], 1, function(x) !all(x==1)),])

df.met.no0 <-  df.met[, colSums(df.met != 0) > 0]
df.met.no0or1 <-  df.met.no0[, colSums(df.met.no0 != 1) < dim(df.met)[1]]

# run regression
df.species <- read.table("data/traits/traits_weibull.txt", 
                         header = TRUE, sep = "\t", 
                         row.names = 1, stringsAsFactors = FALSE)

rownames(df.species) <- df.species$Species
rownames(df.species) <- lapply(rownames(df.species), as.character)
df.species.subset <- df.species[,c("alpha", "beta.log10","A","umax", "Lag")]

maple.dj <- vegdist(df.met.no0or1, method="jaccard", binary=T)

doubs.dbrda <- dbrda(maple.dj ~., as.data.frame(df.species.subset))
drbda.anova <- anova(doubs.dbrda, by="axis")
ordiplot(doubs.dbrda)

#testtt <- phyl.cca(ml.rooted.um, df.met.no0, df.species.subset, fixed=F)

# cca
df.met.no0or1 <- df.met.no0or1[!rownames(df.met.no0or1) %in% c('KBS0812'), ]
df.species.subset <- df.species.subset[!rownames(df.species.subset) %in% c('KBS0812'), ]

maple.cca <- vegan::cca(df.met.no0or1 ~ as.matrix(df.species.subset))
cca.anova <- anova(maple.cca, by="axis")
ordiplot(maple.cca)

m1 <- vegan::cca(df.met.no0or1 ~ ., as.data.frame(df.species.subset))
m0 <- vegan::cca(df.met.no0or1 ~ 1, as.data.frame(df.species.subset))
m <- step(m0, scope=formula(m1), test="perm")
# significance of final model
anova.cca(m, by="term")

# identify significant pathways along CCA axis
# relative abundance of pathways
mapleREL <- df.met.no0or1
for(i in 1:nrow(df.met.no0or1)){
  mapleREL[i, ] = df.met.no0or1[i, ] / sum(df.met.no0or1[i, ])
} 

# permutation test for metabolic pathway composition
pathway.perm <- envfit(m, mapleREL, perm=999)
# correlation cutoff of 0.7
sig.pathways <- pathway.perm$vectors$pvals[pathway.perm$vectors$pvals < 0.05]
pathway.corr <- pathway.perm$vectors$arrows[,1]
corr.pos.pathways <- pathway.corr[pathway.corr > 0.7]
corr.neg.pathways <- pathway.corr[pathway.corr < -0.7]

pos.sig.corr <- intersect( names(corr.pos.pathways), names(sig.pathways))
neg.sig.corr <- intersect( names(corr.neg.pathways), names(sig.pathways))



# make table and make bi plot with labels describing which species is where



rownames(df.PCA.met.merge.subset) <- df.PCA.met.merge.subset$Species
met.explainvar1 <- round(diag(ltde.met.ppca$Eval)[1] / sum(diag(ltde.met.ppca$Eval)), 3) * 100
met.explainvar2 <- round(diag(ltde.met.ppca$Eval)[2] / sum(diag(ltde.met.ppca$Eval)), 3) * 100






# plot PCA
met.PCA <- ggplot(data = df.PCA.met.merge.subset, aes(x = PC1, y = PC2)) +
  geom_point(color='blue', alpha = 0.6, size=10) +
  xlab(paste("MAPLE PC 1 (", met.explainvar1, "%)", sep = "")) + 
  ylab(paste("MAPLE PC 2 (", met.explainvar2, "%)", sep = "")) +
  scale_x_continuous(limits = c(-4, 4)) +
  scale_y_continuous(limits = c(-4, 4)) +
  geom_vline(xintercept=0, linetype="dashed", color = "black") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=22), 
        axis.title.y = element_text(color="black", size=22), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave(file="figs/pPCA.png", met.PCA, width=10,height=10, units='in', dpi=600)



ml.rooted.um.prunned <- drop.tip(ml.rooted.um, ml.rooted.um$tip.label[-match(df.PCA.met.merge.subset$Species, ml.rooted.um$tip.label)])











#g.2 <- ggarrange(pc1.umax.met.plot, pc1.yield.met.plot, pc1.lag.met.plot,pc1.beta.met.plot, pc1.alpha.met.plot,
#               # First row with scatter plot
#               # Second row with box and dot plots
#               ncol = 3, nrow = 2,
#               labels = "auto")#, label.y = c(1, 0.5, 0.25)    


#ggsave(file="figs/pPCR.png", g.2, width=9,height=6, units='in', dpi=600)


