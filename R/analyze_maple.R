rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

set.seed(123456)

library('ggplot2')
library('latex2exp')
library('ggpubr')
library('vegan')
library('ggrepel')

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

#maple.dj <- vegdist(df.met.no0or1, method="jaccard", binary=T)
#doubs.dbrda <- dbrda(maple.dj ~., as.data.frame(df.species.subset))
#drbda.anova <- anova(doubs.dbrda, by="axis")
#ordiplot(doubs.dbrda)

# rda
m1.rda <- vegan::rda(df.met.no0or1 ~ ., as.data.frame(df.species.subset))
m0.rda <- vegan::rda(df.met.no0or1 ~ 1, as.data.frame(df.species.subset))
m.rda <- step(m0.rda, scope=formula(m1.rda), test="perm")
# significance of final model
anova.cca(m.rda, by="term")
ordiplot(m1.rda)



# rda w/out bacillus
df.met.no0or1.noKBS0812 <- df.met.no0or1[!rownames(df.met.no0or1) %in% c('KBS0812'), ]
df.species.subset.noKBS0812 <- df.species.subset[!rownames(df.species.subset) %in% c('KBS0812'), ]

m1.rda.noKBS0812 <- vegan::rda(df.met.no0or1.noKBS0812 ~ ., as.data.frame(df.species.subset.noKBS0812))
m0.rda.noKBS0812 <- vegan::rda(df.met.no0or1.noKBS0812 ~ 1, as.data.frame(df.species.subset.noKBS0812))
m.rda.noKBS0812 <- step(m0.rda.noKBS0812, scope=formula(m1.rda.noKBS0812), test="perm")
# reduced model
m.rda.reduced.noKBS0812 <- vegan::rda(df.met.no0or1.noKBS0812 ~ umax, as.data.frame(df.species.subset.noKBS0812))





# make plot
# the code for this plot was adapted from the GitHub repo HJA-streams
# under a GNU General Public License v3.0
# https://github.com/nwisnoski/HJA-streams
rda.plot <- cbind.data.frame(scores(m1.rda.noKBS0812)$sites)
rda.var1 <- round(eigenvals(m1.rda.noKBS0812)[1] / sum(eigenvals(m1.rda.noKBS0812)) * 100, 2)
rda.var2 <- round(eigenvals(m1.rda.noKBS0812)[2] / sum(eigenvals(m1.rda.noKBS0812)) * 100, 2)
rda.vecs <- as.data.frame(m1.rda.noKBS0812$CCA$biplot)
rda.vecs$predictor <- c("alpha", "beta.log10", "A", "umax", "Lag")
rda.vecs$origin <- 0
scale.arrows = 4

ggplot(data = rda.plot, aes(RDA1, RDA2)) +
  geom_hline(aes(yintercept = 0), color = 'black', linetype = "dashed", alpha = 1, size = 0.25) +
  geom_vline(aes(xintercept = 0),color = 'black',  linetype = "dashed", alpha = 1, size = 0.25) +
  #geom_point(aes(color = 'blue'), size = 2, alpha = 0.8) +
  geom_point(color='blue', alpha = 0.8, size=5)+
  #geom_point(data = subset(rda.plot, group == "Sediment"), shape = 1, color = "black", size = 2) +
  #geom_point(data = subset(rda.plot, group == "Water"), shape = 2, color = "black", size = 2) +
  # geom_path(data = df_ell, 
  #           aes(x = Dim1, y = Dim2, color = group),
  #           size = .8, alpha = 1, linetype = 2) +
  labs(x = paste0("RDA1 (", rda.var1, "%)"),
       y = paste0("RDA2 (", rda.var2, "%)"),
       color = "Habitat", shape = "Habitat") +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_segment(data = rda.vecs, size = .5,
               aes(x = origin, y = origin, 
                   xend = scale.arrows*RDA1, 
                   yend = scale.arrows*RDA2),
               alpha = 1, color = "black",
               arrow = arrow(angle = 20,
                             length = unit(.1, "inches"),
                             type = "open")) +
  annotate("text", x=-2.9, y=1.5, label=TeX("$\\mu_{max}^{*}$"), size = 4) +
  annotate("text", x=1.8, y=1.6, label=TeX("$Yield$"), size = 4) +
  annotate("text", x=-1, y=-2.4, label=TeX("$\\log_{10} (\\lambda) $"), size = 4) +
  annotate("text", x=0.6, y=-2.6, label=TeX("$k$"), size = 4) +
  annotate("text", x=1.5, y=-3.6, label="Lag time", size = 4) +
  ggsave("figs/RDA.png", width = 84, height = 84, units = "mm", dpi = 500)














# cca
#
maple.cca <- vegan::cca(df.met.no0or1 ~ as.matrix(df.species.subset))
cca.anova <- anova(maple.cca, by="axis")
ordiplot(m1)

m1 <- vegan::cca(df.met.no0or1 ~ ., as.data.frame(df.species.subset))
m0 <- vegan::cca(df.met.no0or1 ~ 1, as.data.frame(df.species.subset))
m <- step(m0, scope=formula(m1), test="perm")
# significance of final model
anova.cca(m, by="term")
ordiplot(m)



# cca w/out bacillus
df.met.no0or1.noKBS0812 <- df.met.no0or1[!rownames(df.met.no0or1) %in% c('KBS0812'), ]
df.species.subset.noKBS0812 <- df.species.subset[!rownames(df.species.subset) %in% c('KBS0812'), ]

maple.cca.noKBS0812 <- vegan::cca(df.met.no0or1.noKBS0812 ~ as.matrix(df.species.subset.noKBS0812))
cca.anova.noKBS0812 <- anova(maple.cca.noKBS0812, by="axis")
ordiplot(cca.anova.noKBS0812)

m1.noKBS0812 <- vegan::cca(df.met.no0or1.noKBS0812 ~ ., as.data.frame(df.species.subset.noKBS0812))
m0.noKBS0812 <- vegan::cca(df.met.no0or1.noKBS0812 ~ 1, as.data.frame(df.species.subset.noKBS0812))
m.noKBS0812 <- step(m0.noKBS0812, scope=formula(m1.noKBS0812), test="perm")
# significance of final model
anova.cca(m.noKBS0812, by="term")
#ordiplot(m.noKBS0812)




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
sig.corr <- append(pos.sig.corr, neg.sig.corr)
# make table and make bi plot with labels describing which species is where



cca.res<-summary(m)
cca.taxa <-data.frame(cca.res$sites)
ord_df<-data.frame(CCA1=cca.taxa$CCA1,CCA2=cca.taxa$CCA2)
#ord_df$Year <- factor(ord_df$Year)
#ord_df$Site <- factor(ord_df$Site)
#exp<-cca.res$concont
#exp<-data.frame(exp$importance)
#cca.species<-data.frame(cca.res$species)
#cca.species<-data.frame(Cca1=cca.species$CCA1,Cca2=cca.species$CCA2,species=rownames(cca.species))
#cca.benthos<-data.frame(cca.res$biplot)
#cca.benthos<-data.frame(cca1=cca.benthos$CCA1,cca2=cca.benthos$CCA2, Species = rownames(cca.benthos))








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


