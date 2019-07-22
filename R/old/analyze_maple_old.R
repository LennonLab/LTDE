rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

set.seed(123456)

library('phytools')
library('ape')
library('ggplot2')
library('phylolm')
library('latex2exp')
library('ggpubr')


#source("http://www.phytools.org/utilities/v4.6/utilities.R")
#source("http://www.phytools.org/phyl.pca/v0.5/phyl.pca.R")
#require("corpcor")

df.met <- read.table("data/genomes/genomes_ncbi_maple.txt", 
                     header = TRUE, sep = "\t", row.names = 1)
df.met <- subset(df.met, select = -c( KBS0727))#, KBS0710, KBS0721))
# remove rows with all ones
df.met<- t(df.met[apply(df.met[,-1], 1, function(x) !all(x==1)),])

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

# remove columns with all zeros

#df <- df[,colSums(is.na(df))<nrow(df)]

df.met.no0 <-  df.met[, colSums(df.met != 0) > 0]
df.met.no0or1 <-  df.met.no0[, colSums(df.met.no0 != 1) < dim(df.met)[1]]

# exactly singular
ltde.met.ppca <- phyl.pca(ml.rooted.um, df.met.no0, method = "lambda", mode = "cov")


# run regression
df.species <- read.table("data/traits/traits_weibull.txt", 
                         header = TRUE, sep = "\t", 
                         row.names = 1, stringsAsFactors = FALSE)

rownames(df.species) <- df.species$Species
# remove bacillus
rownames(df.species) <- lapply(rownames(df.species), as.character)
rownames(ltde.met.ppca$S) <- lapply(rownames(ltde.met.ppca$S), as.character)

df.PCA.met.merge <- merge(df.species, ltde.met.ppca$S,by="row.names")

select.me <- c('Species', "alpha", "beta.log10","A","umax", "Lag", 'PC1', "PC2")
df.PCA.met.merge.subset <- df.PCA.met.merge[,select.me]
rownames(df.PCA.met.merge.subset) <- df.PCA.met.merge.subset$Species
met.explainvar1 <- round(diag(ltde.met.ppca$Eval)[1] / sum(diag(ltde.met.ppca$Eval)), 3) * 100
met.explainvar2 <- round(diag(ltde.met.ppca$Eval)[2] / sum(diag(ltde.met.ppca$Eval)), 3) * 100



traits.names <- c('Species', "alpha", "beta.log10","A","umax", "Lag", 'PC1', "PC2")
df.PCA.met.merge.subset <- df.PCA.met.merge[,select.me]


df.species.subset <- df.species[,c("alpha", "beta.log10","A","umax", "Lag")]
testtt <- phyl.cca(ml.rooted.um, df.met.no0, df.species.subset, fixed=F)



ml.rooted.um, df.met.no0, method = "lambda", mode = "cov")
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


phylo.alpha <- phylolm(alpha ~ PC1, data = df.PCA.met.merge.subset, 
                       ml.rooted.um.prunned, model = 'lambda', boot = 0)

phylo.beta <- phylolm(beta.log10 ~ PC1, data = df.PCA.met.merge.subset, 
                      ml.rooted.um.prunned, model = 'lambda', boot = 0)


phylo.umax <- phylolm(log10(umax) ~ PC1, data = df.PCA.met.merge.subset, 
                      ml.rooted.um.prunned, model = 'lambda', boot = 0)

phylo.yield <- phylolm(A ~ PC1, data = df.PCA.met.merge.subset, 
                       ml.rooted.um.prunned, model = 'lambda', boot = 0)

phylo.lag <- phylolm(Lag ~ PC1, data = df.PCA.met.merge.subset, 
                     ml.rooted.um.prunned, model = 'lambda', boot = 0)


pc1.beta.met.plot <- ggplot(data = df.PCA.met.merge.subset, aes(x = PC1, y = 10**beta.log10)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Mean scale paramater, $\\bar{\\lambda}$") ) +
  xlab(paste("MAPLE PC 1 (", met.explainvar1, "%)", sep = "")) + 
  #scale_y_continuous(limits = c(0, 1)) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=12), 
        axis.title.y = element_text(color="black", size=12), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


pc1.alpha.met.plot <- ggplot(data = df.PCA.met.merge.subset, aes(x = PC1, y = alpha)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Mean shape paramater, $\\bar{k}$")) +
  xlab(paste("MAPLE PC 1 (", met.explainvar1, "%)", sep = "")) + 
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=12), 
        axis.title.y = element_text(color="black", size=12), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


phylo.umax.x.line <- seq(-4, 2, length.out = 1000) 
phylo.umax.y.line <- phylo.umax$coefficients[1] + phylo.umax$coefficients[2] * phylo.umax.x.line
pc1.umax.met.plot <- ggplot(data = df.PCA.met.merge.subset, aes(x = PC1, y = umax)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Max. growth rate ($hours^{-1}$), $\\mu_{max}$") ) + 
  xlab(paste("MAPLE PC 1 (", met.explainvar1, "%)", sep = "")) + 
  #scale_y_continuous(limits = c(0, 1)) +
  geom_line(aes(y = y, x = x), size=0.75, data=data.frame(x=phylo.umax.x.line, y=10**phylo.umax.y.line)) +
  
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=12), 
        axis.title.y = element_text(color="black", size=12), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



pc1.yield.met.plot <- ggplot(data = df.PCA.met.merge.subset, aes(x = PC1, y = A)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Yield, $OD_{600}$")) +
  xlab(paste("MAPLE PC 1 (", met.explainvar1, "%)", sep = "")) + 
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=12), 
        axis.title.y = element_text(color="black", size=12), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


pc1.lag.met.plot <- ggplot(data = df.PCA.met.merge.subset, aes(x = PC1, y = Lag)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Lag time (hours)")) +
  xlab(paste("MAPLE PC 1 (", met.explainvar1, "%)", sep = "")) + 
  #scale_y_continuous(limits = c(0, 1)) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=12), 
        axis.title.y = element_text(color="black", size=12), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())






g.2 <- ggarrange(pc1.umax.met.plot, pc1.yield.met.plot, pc1.lag.met.plot,pc1.beta.met.plot, pc1.alpha.met.plot,
               # First row with scatter plot
               # Second row with box and dot plots
               ncol = 3, nrow = 2,
               labels = "auto")#, label.y = c(1, 0.5, 0.25)    


ggsave(file="figs/pPCR.png", g.2, width=9,height=6, units='in', dpi=600)


