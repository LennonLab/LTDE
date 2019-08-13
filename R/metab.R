rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('ggplot2')
library('phylolm')
library("ape")
library('latex2exp')


df.metab <- read.table("data/metab_stain/staining.all.new.txt", sep="\t", 
                                 header=TRUE, stringsAsFactors = FALSE)

df.metab.noNA <- df.metab[complete.cases(df.metab), ]
df.metab.noNA.anc <- df.metab.noNA[df.metab.noNA$hist == "anc", ] 
df.metab.noNA.anc.mean <- aggregate(df.metab.noNA.anc[,4:6], list(df.metab.noNA.anc$strain), mean )
colnames(df.metab.noNA.anc.mean) <- c("Species", "Active.anc", "Dead.anc", "Dormant.anc")
# derived
df.metab.noNA.der <- df.metab.noNA[df.metab.noNA$hist == "der", ] 
df.metab.noNA.der.mean <- aggregate(df.metab.noNA.der[,4:6], list(df.metab.noNA.der$strain), mean )
colnames(df.metab.noNA.der.mean) <- c("Species", "Active.der", "Dead.der", "Dormant.der")
# merge
df.metab.mean.merge <- merge(df.metab.noNA.anc.mean, df.metab.noNA.der.mean, by=c("Species"))
df.metab.mean.merge$dorm.diff <- df.metab.mean.merge$Dormant.der - df.metab.mean.merge$Dormant.anc
df.metab.mean.merge$dead.diff <- df.metab.mean.merge$Dead.der - df.metab.mean.merge$Dead.anc

rownames(df.metab.mean.merge) <- df.metab.mean.merge$Species

# Death curve
df.species <- read.table("data/demography/weibull_results_clean_species.csv", 
                         header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
rownames(df.species) <- df.species$Species
df.species<-df.species[!(df.species$Species=="KBS0727") ,]
#df.species<-df.species[!(df.species$Species=="KBS0812") ,]

df.species.merge <- merge(df.metab.mean.merge, df.species, by="row.names",all.x=TRUE)
rownames(df.species.merge) <- df.species.merge$Row.names
df.species.merge <- df.species.merge[!(df.species.merge$Species.x=="KBS0711W" ),]
#x <- x[!(x$strain=="KBS0812" ),]


plot(df.species.merge$Dead.der, df.species.merge$beta.log10)

# Load ML tree
ml.tree <- read.tree("data/tree/RAxML_bipartitionsBranchLabels.ltde")
# Define the outgroup
outgroup <- match("NC_005042.1.353331-354795", ml.tree$tip.label)
# Create a rooted tree {ape} using the outgroup
ml.rooted <- root(ml.tree, outgroup, resolve.root = TRUE)
# Keep rooted but drop outgroup branch
ml.rooted <- drop.tip(ml.rooted, c("NC_005042.1.353331-354795", "KBS0725", "KBS0727", "KBS0714"))
is.ultrametric(ml.rooted)
ml.rooted.um  <- chronos(ml.rooted)
is.ultrametric(ml.rooted.um)


phylolm.beta <- phylolm(beta.log10 ~ Dormant.anc + dorm.diff,
                                 data= df.species.merge, phy=ml.rooted.um, 
                                 model = 'lambda')

summary(phylolm.beta)

phylolm.alpha <- phylolm(alpha ~ Dormant.anc + dorm.diff,
                                   data= df.species.merge, phy=ml.rooted.um, 
                                 model = 'lambda')
summary(phylolm.alpha)


beta.anc.plot <- ggplot(data = df.species.merge, aes(x = Dormant.anc, y = 10**beta.log10)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Mean scale paramater, $\\bar{\\lambda}$") ) + 
  xlab("Initial proportion of dormant cells") +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw() +
  scale_x_continuous(limits = c(0, 0.8)) +
  theme(axis.title.x = element_text(color="black", size=12), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



beta.rel.plot <- ggplot(data = df.species.merge, aes(x = dorm.diff, y = 10**beta.log10)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Mean scale paramater, $\\bar{\\lambda}$") ) + 
  xlab("Change in proportion of dormant cells") +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.5, 1)) +
  theme(axis.title.x = element_text(color="black", size=12), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


alpha.anc.plot <- ggplot(data = df.species.merge, aes(x = Dormant.anc, y = alpha)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Mean shape paramater, $\\bar{k}$")) +
  xlab("Initial proportion of dormant cells") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 0.8)) +
  theme(axis.title.x = element_text(color="black", size=12), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



alpha.rel.plot <- ggplot(data = df.species.merge, aes(x = dorm.diff, y = alpha)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Mean shape paramater, $\\bar{k}$")) +
  xlab("Change in proportion of dormant cells") +
  theme_bw() +
  scale_x_continuous(limits = c(-0.5, 1)) +
  theme(axis.title.x = element_text(color="black", size=12), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())




g <- ggarrange(beta.anc.plot, beta.rel.plot, alpha.anc.plot, alpha.rel.plot,
               # First row with scatter plot
               # Second row with box and dot plots
               ncol = 2, nrow = 2,
               labels = "auto")#, label.y = c(1, 0.5, 0.25)    

ggsave(file="figs/metab.png", g, width=10,height=10, units='in', dpi=600)




