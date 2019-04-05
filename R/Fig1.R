rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('ggplot2')
library('phylolm')
library('NPCirc')
library('latex2exp')
library('rr2')

df <- read.csv("data/demography/weibull_results_clean.csv", 
               header = TRUE, stringsAsFactors = FALSE)
df.no_812<-df[!(df$strain=="KBS0812"),]

df.species <- read.table("data/demography/weibull_results_clean_species.csv", 
                 header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
rownames(df.species) <- df.species$Species
# remove bacillus
df.species.no_812<-df.species[!(df.species$Species=="KBS0812"),]
rownames(df.species.no_812) <- df.species.no_812$Species

###### ggplot KDE
bw <- bw.CV(log10(df$mttf), method="LCV", lower=0, upper=20)


kde.plot <- ggplot(df.no_812, aes(log10(mttf))) +
            #scale_x_log10() + 
            xlab(TeX("$\\T_{death}$ (days),  $\\log_{10}$") ) + 
            ylab('Density') +
            geom_density(fill = "blue", alpha = 0.2) +
            geom_vline(xintercept=mean(log10(df.no_812$mttf)), linetype = "longdash") + 
            theme_bw()
            
kde.plot = kde.plot + theme(axis.title.x = element_text(color="black", size=16), 
                            axis.title.y = element_text(color="black", size=16)) +
          scale_x_continuous(breaks = c(0, 1, 2), labels = c(1, 10, 100))

##### ggplot parameter regression
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
# Run a phylogeny-corrected regression with no bootstrap replicates
fit.phy <- phylolm(alpha.mean  ~ log10(beta.mean), data = df.species.no_812, 
                   re.ml.rooted.um.prunned, model = 'lambda', boot = 0)

#r2 = format(R2(fit.phy, phy=re.ml.rooted.um.prunned)[2], digits = 3)
phylo.params <- ggplot(data = df.species.no_812, aes(x = log10(beta.mean), y = alpha.mean)) +
                geom_point(color='blue') +
                xlab(TeX("$\\bar{\\beta}$,  $\\log_{10}$") ) + 
                ylab(TeX("$\\bar{\\alpha}$")) +
                stat_function(fun = function(x) fit.phy$coefficients[1] + fit.phy$coefficients[2] * x) + 
                theme_bw()
phylo.params <- phylo.params + theme(axis.title.x = element_text(color="black", size=16), 
                            axis.title.y = element_text(color="black", size=16))


#### make boxplot ggplot

ggplot(data = df.species.no_812) +
  geom_point(aes(x = Species, y = mttf.mean)) +
  scale_y_log10() + 
  coord_flip() +
  theme_bw()




grid.arrange(bp, arrangeGrob(dp, sc), ncol = 2)


grid.arrange(bp, dp, sc, vp, ncol = 2, 
             layout_matrix = cbind(c(1,1,1), c(2,3,4)))

