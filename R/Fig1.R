rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('ggplot2')
library('phylolm')
library('NPCirc')
library('latex2exp')
library('rr2')
library('dplyr')
library('gridExtra')
library('ggpubr')
library('grid')

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
            
kde.plot = kde.plot + theme(axis.title.x = element_text(color="black", size=14), 
                            axis.title.y = element_text(color="black", size=14)) +
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
                geom_point(color='blue', alpha = 0.6, size=4) +
                xlab(TeX("$\\bar{\\beta}$,  $\\log_{10}$") ) + 
                ylab(TeX("$\\bar{\\alpha}$")) +
                stat_function(fun = function(x) fit.phy$coefficients[1] + fit.phy$coefficients[2] * x) + 
                theme_bw()
phylo.params <- phylo.params + theme(axis.title.x = element_text(color="black", size=14), 
                            axis.title.y = element_text(color="black", size=14))


#### make boxplot ggplot
boxplot <- ggplot(data = df.species.no_812) +
  geom_point(aes(x = reorder(Species, -mttf.mean), y = mttf.mean), color='blue', alpha = 0.6, size=2.2) +
  geom_point(aes(x = reorder(Species, -mttf.mean), y = mttf.CI.2.5), shape=124,size=2.5 ) +
  geom_point(aes(x = reorder(Species, -mttf.mean), y = mttf.CI.97.5), shape=124,size=2.5 ) +
  geom_segment(aes(x = Species, y = mttf.mean, xend = Species, yend = mttf.CI.2.5), size = 0.5) +
  geom_segment(aes(x = Species, y = mttf.mean, xend = Species, yend = mttf.CI.97.5), size = 0.5) +
  ylab(TeX("$\\T_{death}$ (days),  $\\log_{10}$") ) + 
  scale_y_log10() + 
  coord_flip() +
  theme_bw() + theme(axis.title.y=element_blank(), axis.text.y=element_text(size = 7), axis.title.x = element_text(color="black", size=14)) +
  scale_x_discrete( position = "top",
    labels=c("KBS0707" = expression(paste(italic("Pseudomonas"), " sp. KBS0707")), 
             "KBS0702" = expression(paste(italic("Arthrobacter"), " sp. KBS0702")),
             "ATCC13985" = expression(paste(italic("Pseudomonas"), " sp. ATCC13985")),
             "ATCC43928" = expression(paste(italic("Pseudomonas"), " sp. ATCC43928")),
             "KBS0701" = expression(paste(italic("Pedobacter"), " sp. KBS0701")),
             "KBS0703" = expression(paste(italic("Arthrobacter"), " sp. KBS0703")),
             "KBS0705" = expression(paste(italic("Inquilinus"), " sp. KBS0705")),
             "KBS0706" = expression(paste(italic("Mycobacterium"), " sp. KBS0706")),
             "KBS0710" = expression(paste(italic("Pseudomonas"), " sp. KBS0710")),
             "KBS0711" = expression(paste(italic("Janthinobacterium"), " sp. KBS0711")),
             "KBS0712" = expression(paste(italic("Variovorax"), " sp. KBS0712")),
             "KBS0713" = expression(paste(italic("Yersinia"), " sp. KBS0713")),
             "KBS0714" = expression(paste(italic("Arthrobacter"), " sp. KBS0714")),
             "KBS0715" = expression(paste(italic("Curtobacterium"), " sp. KBS0715")),
             "KBS0721" = expression(paste(italic("Flavobacterium"), " sp. KBS0721")),
             "KBS0722" = expression(paste(italic("Oerskovia"), " sp. KBS0722")),
             "KBS0724" = expression(paste(italic("Rhodococcus"), " sp. KBS0724")),
             "KBS0725" = expression(paste(italic("Bradyrhizobium"), " sp. KBS0725")),
             "KBS0727" = expression(paste(italic("Bradyrhizobium"), " sp. KBS0727")),
             "KBS0801" = expression(paste(italic("Burkholderia"), " sp. KBS0801")),
             "KBS0802" = expression(paste(italic("Pseudomonas"), " sp. KBS0802"))) )
  


g <- ggarrange(boxplot,                                                 # First row with scatter plot
          ggarrange(kde.plot, phylo.params, ncol = 2, labels = c("B", "C")), # Second row with box and dot plots
          nrow = 2, 
          labels = "A"                                        # Labels of the scatter plot
) 


ggsave(file="figs/Fig1.png", g, units='in', dpi=600)



#grid.arrange(kde.plot, phylo.params, boxplot, ncol=2, nrow =2)


# Move to a new page
#grid.newpage()
# Create layout : nrow = 3, ncol = 2
#pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 2)))
# A helper function to define a region on the layout
#define_region <- function(row, col){
#  viewport(layout.pos.row = row, layout.pos.col = col)
#} 
# Arrange the plots
#print(boxplot, vp = define_region(row = 1, col = 1:2))   # Span over two columns
#print(kde.plot, vp = define_region(row = 2, col = 1))
#print(phylo.params, vp = define_region(row = 2, col = 2))

