#### trait regression
rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('phylolm')
library('ggplot2')
library("ape")
library('latex2exp')
library('ggpubr')
library('plotrix')

df.species <- read.table("data/demography/weibull_results_clean_species.csv", 
                         header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
rownames(df.species) <- df.species$Species
div <-  read.table("data/breseq/genetic_diversity.txt",sep = "\t",header = TRUE)

# get mean of log 10 transformed observations
div.mean <- aggregate(div[, c('dN_dS_total', 'dN_dS_fixed', 'mean_freq')], list(div$Species), mean)
colnames(div.mean)[1] <- "Species"
colnames(div.mean)[2] <- "dN_dS_total"
colnames(div.mean)[3] <- "dN_dS_fixed"
colnames(div.mean)[4] <- "mean_freq"


div.dNdS.merge <- merge(df.species, div.mean, by="Species")


dN_dS.se <- aggregate(div[, c('dN_dS_total', 'dN_dS_fixed')], list(div$Species), std.error)
colnames(dN_dS.se) <- c("Species", "dN_dS_total.se", "dN_dS_fixed.se" )
div.dNdS.merge <- merge(div.dNdS.merge, dN_dS.se, by="Species")


# make dN/dS plot
boxplot.dNdS <- ggplot(data = div.dNdS.merge) +
  geom_point(aes(x = reorder(Species, -dN_dS_total), y = dN_dS_total), color='blue', alpha = 0.6, size=8) +
  geom_point(aes(x = reorder(Species, -dN_dS_total), y = (dN_dS_total-dN_dS_total.se)), shape=124,size=8 ) +
  geom_point(aes(x = reorder(Species, -dN_dS_total), y = (dN_dS_total+dN_dS_total.se)), shape=124,size=8 ) +
  geom_segment(aes(x = Species, y = dN_dS_total, xend = Species, yend = (dN_dS_total-dN_dS_total.se)), size = 1.5) +
  geom_segment(aes(x = Species, y = dN_dS_total, xend = Species, yend = (dN_dS_total+dN_dS_total.se)), size = 1.5) +
  geom_hline(yintercept= 1, linetype = "longdash", size=2) +
  ylab(TeX("Ratio of nonsynonymous to synonymous mutations, $dN/dS$") ) + 
  coord_flip() +
  theme_bw() + 
  theme(axis.title.y=element_blank(), 
        axis.text.y=element_text(size = 18), 
        axis.title.x = element_text(color="black", size=24, vjust=0, hjust=0.5), 
        axis.text.x = element_text(size=18),
        #axis.title.y = theme_text(vjust=-0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none") +

  scale_x_discrete( position = "top",
                    labels=c("KBS0707" = expression(paste(italic("Pseudomonas"), " sp. KBS0707")), 
                             "KBS0702" = expression(paste(italic("Arthrobacter"), " sp. KBS0702")),
                             "ATCC13985" = expression(paste(italic("Pseudomonas"), " sp. ATCC13985")),
                             "KBS0711" = expression(paste(italic("Janthinobacterium"), " sp. KBS0711")),
                             "KBS0712" = expression(paste(italic("Variovorax"), " sp. KBS0712")),
                             "KBS0713" = expression(paste(italic("Yersinia"), " sp. KBS0713")),
                             "KBS0715" = expression(paste(italic("Curtobacterium"), " sp. KBS0715")),
                             "KBS0721" = expression(paste(italic("Flavobacterium"), " sp. KBS0721")),
                             "KBS0801" = expression(paste(italic("Burkholderia"), " sp. KBS0801")),
                             "KBS0802" = expression(paste(italic("Pseudomonas"), " sp. KBS0802")),
                             "KBS0812" = expression(paste(italic("Bacillus"), " sp. KBS0812"))))


ggsave(file="figs/dNdS_total.png", boxplot.dNdS, width=13,height=10, units='in', dpi=600)

