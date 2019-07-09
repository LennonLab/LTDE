#### trait regression
rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('phylolm')
library('ggplot2')
library('ggpubr')

df.species <- read.table("data/demography/weibull_results_clean_species.csv", 
                         header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
rownames(df.species) <- df.species$Species
# remove bacillus
df.species.no_812<-df.species[!(df.species$Species=="KBS0812"),]
df.species.no_812<-df.species.no_812[!(df.species.no_812$Species=="KBS0727"),]

df.iRep <- read.table("data/iRep_clean.txt", 
                   header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
df.iRep.species <- aggregate(df.iRep[, c('iRep')], list(df.iRep$Species), mean)
colnames(df.iRep.species) <- c("Species","iRep")
rownames(df.iRep.species) <- df.iRep.species$Species
df.iRep.merge <- merge(df.species.no_812, df.iRep.species, by=0)
rownames(df.iRep.merge) <- df.iRep.merge$Species.x

summary(lm(df.iRep.merge$iRep ~ df.iRep.merge$alpha))

beta.plot <- ggplot(data = df.iRep.merge, aes(x = 10** beta.log10, y = iRep)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  xlab(TeX("Mean scale paramater, $\\bar{\\lambda}$") ) + 
  ylab(TeX("Mean index of replication, $\\bar{\\iRep}$")) +
  #scale_y_continuous(limits = c(0, 1)) +
  scale_x_log10(
    limits = c(0.0001, 1000),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14), 
          axis.title.y = element_text(color="black", size=14), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + scale_y_reverse()




alpha.plot <- ggplot(data = df.iRep.merge, aes(x = alpha, y = iRep)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  xlab(TeX("Mean shape paramater, $\\bar{k}$")) +
  ylab(TeX("Mean index of replication, $\\bar{\\iRep}$")) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + scale_y_reverse()




g <- ggarrange(beta.plot, alpha.plot,                                              # First row with scatter plot
               # Second row with box and dot plots
               ncol = 2, nrow = 1,
               labels = "auto")#, label.y = c(1, 0.5, 0.25)                                     # Labels of the scatter plot


ggsave(file="figs/iRep.png", g,width=10,height=5, units='in', dpi=600)

