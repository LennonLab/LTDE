#### trait regression
rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('phylolm')
library('ggplot2')
library('ggpubr')
library('latex2exp')

df.species <- read.table("data/demography/weibull_results_clean_species.csv", 
                         header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
rownames(df.species) <- df.species$Species
# remove bacillus
#df.species.no_812<-df.species[!(df.species$Species=="KBS0812"),]
df.species<-df.species[!(df.species$Species=="KBS0727"),]

df.iRep <- read.table("data/iRep_clean.txt", 
                   header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
df.iRep.species <- aggregate(df.iRep[, c('iRep')], list(df.iRep$Species), mean)
colnames(df.iRep.species) <- c("Species","iRep")
rownames(df.iRep.species) <- df.iRep.species$Species
df.iRep.merge <- merge(df.species, df.iRep.species, by=0)
rownames(df.iRep.merge) <- df.iRep.merge$Species.x

summary(lm(df.iRep.merge$beta.log10 ~ df.iRep.merge$iRep))

beta.plot <- ggplot(data = df.iRep.merge, aes(x = iRep, y = 10** beta.log10)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  xlab(TeX("Index of replication, $\\iRep$")) +
  ylab(TeX("Scale paramater, $\\lambda$") ) + 
  #scale_y_continuous(limits = c(0, 1)) +
  scale_y_log10(
    limits = c(0.000001, 2000),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14), 
          axis.title.y = element_text(color="black", size=14), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) #+ scale_y_reverse()



summary(lm(df.iRep.merge$alpha ~ df.iRep.merge$iRep ))

alpha.plot <- ggplot(data = df.iRep.merge, aes(x = iRep, y = alpha)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Shape paramater, $\\k$")) +
  xlab(TeX("Index of replication, $\\iRep$")) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())




g <- ggarrange(beta.plot, alpha.plot,                                              # First row with scatter plot
               # Second row with box and dot plots
               ncol = 2, nrow = 1,
               labels = "auto")#, label.y = c(1, 0.5, 0.25)                                     # Labels of the scatter plot


ggsave(file="figs/iRep.png", g,width=10,height=5, units='in', dpi=600)

