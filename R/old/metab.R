rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('ggplot2')
library('latex2exp')


df.metab <- read.table("data/metab_stain/staining.all.new.txt", sep="\t", 
                                 header=TRUE, stringsAsFactors = FALSE)

df.metab.noNA <- df.metab[complete.cases(df.metab), ]
df.metab.noNA.anc <- df.metab.noNA[df.metab.noNA$hist == "anc", ] 
df.metab.noNA.anc.mean <- aggregate(df.metab.noNA.anc[,4:6], list(df.metab.noNA.anc$strain), mean )
colnames(df.metab.noNA.anc.mean) <- c("Species", "Active.anc", "Dead.anc", "Dormant.anc")
# derived
df.metab.noNA.der <- df.metab.noNA[df.metab.noNA$hist == "der", ] 
colnames(df.metab.noNA.der)[1] <- "Species"
df.metab.merge <- merge(df.metab.noNA.der, df.metab.noNA.anc.mean, by ="Species")
df.metab.merge$dead.diff <- df.metab.merge$dead - df.metab.merge$Dead.anc
df.metab.merge.mean <- aggregate(df.metab.merge[,10], list(df.metab.merge$Species), mean )
colnames(df.metab.merge.mean) <- c("Species", "change_fraction_dead")

t.test(df.metab.merge.mean$change_fraction_dead, mu=0, alternative="greater")


kde.plot <- ggplot(df.metab.merge.mean, aes(df.metab.merge.mean$change_fraction_dead)) +
  xlab(TeX("Change in proportion of dead cells") ) + 
  ylab('Density') +
  geom_density(fill = "blue", alpha = 0.2) +
  theme_bw() +
  geom_vline(xintercept= 0, linetype = "longdash", size=1.5)+
  xlim(-1, 1)


kde.plot <- kde.plot + theme(axis.title.x = element_text(color="black", size=24), 
                             axis.title.y = element_text(color="black", size=24), 
                             panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank())


ggsave(file="figs/change_fract_dead.png", kde.plot, width=10,height=10, units='in', dpi=600)




