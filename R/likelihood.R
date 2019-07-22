rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('ggplot2')
library('latex2exp')
library('reshape2')
#library('cowplot')
library('ggpubr')
df.weib <- read.table("data/demography/weibull_results_clean.csv", 
                      header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)

df.weib$p.value.BH.bin <- df.weib$p.value.BH < 0.05
lr.plot <- ggplot(data = df.weib) +
  geom_hline(yintercept=0, linetype = "longdash", size=2) +
  geom_point(aes(x = reorder(strain, LR), y = LR, color = p.value.BH.bin),  alpha = 0.6, size=8) +
  scale_color_manual(values=c("red", "blue")) +
  ylab(TeX("Log-Likelihood ratio of the Weibull vs. the exponential distribution") ) + 
  ylim(-10, 250) +
  coord_flip() +
  theme_bw() + 
  #annotate("text", x=22, y=250, label= "a", size = 12) +
  theme(axis.title.y=element_blank(), 
        axis.text.y=element_text(size = 18), 
        axis.title.x = element_text(color="black", size=25, vjust=0, hjust=0.5), 
        axis.text.x = element_text(size=18),
        #axis.title.y = theme_text(vjust=-0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none") +
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
                             "KBS0714" = expression(paste(italic("Micrococcus"), " sp. KBS0714")),
                             "KBS0715" = expression(paste(italic("Curtobacterium"), " sp. KBS0715")),
                             "KBS0721" = expression(paste(italic("Flavobacterium"), " sp. KBS0721")),
                             "KBS0722" = expression(paste(italic("Oerskovia"), " sp. KBS0722")),
                             "KBS0724" = expression(paste(italic("Rhodococcus"), " sp. KBS0724")),
                             "KBS0725" = expression(paste(italic("Bradyrhizobium"), " sp. KBS0725")),
                             "KBS0727" = expression(paste(italic("Bradyrhizobium"), " sp. KBS0727")),
                             "KBS0801" = expression(paste(italic("Burkholderia"), " sp. KBS0801")),
                             "KBS0802" = expression(paste(italic("Pseudomonas"), " sp. KBS0802")),
                             "KBS0812" = expression(paste(italic("Bacillus"), " sp. KBS0812"))))

ggsave(file="figs/exp_vs_weib.png", lr.plot, width=15,height=10,units='in', dpi=600)









#theme(axis.title.y=element_blank(), 
#      axis.ticks.y=element_blank(),
#      #axis.text.y=element_text(size = 9), 
#      axis.text.y=element_blank(),
#      axis.text.x=element_text(size = 18),
#      #axis.title.x = element_blank(),
#      axis.title.x = element_text(color="black", size=20),#, size=11, vjust=0, hjust=0.5), 
#      #axis.title.y = theme_text(vjust=-0.5),
#      panel.grid.major = element_blank(), 
#      panel.grid.minor = element_blank()) +





