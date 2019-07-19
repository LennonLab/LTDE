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

lr.plot <- ggplot(data = df.weib) +
  geom_hline(yintercept=1, linetype = "longdash", size=2) +
  geom_point(aes(x = reorder(strain, LR), y = LR), color='blue', alpha = 0.6, size=8) +
  ylab(TeX("Likelihood-ratio") ) + 
  ylim(-10, 250) +
  coord_flip() +
  theme_bw() + 
  annotate("text", x=22, y=250, label= "a", size = 12) +
  theme(axis.title.y=element_blank(), 
        axis.ticks.y=element_blank(),
        #axis.text.y=element_text(size = 9), 
        axis.text.y=element_blank(),
        axis.text.x=element_text(size = 18),
        #axis.title.x = element_blank(),
        axis.title.x = element_text(color="black", size=20),#, size=11, vjust=0, hjust=0.5), 
        #axis.title.y = theme_text(vjust=-0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

p.value.plot <- ggplot(data = df.weib) +
  geom_hline(yintercept= 0.05, linetype = "longdash", size=2) +
  geom_point(aes(x = reorder(strain, LR), y = p.value.BH), color='blue', alpha = 0.6, size=8) +
  ylab(TeX("Benjamini-Hochberg corrected p-value") ) + 
  #ylim(-10, 250) +
  
  coord_flip() +
  theme_bw() + 
  theme(axis.title.y=element_blank(), 
        axis.text.y=element_text(size = 18), 
        axis.title.x = element_text(color="black", size=20, vjust=0, hjust=0.5), 
        axis.text.x = element_text(size=16),
        #axis.title.y = theme_text(vjust=-0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  annotate("text", x=22, y=1e-5, label= "b", size = 12) +
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


g.2 <- plot_grid(lr.plot, p.value.plot, labels = c('a', 'b'), 
                 label_size = 12, label_x = c(40,100), 
                 label_y = c(1000,2),align = "h", 
                 nrow = 1, rel_widths = c(0.4, 0.6))


ggsave(file="figs/exp_vs_weib.png", g.2, width=15,height=10,units='in', dpi=600)




