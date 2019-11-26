rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('ggplot2')
library('ggpubr')
library('latex2exp')
library(scales)




df <- read.table("data/demography/weibull_results_clean.csv", 
                 header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)

d <- density(df$alpha)
plot(d, main="",xlab="Shape parameter")
abline(v = 1, col="red", lwd=3, lty=2)


df$p.value.BH.bin <- df$p.value.BH < 0.05
lr.plot <- ggplot(data = df) +
  geom_hline(yintercept=0, linetype = "longdash", size=2) +
  geom_point(aes(x = reorder(strain, LR), y = LR, color = p.value.BH.bin),  alpha = 0.6, size=5.5) +
  scale_color_manual(values=c("red", "blue")) +
  ylab(TeX("Log-Likelihood ratio of the Weibull vs. the exponential distribution") ) + 
  ylim(-10, 250) +
  coord_flip() +
  theme_bw() + 
  annotate("text", x=20, y=-10, label= "d", size = 8, fontface =2) +
  theme(axis.title.y=element_blank(), 
        axis.text.y=element_text(size = 13), 
        axis.title.x = element_text(color="black", size=20, vjust=0, hjust=0.5), 
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


#g <- ggarrange(KBS0714.R4.plot, KBS0703.R4.plot, KBS0812.R4.plot,                                              # First row with scatter plot
#               # Second row with box and dot plots
#               ncol = 3, nrow = 1,
#               labels = "auto")

g <- ggarrange(ggarrange(KBS0714.R4.plot, KBS0703.R4.plot, KBS0812.R4.plot, ncol = 3, labels = c("a", "b", "c"), font.label = list(size = 22, color = "black")),
               # First row with scatter plot
               lr.plot,
               # Second row with box and dot plots
               nrow = 2, 
               labels = NULL)



ggsave(file="figs/likelihood.png", g, width=12,height=9, units='in', dpi=600)


