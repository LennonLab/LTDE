rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library(ggplot2)
library('latex2exp')

df.weib <- read.table("data/demography/weibull_results_clean.csv", 
                      header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
colnames(df.weib)[6] <- "AIC_weib"
df.weib <- df.weib[c("strain", "rep", "AIC_weib")]

df.gomp <- read.table("data/demography/gompertz_results_clean.csv", 
                 header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
colnames(df.gomp)[6] <- "AIC_gomp"
df.gomp <- df.gomp[c("strain", "rep", "AIC_gomp")]

df.merge <- merge(df.weib, df.gomp, by=c("strain","rep")) # NA's match
df.merge <- transform(df.merge, min = pmin(AIC_weib, AIC_gomp))
df.merge$aic_log_likelihood <- (df.merge$min - df.merge$AIC_gomp ) /2

df.merge <- df.merge[!(df.merge$strain == "KBS0711" & df.merge$rep == 1 ),] 
df.merge <- df.merge[!(df.merge$strain == "KBS0711" & df.merge$rep == 2 ),] 
df.merge <- df.merge[!(df.merge$strain == "KBS0711" & df.merge$rep == 3 ),] 



aic.plot <- ggplot(data = df.merge) +
  geom_point(aes(x = reorder(strain, aic_log_likelihood), y = aic_log_likelihood), color='blue', alpha = 0.6, size=2.2) +
  ylab(TeX("Relative log likelihood \n of Gompertz distribution") ) + 
  ylim(-275, 10) +
  geom_hline(yintercept= 0, linetype = "longdash") +
  coord_flip() +
  theme_bw() + 
  theme(axis.title.y=element_blank(), 
        axis.text.y=element_text(size = 9), 
        axis.title.x = element_text(color="black", size=11, vjust=-2.5), 
        #axis.title.y = theme_text(vjust=-0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
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

ggsave(file="figs/aic.png", aic.plot, units='in', dpi=600)


