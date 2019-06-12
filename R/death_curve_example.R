rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('ggplot2')
library('ggpubr')
library('latex2exp')

obs <- read.csv("data/demography/longtermdormancy_20190528_nocomments.csv", 
                header = TRUE, stringsAsFactors = FALSE)
## Adding 1 to deal with log(0) observations
obs$Abund <- (as.numeric(obs$Colonies) +1)* (1000 / as.numeric(obs$Inoculum )) * ( 10 ^  as.numeric(obs$Dilution) )
strains <- sort(unique(obs$Strain))
strains <- strains[table(obs$Strain)>10]

KBS0812.R4 <- obs[(obs$Strain =='KBS0812') & (obs$Rep ==4),]
KBS0714.R3 <- obs[(obs$Strain =='KBS0714') & (obs$Rep ==4),]
KBS0715.R4 <- obs[(obs$Strain =='KBS0715') & (obs$Rep ==4),]

KBS0812.R4.time <-(as.numeric(strptime(KBS0812.R4$Firstread_date,format="%d-%b-%y",tz="EST"))-
                     as.numeric(strptime(KBS0812.R4$Dormstart_date,format="%d-%b-%y",tz="EST")))/(3600*24)

KBS0714.R3.time <-(as.numeric(strptime(KBS0714.R3$Firstread_date,format="%d-%b-%y",tz="EST"))-
                     as.numeric(strptime(KBS0714.R3$Dormstart_date,format="%d-%b-%y",tz="EST")))/(3600*24)

KBS0715.R4.time <-(as.numeric(strptime(KBS0715.R4$Firstread_date,format="%d-%b-%y",tz="EST"))-
                     as.numeric(strptime(KBS0715.R4$Dormstart_date,format="%d-%b-%y",tz="EST")))/(3600*24)

KBS0812.R4.time.logRel <- log(KBS0812.R4$Abund / max(KBS0812.R4$Abund))
KBS0714.R3.time.logRel <- log(KBS0714.R3$Abund / max(KBS0714.R3$Abund))
KBS0715.R4.time.logRel <- log(KBS0715.R4$Abund / max(KBS0715.R4$Abund))

KBS0812.R4.df <- do.call(rbind, Map(data.frame, time=KBS0812.R4.time, logRel=KBS0812.R4.time.logRel))
KBS0714.R3.df <- do.call(rbind, Map(data.frame, time=KBS0714.R3.time, logRel=KBS0714.R3.time.logRel))
KBS0715.R4.df <- do.call(rbind, Map(data.frame, time=KBS0715.R4.time, logRel=KBS0715.R4.time.logRel))

df <- read.table("data/demography/weibull_results_clean.csv", 
                 header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
KBS0812.R4.alpha <- df[(df$strain =='KBS0812') & (df$rep ==4), "alpha"]
KBS0812.R4.beta <- df[(df$strain =='KBS0812') & (df$rep ==4), "beta"]

KBS0714.R3.alpha <- df[(df$strain =='KBS0714') & (df$rep ==4), "alpha"]
KBS0714.R3.beta <- df[(df$strain =='KBS0714') & (df$rep ==4), "beta"]

KBS0715.R4.alpha <- df[(df$strain =='KBS0715') & (df$rep ==4), "alpha"]
KBS0715.R4.beta <- df[(df$strain =='KBS0715') & (df$rep ==4), "beta"]

KBS0812.R4.fxn <- function(x) -1 * ((x/KBS0812.R4.beta) ** KBS0812.R4.alpha)
KBS0714.R3.fxn <- function(x) -1 * ((x/KBS0714.R3.beta) ** KBS0714.R3.alpha)
KBS0715.R4.fxn <- function(x) -1 * ((x/KBS0715.R4.beta) ** KBS0715.R4.alpha)

#stat_function(fun = KBS0812.R4.fxn, size=1.5, lty=2, col = "black") +
KBS0812.R4.plot <- ggplot(KBS0812.R4.df, aes(x=time, y=logRel)) + 
                  geom_point(color='blue', alpha = 0.6, size=4) +
                  ylab(TeX("Log Survivorship, $ln(S(t))$")) +
                  xlab(TeX("Days, $t$")) + 
                  #scale_y_continuous(limits = c(0, 1)) +
                  theme_bw() +
                  annotate("text", x=300, y=1, label=expression(paste(italic("Bacillus"), " sp. KBS0812"))) +
                  theme(axis.title.x = element_text(color="black", size=14), 
                        axis.title.y = element_text(color="black", size=14), 
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank())
                  
                  
                  
KBS0714.R3.plot <- ggplot(KBS0714.R3.df, aes(x=time, y=logRel)) + 
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Log Survivorship, $ln(S(t))$")) +
  xlab(TeX("Days, $t$")) + 
  theme_bw() +
  annotate("text", x=100, y=1, label=expression(paste(italic("Arthrobacter"), " sp. KBS0714"))) +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



KBS0715.R4.plot <- ggplot(KBS0715.R4.df, aes(x=time, y=logRel)) + 
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Log Survivorship, $ln(S(t))$")) +
  xlab(TeX("Days, $t$")) + 
  theme_bw() +
  annotate("text", x=450, y=1, label=expression(paste(italic("Curtobacterium"), " sp. KBS0715"))) +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


g <- ggarrange(KBS0812.R4.plot, KBS0714.R3.plot, KBS0715.R4.plot,                                              # First row with scatter plot
               # Second row with box and dot plots
               ncol = 3, nrow = 1,
               labels = "auto")


ggsave(file="figs/death_curve_example.png", g, width=9.3,height=3.1, units='in', dpi=600)



