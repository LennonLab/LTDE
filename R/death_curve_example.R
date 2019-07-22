rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('ggplot2')
library('ggpubr')
library('latex2exp')
library(scales)

obs <- read.csv("data/demography/longtermdormancy_20190528_nocomments.csv", 
                header = TRUE, stringsAsFactors = FALSE)
## Adding 1 to deal with log(0) observations
obs$Abund <- (as.numeric(obs$Colonies) +1)* (1000 / as.numeric(obs$Inoculum )) * ( 10 ^  as.numeric(obs$Dilution) )
strains <- sort(unique(obs$Strain))
strains <- strains[table(obs$Strain)>10]

KBS0714.R4 <- obs[(obs$Strain =='KBS0714') & (obs$Rep ==4),]
KBS0703.R4 <- obs[(obs$Strain =='KBS0703') & (obs$Rep ==4),]
KBS0812.R4 <- obs[(obs$Strain =='KBS0812') & (obs$Rep ==4),]

KBS0714.R4.time <-(as.numeric(strptime(KBS0714.R4$Firstread_date,format="%d-%b-%y",tz="EST"))-
                     as.numeric(strptime(KBS0714.R4$Dormstart_date,format="%d-%b-%y",tz="EST")))/(3600*24)

KBS0703.R4.time <-(as.numeric(strptime(KBS0703.R4$Firstread_date,format="%d-%b-%y",tz="EST"))-
                     as.numeric(strptime(KBS0703.R4$Dormstart_date,format="%d-%b-%y",tz="EST")))/(3600*24)

KBS0812.R4.time <-(as.numeric(strptime(KBS0812.R4$Firstread_date,format="%d-%b-%y",tz="EST"))-
                     as.numeric(strptime(KBS0812.R4$Dormstart_date,format="%d-%b-%y",tz="EST")))/(3600*24)
KBS0714.R4.N <- KBS0714.R4$Abund
KBS0703.R4.N <- KBS0703.R4$Abund
KBS0812.R4.N <- KBS0812.R4$Abund

df <- read.table("data/demography/weibull_results_clean.csv", 
                 header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)

KBS0714.R4.N0 <- df[(df$strain =='KBS0714') & (df$rep ==4), "N_0"]
KBS0714.R4.alpha <- df[(df$strain =='KBS0714') & (df$rep ==4), "alpha"]
KBS0714.R4.beta <- df[(df$strain =='KBS0714') & (df$rep ==4), "beta"]

KBS0703.R4.N0 <- df[(df$strain =='KBS0703') & (df$rep ==4), "N_0"]
KBS0703.R4.alpha <- df[(df$strain =='KBS0703') & (df$rep ==4), "alpha"]
KBS0703.R4.beta <- df[(df$strain =='KBS0703') & (df$rep ==4), "beta"]

KBS0812.R4.N0 <- df[(df$strain =='KBS0812') & (df$rep ==4), "N_0"]
KBS0812.R4.alpha <- df[(df$strain =='KBS0812') & (df$rep ==4), "alpha"]
KBS0812.R4.beta <- df[(df$strain =='KBS0812') & (df$rep ==4), "beta"]

x.range <- seq(0, 1000, length.out = 1000) 
x.range.KBS0714 <- seq(0, 250, length.out = 1000) 
KBS0714.R4.weib.y <- log(KBS0714.R4.N0) - ((x.range.KBS0714/KBS0714.R4.beta) ^ KBS0714.R4.alpha)
KBS0703.R4.weib.y <- max(KBS0714.R4.N) *  exp(-1* ((x.range/KBS0703.R4.beta) ^ KBS0703.R4.alpha) )
KBS0703.R4.weib.y <- log(KBS0812.R4.N0) - ((x.range/KBS0812.R4.beta) ^ KBS0812.R4.alpha)

KBS0714.R4.exp.y <- log(KBS0714.R4.N0) - (x.range.KBS0714/KBS0714.R4.beta)
KBS0703.R4.exp.y <- log(KBS0703.R4.N0) - (x.range/KBS0703.R4.beta) 
KBS0812.R4.exp.y <- log(KBS0812.R4.N0) - (x.range/KBS0812.R4.beta)

KBS0714.R4.df <- do.call(rbind, Map(data.frame, time=KBS0714.R4.time, N_time=KBS0714.R4.N))
KBS0703.R4.df <- do.call(rbind, Map(data.frame, time=KBS0703.R4.time, N_time=KBS0703.R4.N))
KBS0812.R4.df <- do.call(rbind, Map(data.frame, time=KBS0812.R4.time, N_time=KBS0812.R4.N))


KBS0714.R4.plot <- ggplot(KBS0714.R4.df, aes(x=time, y=N_time)) + 
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Number of individuals, $N(t)$")) +
  xlab(TeX("Days, $t$")) + 
  stat_function(fun = function(x) log(KBS0714.R4.N0) + -1* x/KBS0714.R4.beta, linetype="dashed", color="darkgrey", size=1.2) +
  stat_function(fun = function(x) log(KBS0714.R4.N0) + -1* ((x/KBS0714.R4.beta) ** KBS0714.R4.alpha), size=1.2 ) +
  scale_y_continuous(trans = log_trans(), limits = c(1, 2500000),
                     breaks = trans_breaks("log", function(x) exp(x)),
                     labels = trans_format("log", math_format(e^.x))) +
  xlim(-10, 300) +
  theme_bw() +
  annotate("text", x=250, y=2000000, label=TeX(sprintf("$\\lambda = %g$", round(KBS0714.R4.beta,2)))) +
  annotate("text", x=245, y=550000, label=TeX(sprintf("$k = %g$", round(KBS0714.R4.alpha,2)))) +
  
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size=14)) +
  ggtitle(expression(paste(italic("Micrococcus"), " sp. KBS0714")))



KBS0703.R4.plot <- ggplot(KBS0703.R4.df, aes(x=time, y=N_time)) + 
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Number of individuals, $N(t)$")) +
  xlab(TeX("Days, $t$")) + 
  stat_function(fun = function(x) log(KBS0703.R4.N0) + -1* x/KBS0703.R4.beta, linetype="dashed", color="darkgrey", size=1.2) +
  stat_function(fun = function(x) log(KBS0703.R4.N0) + -1* ((x/KBS0703.R4.beta) ** KBS0703.R4.alpha), size=1.2 ) +
  scale_y_continuous(trans = log_trans(), limits = c(600000, 9e+08),
                     breaks = trans_breaks("log", function(x) exp(x)),
                     labels = trans_format("log", math_format(e^.x))) +
  xlim(-10, 1010) +
  theme_bw() +
  annotate("text", x=850, y=6.5e+08, label=TeX(sprintf("$\\lambda = %g$", round(KBS0703.R4.beta,2)))) +
  annotate("text", x=830, y=3.5e+08, label=TeX(sprintf("$k = %g$", round(KBS0703.R4.alpha,2)))) +
  
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size=14)) +
  ggtitle(expression(paste(italic("Arthrobacter"), " sp. KBS0703")))


KBS0812.R4.plot <- ggplot(KBS0812.R4.df, aes(x=time, y=N_time)) + 
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Number of individuals, $N(t)$")) +
  xlab(TeX("Days, $t$")) + 
  # rate of decay is so high, extinction would almost be instantaneous, so just plot line at zero
  geom_vline(xintercept = 0, linetype="dashed", color="darkgrey", size=1.2) +
  #stat_function(fun = function(x) log(KBS0812.R4.N0) + -1* x/KBS0812.R4.beta, linetype="dashed", color="darkgrey", size=1.2) +
  stat_function(fun = function(x) log(KBS0812.R4.N0) + -1* ((x/KBS0812.R4.beta) ** KBS0812.R4.alpha), size=1.2 ) +
  scale_y_continuous(trans = log_trans(), limits = c(600000, 9e+08),
                     breaks = trans_breaks("log", function(x) exp(x)),
                     labels = trans_format("log", math_format(e^.x))) +
  xlim(-10, 1010) +
  theme_bw() +
  annotate("text", x=850, y=6.5e+08, label=TeX(sprintf("$\\lambda = 2.68e-07$"))) +
  annotate("text", x=780, y=3.5e+08, label=TeX(sprintf("$k = %g$", round(KBS0812.R4.alpha,2)))) +
  
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size=14)) +
  ggtitle(expression(paste(italic("Bacillus"), " sp. KBS0812")))








g <- ggarrange(KBS0714.R4.plot, KBS0703.R4.plot, KBS0812.R4.plot,                                              # First row with scatter plot
               # Second row with box and dot plots
               ncol = 3, nrow = 1,
               labels = "auto")


ggsave(file="figs/death_curve_example.png", g, width=9.3,height=3.1, units='in', dpi=600)


