rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('ggplot2')
library('latex2exp')
library('gridExtra')
library('ggpubr')
library('viridis')

obs <- read.csv("data/demography/longtermdormancy_20190528_nocomments.csv", 
                header = TRUE, stringsAsFactors = FALSE)
df <- read.csv("data/demography/weibull_results_clean.csv", 
               header = TRUE, stringsAsFactors = FALSE)
df.species <- read.table("data/demography/weibull_results_clean_species.csv", 
                 header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
rownames(df.species) <- df.species$Species
df.species <- df.species[!(df.species$Species=="KBS0727"),]
df.species <- df.species[order(df.species$mttf.log10),]
# get colors
colors <- viridis(dim(df.species)[1], option = 'viridis')
col.df <- data.frame("Color" =colors)
col.df$strain <- df.species$Species
rownames(col.df) <- df.species$Species
df.species <- merge(df.species, col.df, by="row.names",all.x=TRUE)
df <- merge(df, col.df, by="strain",all.x=TRUE)


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
  geom_point(color='#1E9C89FF', alpha = 0.6, size=4) +
  ylab(TeX("Number of individuals, $N(t)$")) +
  xlab(TeX("Days, $t$")) + 
  stat_function(fun = function(x) log10(KBS0714.R4.N0 * exp(-1* x/KBS0714.R4.beta)), linetype="dashed", color="darkgrey", size=1.2) +
  stat_function(fun = function(x) log10(KBS0714.R4.N0 * exp( -1* ((x/KBS0714.R4.beta) ** KBS0714.R4.alpha))), size=1.2 ) +
  scale_y_log10(
    limits = c(10**1, 10**7),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  
  xlim(-10, 1010) +
  theme_bw() +
  annotate("text", x=850, y=4000000, label=TeX(sprintf("$\\mathit{\\lambda} = %g$", round(KBS0714.R4.beta,2))), size = 5) +
  annotate("text", x=830, y=950000, label=TeX(sprintf("$\\mathit{k} = %g$", round(KBS0714.R4.alpha,2))), size = 5) +
  
  theme(axis.title.x = element_text(color="black", size=18), 
        axis.title.y = element_text(color="black", size=18), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size=18)) +
  ggtitle(expression(paste(italic("Micrococcus"), " sp. KBS0714")))


KBS0703.R4.plot <- ggplot(KBS0703.R4.df, aes(x=time, y=N_time)) + 
  geom_point(color='#2FB47CFF', alpha = 0.6, size=4) +
  ylab(TeX("Number of individuals, $N(t)$")) +
  xlab(TeX("Days, $t$")) + 
  stat_function(fun = function(x) log10(KBS0703.R4.N0 * exp(-1* x/KBS0703.R4.beta)), linetype="dashed", color="darkgrey", size=1.2) +
  stat_function(fun = function(x) log10(KBS0703.R4.N0 * exp( -1* ((x/KBS0703.R4.beta) ** KBS0703.R4.alpha))), size=1.2 ) +
  #scale_y_continuous(trans = log_trans(), limits = c(600000, 9e+08),
  #                   breaks = trans_breaks("log", function(x) exp(x)),
  #                   labels = trans_format("log", math_format(e^.x))) +
  scale_y_log10(
    limits = c(10**6, 10**9),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  
  xlim(-10, 1010) +
  theme_bw() +
  annotate("text", x=850, y=6.5e+08, label=TeX(sprintf("$\\mathit{\\lambda} = %g$", round(KBS0703.R4.beta,2))), size = 5) +
  annotate("text", x=830, y=3.5e+08, label=TeX(sprintf("$\\mathit{k} = %g$", round(KBS0703.R4.alpha,2))), size = 5) +
  
  theme(axis.title.x = element_text(color="black", size=18), 
        axis.title.y = element_text(color="black", size=18), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size=18)) +
  ggtitle(expression(paste(italic("Arthrobacter"), " sp. KBS0703")))


KBS0812.R4.plot <- ggplot(KBS0812.R4.df, aes(x=time, y=N_time)) + 
  geom_point(color='#FDE725FF', alpha = 0.6, size=4) +
  ylab(TeX("Number of individuals, $N(t)$")) +
  xlab(TeX("Days, $t$")) + 
  # rate of decay is so high, extinction would almost be instantaneous, so just plot line at zero
  geom_vline(xintercept = 0, linetype="dashed", color="darkgrey", size=1.2) +
  stat_function(fun = function(x) log10(KBS0812.R4.N0 * exp( -1* ((x/KBS0812.R4.beta) ** KBS0812.R4.alpha))), size=1.2 ) +
  scale_y_log10(
    limits = c(10**6, 10**9),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +

  xlim(-10, 1010) +
  theme_bw() +
  annotate("text", x=850, y=6.5e+08, label=TeX(sprintf("$\\mathit{\\lambda} = 2.68e-07$")), size = 5) +
  annotate("text", x=780, y=3.5e+08, label=TeX(sprintf("$\\mathit{k} = %g$", round(KBS0812.R4.alpha,2))), size = 5) +
  
  theme(axis.title.x = element_text(color="black", size=18), 
        axis.title.y = element_text(color="black", size=18), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size=18)) +
  ggtitle(expression(paste(italic("Bacillus"), " sp. KBS0812")))




#### make boxplot ggplot
boxplot <- ggplot(data = df.species) +
  geom_point(aes(x = reorder(Species, -mttf.log10), y = 10**mttf.log10), color=df.species$Color, alpha = 1, size=5.5) +
  geom_point(aes(x = reorder(Species, -mttf.log10), y = 10**(mttf.log10-mttf.log10.se)), shape=124,size=4.5) +
  geom_point(aes(x = reorder(Species, -mttf.log10), y = 10**(mttf.log10+mttf.log10.se)), shape=124,size=4.5) +
  geom_segment(aes(x = Species, y = 10**mttf.log10, xend = Species, yend = 10**(mttf.log10-mttf.log10.se)), size = 1) +
  geom_segment(aes(x = Species, y = 10**mttf.log10, xend = Species, yend = 10**(mttf.log10+mttf.log10.se)), size = 1) +
  ylab(TeX("Mean time to death, $\\mathit{\\bar{T}_{d}}$ (days)") ) + 
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  coord_flip() +
  geom_hline(yintercept= 10**mean(df.species$mttf.log10), linetype = "longdash", color="black") +
  theme_bw() + 
  annotate("text", x=19, y=0.9, label= "d", size = 9, fontface =2) +
  theme(axis.title.y=element_blank(), 
        axis.text.y=element_text(size = 13), 
        axis.title.x = element_text(color="black", size=24), 
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



g <- ggarrange(ggarrange(KBS0714.R4.plot, KBS0703.R4.plot, KBS0812.R4.plot, ncol = 3, labels = c("a", "b", "c"), font.label = list(size = 22, color = "black")),
               # First row with scatter plot
               boxplot,
               # Second row with box and dot plots
               nrow = 2, 
               labels = NULL)


ggsave(file="figs/mttd_fig1.pdf", g, width=12,height=9, device='pdf', units='in', dpi=180)



######
# time to extinction plot
######





#### make boxplot ggplot
boxplot.T_ext <- ggplot(data = df.species) +
  geom_point(aes(x = reorder(Species, -T_ext.log10), y = (10**T_ext.log10)/365), color=df.species$Color, alpha = 1, size=5.5) +
  geom_point(aes(x = reorder(Species, -T_ext.log10), y = (10**(T_ext.log10-T_ext.log10.se))/365), shape=124,size=4.5) +
  geom_point(aes(x = reorder(Species, -T_ext.log10), y = (10**(T_ext.log10+T_ext.log10.se))/365), shape=124,size=4.5) +
  geom_segment(aes(x = Species, y = (10**T_ext.log10)/365, xend = Species, yend = (10**(T_ext.log10-T_ext.log10.se))/365), size = 1) +
  geom_segment(aes(x = Species, y = (10**T_ext.log10)/365, xend = Species, yend = (10**(T_ext.log10+T_ext.log10.se))/365), size = 1) +
  ylab(TeX("Time to extinction, $\\mathit{T_{d}}$ (years)") ) + 
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  coord_flip() +
  geom_hline(yintercept= (10**mean(df.species$T_ext.log10))/365, linetype = "longdash", color="black") +
  geom_hline(yintercept= 1000/365, linetype = "longdash", color = 'grey') +
  theme_bw() + 
  theme(axis.title.y=element_blank(), 
        axis.text.y=element_text(size = 13), 
        axis.title.x = element_text(color="black", size=24), 
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



ggsave(file="figs/T_ext.pdf", boxplot.T_ext, device='pdf', width=12,height=9, units='in', dpi=600)


#ggsave(file="figs/dead_cells_k_fig2.eps", merge_plots, device='eps', width=9,height=9, units='in', dpi=180)


