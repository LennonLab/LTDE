#### trait regression
rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('phylolm')
library('ggplot2')
library("ape")
library('latex2exp')
library('ggpubr')
library('plotrix')

plot(c(1,2,3,4), c(3,4,5,6))

df.species <- read.table("data/demography/weibull_results_clean_species.csv", 
                         header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
rownames(df.species) <- df.species$Species
div <-  read.table("data/breseq/genetic_diversity_taxa.txt",sep = "\t",header = TRUE)
df.merge <- merge(df.species, div,by="Species")

# two sample t-test for beta and alpha b/w taxa with enough and not enough mutations
rownames(df.merge) <- df.merge$Species
nonsyn.f <- df.merge[c("KBS0802", "KBS0721", "KBS0812", "KBS0801", "KBS0707", "KBS0713", "KBS0715"),]
nonsyn.t <- df.merge[c("ATCC13985", "KBS0702", "KBS0711", "KBS0712"),]

t.test(nonsyn.t$beta.log10, nonsyn.f$beta.log10, alternative = "greater", var.equal = FALSE)
t.test(nonsyn.t$alpha, nonsyn.f$alpha, alternative = "greater", var.equal = FALSE)


summary(lm(df.merge$beta.log10 ~ df.merge$mean_freq))
summary(lm(df.merge$alpha ~ df.merge$mean_freq))

summary(lm(df.merge$beta.log10 ~ df.merge$Theta))
summary(lm(df.merge$alpha ~ df.merge$Theta))

summary(lm(df.merge$beta.log10 ~ log10(df.merge$Pi)))
summary(lm(df.merge$alpha ~ log10(df.merge$Pi)))

summary(lm(df.merge$beta.log10 ~ log10(df.merge$Tajimas_D)))
summary(lm(df.merge$alpha ~ log10(df.merge$Tajimas_D)))


f.beta.plot <- ggplot(data = df.merge, aes(x = mean_freq, y = 10** beta.log10)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  xlab(TeX("Mean mutation frequency, $\\f$")) +
  ylab(TeX("Scale paramater, $\\lambda$") ) + 
  #scale_y_continuous(limits = c(0, 1)) +
  scale_y_log10(
    limits = c(0.000001, 2000),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=11), 
        axis.title.y = element_text(color="black", size=12), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) #+ scale_y_reverse()


f.alpha.plot <- ggplot(data = df.merge, aes(x = mean_freq, y = alpha)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Shape paramater, $\\k$")) +
  xlab(TeX("Mean mutation frequency, $\\f$")) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=11), 
        axis.title.y = element_text(color="black", size=12), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())




theta.beta.plot <- ggplot(data = df.merge, aes(x = Theta, y = 10** beta.log10)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  xlab(TeX("Number of mutated sites, $\\theta_{W} \\, (\\mathrm{bp}^{-1})$")) +
  ylab(TeX("Scale paramater, $\\lambda$") ) + 
  #scale_y_continuous(limits = c(0, 1)) +
  scale_y_log10(
    limits = c(0.000001, 2000),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=11), 
        axis.title.y = element_text(color="black", size=12), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) #+ scale_y_reverse()


theta.alpha.plot <- ggplot(data = df.merge, aes(x = Theta, y = alpha)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Shape paramater, $\\k$")) +
  xlab(TeX("Number of mutated sites, $\\theta_{W} \\, (\\mathrm{bp}^{-1})$")) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=11), 
        axis.title.y = element_text(color="black", size=12), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



pi.beta.plot <- ggplot(data = df.merge, aes(x = Pi, y = 10** beta.log10)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  xlab(TeX("Nucleotide diversity, $\\theta_{\\pi} \\, (\\mathrm{bp}^{-1})$")) +
  ylab(TeX("Scale paramater, $\\lambda$") ) + 
  #scale_y_continuous(limits = c(0, 1)) +
  scale_x_log10(
    limits = c(0.0000001, 0.0001),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    limits = c(0.000001, 2000),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=11), 
        axis.title.y = element_text(color="black", size=12), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) #+ scale_y_reverse()


pi.alpha.plot <- ggplot(data = df.merge, aes(x = Pi, y = alpha)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Shape paramater, $\\k$")) +
  xlab(TeX("Nucleotide diversity, $\\theta_{\\pi} \\, (\\mathrm{bp}^{-1})$")) +
  scale_x_log10(
    limits = c(0.0000001, 0.0001),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=11), 
        axis.title.y = element_text(color="black", size=12), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



td.beta.plot <- ggplot(data = df.merge, aes(x = Tajimas_D, y = 10** beta.log10)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  xlab(TeX("Difference between $\\theta_{\\pi}$ and  $\\theta_{W}$, $T_{D}$")) +
  ylab(TeX("Scale paramater, $\\lambda$") ) + 
  #scale_y_continuous(limits = c(0, 1)) +
  scale_x_log10(
    limits = c(0.0000005, 0.001),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    limits = c(0.000001, 2000),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=11), 
        axis.title.y = element_text(color="black", size=12), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) #+ scale_y_reverse()


td.alpha.plot <- ggplot(data = df.merge, aes(x = Tajimas_D, y = alpha)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Shape paramater, $\\k$")) +
  xlab(TeX("Difference between $\\theta_{\\pi}$ and  $\\theta_{W}$, $T_{D}$")) +
  scale_x_log10(
    limits = c(0.0000005, 0.001),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=11), 
        axis.title.y = element_text(color="black", size=12), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



g <- ggarrange(f.beta.plot, theta.beta.plot, pi.beta.plot, td.beta.plot, 
               f.alpha.plot, theta.alpha.plot, pi.alpha.plot, td.alpha.plot,     
               # Second row with box and dot plots
               ncol = 4, nrow = 2,
               labels = "auto")#, label.y = c(1, 0.5, 0.25)                                     # Labels of the scatter plot

ggsave(file="figs/genetic_diversity.png", g, width=12,height=8, units='in', dpi=600)

