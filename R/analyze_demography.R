rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('ggplot2')
library('phylolm')
library('NPCirc')
library('latex2exp')
library('rr2')
library('dplyr')
library('gridExtra')
library('ggpubr')
library('grid')

library('scales')

library('viridis')
library('lmerTest')

library('MuMIn')

library('merTools')
library('lme4')


library('nlme')

df <- read.csv("data/demography/weibull_results_clean.csv", 
               header = TRUE, stringsAsFactors = FALSE)
df.species <- read.table("data/demography/weibull_results_clean_species.csv", 
                 header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
rownames(df.species) <- df.species$Species
df.species <- df.species[!(df.species$Species=="KBS0727"),]
df.species <- df.species[order(df.species$mttf.log10),]
# get colors

colors <- viridis(dim(df.species)[1], option = 'viridis')
#df.species$cols <- colors

col.df <- data.frame("Color" =colors)
col.df$strain <- df.species$Species
rownames(col.df) <- df.species$Species

df.species <- merge(df.species, col.df, by="row.names",all.x=TRUE)

df <- merge(df, col.df, by="strain",all.x=TRUE)

#### make boxplot ggplot
boxplot <- ggplot(data = df.species) +
  geom_point(aes(x = reorder(Species, -mttf.log10), y = 10**mttf.log10), color=df.species$Color, alpha = 1, size=2.2) +
  geom_point(aes(x = reorder(Species, -mttf.log10), y = 10**(mttf.log10-mttf.log10.se)), shape=124,size=2.5 ) +
  geom_point(aes(x = reorder(Species, -mttf.log10), y = 10**(mttf.log10+mttf.log10.se)), shape=124,size=2.5 ) +
  geom_segment(aes(x = Species, y = 10**mttf.log10, xend = Species, yend = 10**(mttf.log10-mttf.log10.se)), size = 0.5) +
  geom_segment(aes(x = Species, y = 10**mttf.log10, xend = Species, yend = 10**(mttf.log10+mttf.log10.se)), size = 0.5) +
  ylab(TeX("Mean time to death, $\\mathit{\\bar{T}_{d}}$ (days)") ) + 
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  coord_flip() +
  geom_hline(yintercept= 10**mean(df.species$mttf.log10), linetype = "longdash") +
  theme_bw() + 
  annotate("text", x=19, y=1, label= "a", size = 5, fontface =2) +
  theme(axis.title.y=element_blank(), 
        axis.text.y=element_text(size = 7), 
        axis.title.x = element_text(color="black", size=14), 
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

#### examine relationship between number of dead cells and shape paraemter
df$N_0_beta <- df$N_0/df$beta
df$N_0_beta_log10 <- log10(df$N_0_beta)

#ols <- lm(log10(alpha) ~ log10(N_0_beta), data = df)
mixed.b0 <- lmer(log10(alpha) ~ log10(N_0_beta) + (1|strain), data = df)
#AIC(ols, mixed.inter)
mixed.b0_and_b1 <- lmer(log10(alpha) ~ log10(N_0_beta) + (log10(N_0_beta)|strain), data = df)
# test whether random slopes explain enough variation
anova(mixed.b0, mixed.b0_and_b1)
# test whether the fixed effect alpha is significant (slope test)
summary(mixed.b0_and_b1)

print(VarCorr(mixed.b0_and_b1),comp="Variance")

coef(mixed.b0_and_b1)$strain

mixed.b0_and_b1.ci <- confint(mixed.b0_and_b1)
#bands <- data.frame(cbind(mixed.b0_and_b1.ci[6,1], mixed.b0_and_b1.ci[6,2]))

#lower.ci <-  mixed.b0_and_b1@beta[1] * 1:max(x)

r.squaredGLMM(mixed.b0_and_b1)

lower.ci <- mixed.b0_and_b1@beta[1] + mixed.b0_and_b1.ci[6,1] * log10(df$N_0_beta)
upper.ci <- mixed.b0_and_b1@beta[1] + mixed.b0_and_b1.ci[6,2] * log10(df$N_0_beta)

#pred <- cbind(df, predictInterval(mixed.b0_and_b1, df, level=0.95))


merBoot <- bootMer(mixed.b0_and_b1, predict, nsim = 1000, re.form = NA)
pred <- predict(mixed.b0_and_b1,re.form = NA)

std.err <- apply(merBoot$t, 2, sd)
CI.lower <- pred - std.err*1.96
CI.upper <- pred + std.err*1.96


cell.shape.plot <- ggplot(data = df, aes(x = N_0_beta, y = alpha)) +
  geom_point(color=df$Color, alpha = 0.9, size=4) +
  xlab(TeX("Initial number of dead cells, $\\mathit{N}(0) / \\mathit{\\lambda}$") ) + 
  ylab(TeX("Shape paramater, $\\mathit{\\k}$")) +
  #geom_abline(intercept = mixed.b0_and_b1@beta[1], slope = mixed.b0_and_b1@beta[2]) + 
  geom_segment(aes(x = 10**6, xend = 10**14, 
                   y = 10**(mixed.b0_and_b1@beta[1] + log10(10**6)* mixed.b0_and_b1@beta[2]), 
                   yend = 10**(mixed.b0_and_b1@beta[1] + log10(10**14)* mixed.b0_and_b1@beta[2])
                   ), linetype = 'dashed',size=1) +
  
  #geom_segment(aes(x = 10**6, xend = 10**14, 
  #                 y = 10**(mixed.b0_and_b1@beta[1] + log10(10**6)* mixed.b0_and_b1.ci[6,1]), 
  #                 yend = 10**(mixed.b0_and_b1@beta[1] + log10(10**14)* mixed.b0_and_b1.ci[6,1])
  #), linetype = 'dotted',size=1, color='blue') +
  #
  #geom_segment(aes(x = 10**6, xend = 10**14, 
  #                 y = 10**(mixed.b0_and_b1@beta[1] + log10(10**6)* mixed.b0_and_b1.ci[6,2]), 
  #                 yend = 10**(mixed.b0_and_b1@beta[1] + log10(10**14)* mixed.b0_and_b1.ci[6,2])
  #), linetype = 'dotted',size=1, color='blue') +
  
  #geom_ribbon(aes(ymin = 10**lower.ci, ymax = 10**upper.ci), alpha = .1) +
  
  geom_ribbon(aes(ymin = 10**CI.lower, ymax = 10**CI.upper),size=2, alpha = 0.2, linetype = 0) +

  scale_x_log10(
    #limits = c(0.000001, 200),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    limits = c(10**-2, 10**0.5),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  geom_hline(yintercept= 1, linetype = "longdash") +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=10), 
        axis.title.y = element_text(color="black", size=10), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())







ggsave(file="figs/demography.png", cell.shape.plot,  width=4.2,height=4.9, units='in', dpi=600)












g <- ggarrange(boxplot,                                                 # First row with scatter plot
          ggarrange(kde.plot, phylo.params, ncol = 2, labels = c("b", "c")), # Second row with box and dot plots
          nrow = 2, 
          labels = NULL)







############ old code
