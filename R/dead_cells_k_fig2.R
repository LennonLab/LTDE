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

library('lattice')

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


r2 <- r.squaredGLMM(mixed.b0_and_b1)

lower.ci <- mixed.b0_and_b1@beta[1] + mixed.b0_and_b1.ci[6,1] * log10(df$N_0_beta)
upper.ci <- mixed.b0_and_b1@beta[1] + mixed.b0_and_b1.ci[6,2] * log10(df$N_0_beta)

#pred <- cbind(df, predictInterval(mixed.b0_and_b1, df, level=0.95))


merBoot <- bootMer(mixed.b0_and_b1, predict, nsim = 10000, re.form = NA)
pred <- predict(mixed.b0_and_b1, re.form = NA)

std.err <- apply(merBoot$t, 2, sd)
CI.lower <- pred - std.err*1.96
CI.upper <- pred + std.err*1.96


cell.shape.plot <- ggplot(data = df, aes(x = N_0_beta, y = alpha)) +
  geom_point(color=df$Color, alpha = 0.9, size=6) +
  xlab(TeX("Initial number of dead cells, $N(0) / \\mathit{\\lambda}$") ) + 
  ylab(TeX("Shape paramater, $\\mathit{\\k}$")) +
  #geom_abline(intercept = mixed.b0_and_b1@beta[1], slope = mixed.b0_and_b1@beta[2]) + 
  geom_segment(aes(x = 10**6, xend = 10**14, 
                   y = 10**(mixed.b0_and_b1@beta[1] + log10(10**6)* mixed.b0_and_b1@beta[2]), 
                   yend = 10**(mixed.b0_and_b1@beta[1] + log10(10**14)* mixed.b0_and_b1@beta[2])
  ), linetype = 'dashed',size=2) +
  annotate("text", x=10**6, y=10**-1.3, label=TeX(sprintf("$\\mathit{r}^{2}_{m} = %g$", round(r2[1],2))), size = 6) +
  #annotate("text", x=10**14, y=10**0.3,
  #         label=TeX(sprintf("$\\mathit{p} = %e$", 
  #                           round(summary(mixed.b0_and_b1)$coefficients[2,5],9))), 
  #         size = 6) +
  annotate("text", x=10**6, y=10**-1.4,
           label=TeX("$\\mathit{p}\\, < \\,10^{-5}$"), 
           size = 6) +
  
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
    limits = c(10**5, 10**16),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    limits = c(10**-1.5, 10**0.2),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  geom_hline(yintercept= 1, linetype = "longdash", size=1.5, color='grey') +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=22), 
        axis.title.y = element_text(color="black", size=22), 
        axis.text.x=element_text(size = 16), 
        axis.text.y=element_text(size = 16), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


ggsave(file="figs/dead_cells_k_fig2.png", cell.shape.plot, width=9,height=9, units='in', dpi=600)



# make and save simple histogram of slopes? 
#dotplot(ranef(mixed.b0_and_b1,condVar=TRUE))

#dotplot(ranef(fit, condVar=TRUE))


