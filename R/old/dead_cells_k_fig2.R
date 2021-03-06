rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('ggplot2')
library('phylolm')
library('latex2exp')
library('viridis')
library('MuMIn')
library('lme4')


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

mixed.b0 <- lmer(log10(alpha) ~ log10(N_0_beta) + (1|strain), data = df)
#AIC(ols, mixed.inter)
mixed.b0_and_b1 <- lmer(log10(alpha) ~ log10(N_0_beta) + (log10(N_0_beta)|strain), data = df)
qqnorm(resid(mixed.b0_and_b1))
qqline(resid(mixed.b0_and_b1))



S_1000 <- df$N_final/df$N_0

log10(df$N_0 - df$N_final)

plot(log10(df$alpha), log10(df$N_0 - log10((df$N_0 - df$N_final) / df$N_0 ) ))

summary(lm(log10((df$N_0 - df$N_final) / df$N_0 )  ~  log10(df$alpha) ))



# test whether random slopes explain enough variation
anova(mixed.b0, mixed.b0_and_b1)
# test whether the fixed effect alpha is significant (slope test)
summary(mixed.b0_and_b1)

print(VarCorr(mixed.b0_and_b1),comp="Variance")

# get fit for plm
# Load ML tree
ml.tree <- read.tree("data/tree/RAxML_bipartitionsBranchLabels.ltde")
# Define the outgroup
outgroup <- match("NC_005042.1.353331-354795", ml.tree$tip.label)
# Create a rooted tree {ape} using the outgroup
ml.rooted <- root(ml.tree, outgroup, resolve.root = TRUE)
# Keep rooted but drop outgroup branch
ml.rooted <- drop.tip(ml.rooted, c("NC_005042.1.353331-354795"))
is.ultrametric(ml.rooted)
ml.rooted.um  <- chronos(ml.rooted)
is.ultrametric(ml.rooted.um)

df.species <- read.table("data/demography/weibull_results_clean_species.csv", 
                         header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
rownames(df.species) <- df.species$Species

phylolm.fit <- phylolm(log10(alpha) ~ N_0_beta.log10, data=df.species, phy=ml.rooted.um, model = 'lambda', boot=10)

# plot distribution of random slopes 
kde.plot <- ggplot(coef(mixed.b0_and_b1)$strain, aes(coef(mixed.b0_and_b1)$strain[,2])) +
  xlab(TeX("Random slope coefficient") ) + 
  ylab('Density') +
  geom_density(fill = "blue", alpha = 0.2) +
  theme_bw() +
  geom_vline(xintercept=0, linetype = "longdash", size=1.5, colour = 'grey')+
  geom_vline(xintercept=mixed.b0_and_b1@beta[2], linetype = "longdash", size=1.5)+
  xlim(-0.3, 0.1)


kde.plot <- kde.plot + theme(axis.title.x = element_text(color="black", size=16), 
                             axis.title.y = element_text(color="black", size=16), 
                             panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank())


r2 <- r.squaredGLMM(mixed.b0_and_b1)

merBoot <- bootMer(mixed.b0_and_b1, predict, nsim = 10000, re.form = NA)
pred <- predict(mixed.b0_and_b1, re.form = NA)

std.err <- apply(merBoot$t, 2, sd)
CI.lower <- pred - std.err*1.96
CI.upper <- pred + std.err*1.96


cell.shape.plot <- ggplot(data = df, aes(x = N_0_beta, y = alpha)) +
  geom_point(color=df$Color, alpha = 0.9, size=7) +
  xlab(TeX("Initial number of dead cells, $N(0) / \\mathit{\\lambda}$") ) + 
  ylab(TeX(" Degree of change in growth rate, $\\mathit{\\k}$")) +
  #xlab(expression(atop("A long string of text for the purpose", paste("of illustrating my point" [reported]))))

  #geom_abline(intercept = mixed.b0_and_b1@beta[1], slope = mixed.b0_and_b1@beta[2]) + 
  geom_segment(aes(x = 10**6, xend = 10**14, 
                   y = 10**(mixed.b0_and_b1@beta[1] + log10(10**6)* mixed.b0_and_b1@beta[2]), 
                   yend = 10**(mixed.b0_and_b1@beta[1] + log10(10**14)* mixed.b0_and_b1@beta[2])
                   ), linetype = 'dashed',size=2) +
  geom_segment(aes(x = 10**6, xend = 10**14, 
                   y = 10**(phylolm.fit$coefficients[1] + log10(10**6)* phylolm.fit$coefficients[2]), 
                   yend = 10**(phylolm.fit$coefficients[1] + log10(10**14)* phylolm.fit$coefficients[2])
  ), linetype = 'dashed',size=2,color='grey') +
  annotate("text", x=10**14.5, y=10**-0.3, label=TeX(sprintf("$\\mathit{r}^{2}_{m} = %g$", round(r2[1],2))), size = 6) +
  annotate("text", x=10**14.5, y=10**-0.4,
           label=TeX("$\\mathit{p}\\, < \\,10^{-5}$"), 
           size = 6) +
  

  #geom_segment(aes(x = 10**6, xend = 10**14, 
  #                 y = 10**(mixed.b0_and_b1@beta[1] + log10(10**6)* mixed.b0_and_b1.ci[6,2]), 
  #                 yend = 10**(mixed.b0_and_b1@beta[1] + log10(10**14)* mixed.b0_and_b1.ci[6,2])
  #), linetype = 'dotted',size=1, color='blue') +
  #geom_ribbon(aes(ymin = 10**lower.ci, ymax = 10**upper.ci), alpha = .1) +
  #geom_ribbon(aes(ymin = 10**CI.lower, ymax = 10**CI.upper),size=2, alpha = 1, linetype = "dashed", color="black") +
  geom_line(aes(10**merBoot$data$`log10(N_0_beta)`, 10**CI.lower), col = 'black', linetype = 'dashed', size = 1.25) +
  geom_line(aes(10**merBoot$data$`log10(N_0_beta)`, 10**CI.upper), col = 'black', linetype = 'dashed', size = 1.25) +
  
  scale_x_log10(
    limits = c(10**5, 10**16),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    limits = c(10**-1.5, 10**0.1),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  geom_hline(yintercept= 1, linetype = "longdash", size=1.5, color='grey') +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=26), 
        axis.title.y = element_text(color="black", size=26), 
        axis.text.x=element_text(size = 18), 
        axis.text.y=element_text(size = 18), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


merge_plots <- cell.shape.plot + annotation_custom(ggplotGrob(kde.plot), xmin = 4.5, xmax = 9.5, 
                       ymin = -1.6, ymax = -0.8)



ggsave(file="figs/dead_cells_k_fig2.pdf", merge_plots, device='pdf', width=9,height=9, units='in', dpi=600)

#ggsave(file="figs/dead_cells_k_fig2_no_insert.png", cell.shape.plot, width=9,height=9, units='in', dpi=600)
