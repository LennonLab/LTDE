#### trait regression
rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('phylolm')
library('ggplot2')
library("ape")
library('latex2exp')
library('ggpubr')

df.species <- read.table("data/demography/weibull_results_clean_species.csv", 
                         header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
rownames(df.species) <- df.species$Species
traits <-  read.table("data/traits/traits.txt", 
                      header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

traits.merge <- merge(df.species, traits, by="row.names")
rownames(traits.merge) <- traits.merge$Species
# save table
traits.merge <- traits.merge[, !(colnames(traits.merge) %in% c("Row.names"))]

traits.irep.merge <- traits.merge[, !(colnames(traits.merge) %in% c("Species.x", "Row.names", "Species.y"))]
write.table(traits.irep.merge, "data/traits/traits_weibull.txt", sep="\t")

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



fit.trait.beta.select <- phylostep(beta.log10 ~ A + log10(umax) + log10(Lag), starting.formula=NULL, 
                                   data= traits.merge, phy=ml.rooted.um, model = 'lambda', 
                                   lower.bound = 0, upper.bound = 1)
fit.trait.alpha.select <- phylostep(alpha ~ A + log10(umax) + log10(Lag), starting.formula=NULL, 
                                    data= traits.merge, phy=ml.rooted.um, model = 'lambda', 
                                    lower.bound = 0, upper.bound = 1)

summary(fit.trait.beta.select)
summary(fit.trait.alpha.select)



phy.beta.umax.plot <- ggplot(data = traits.merge, aes(x = umax, y = 10** beta.log10)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Mean scale paramater, $\\bar{\\lambda}$") ) + 
  xlab(TeX("Max. growth rate ($hours^{-1}$), $\\mu_{max}$") ) + 
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme(axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())#,
        #axis.text.x=element_blank(),
        #axis.title.x=element_blank())


phy.beta.yield.plot <- ggplot(data = traits.merge, aes(x = A, y = 10**beta.log10)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Mean scale paramater, $\\bar{\\lambda}$") ) +  
  xlab(TeX("Yield, $OD_{600}$")) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())#,
        #axis.text.x=element_blank(),
        #axis.title.x=element_blank())



phy.beta.lag.plot <- ggplot(data = traits.merge, aes(x = Lag, y = 10**beta.log10)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Mean scale paramater, $\\bar{\\lambda}$") ) + 
  xlab(TeX("Lag time (hours)")) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



phy.alpha.umax.plot <- ggplot(data = traits.merge, aes(x = umax, y = alpha)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Mean shape paramater, $\\bar{k}$")) +
  xlab(TeX("Max. growth rate ($hours^{-1}$), $\\mu_{max}$") ) + 
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# add slope
phylo.alpha.yield <- phylolm(alpha ~ A, data= traits.merge, phy=ml.rooted.um, model = 'lambda')
# add slope
phylo.alpha.yield.x.line <- seq(0.5, 3.5, length.out = 1000) 
phylo.alpha.yield.y.line <- phylo.alpha.yield$coefficients[1] + phylo.alpha.yield$coefficients[2] * phylo.alpha.yield.x.line
phy.alpha.yield.plot <- ggplot(data = traits.merge, aes(x = A, y = alpha)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Mean shape paramater, $\\bar{k}$")) +
  xlab(TeX("Yield, $OD_{600}$")) +
  scale_y_continuous(limits = c(0, 1)) +
  #stat_function(fun = function(x) fit.trait.alpha.select$coefficients[1] + fit.trait.alpha.select$coefficients[2] * x) + 
  geom_line(aes(y = y, x = x), size=0.75, data=data.frame(x=phylo.alpha.yield.x.line, y=phylo.alpha.yield.y.line)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


phylo.alpha.lag <- phylolm(alpha ~ log10(Lag), data= traits.merge, phy=ml.rooted.um, model = 'lambda')
# add slope
phylo.alpha.lag.x.line <- seq(0, 1.5, length.out = 1000) 
phylo.alpha.lag.y.line <- phylo.alpha.lag$coefficients[1] + phylo.alpha.lag$coefficients[2] * phylo.alpha.lag.x.line
phy.alpha.lag.plot <- ggplot(data = traits.merge, aes(x = Lag, y = alpha)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Mean shape paramater, $\\bar{k}$")) +
  xlab(TeX("Lag time (hours)")) +
  #stat_function(fun = function(x) fit.trait.alpha.select$coefficients[1] + fit.trait.alpha.select$coefficients[4] * x) + 
  #geom_line(aes(y = y, x = x), size=0.75, data=data.frame(x=10**phylo.alpha.lag.x.line, y=phylo.alpha.lag.y.line)) +

  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  stat_smooth(method = "lm") +
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


g <- ggarrange(phy.beta.umax.plot, phy.beta.yield.plot, phy.beta.lag.plot,
               phy.alpha.umax.plot, phy.alpha.yield.plot, phy.alpha.lag.plot,
               # First row with scatter plot
               # Second row with box and dot plots
               ncol = 3, nrow = 2,
               labels = "auto")#, label.y = c(1, 0.5, 0.25)    


ggsave(file="figs/traits.png", g, width=15,height=10, units='in', dpi=600)





# main figure
phy.beta.lag.plot.main <- ggplot(data = traits.merge, aes(x =Lag, y = 10**beta.log10.se)) +
  geom_point(color='blue', alpha = 0.6, size=6) +
  ylab(TeX("Mean scale paramater, $\\bar{\\lambda}$") ) + 
  xlab(TeX("Lag time (hours)")) +
  #stat_function(fun = function(x) fit.trait.alpha.select$coefficients[1] + fit.trait.alpha.select$coefficients[4] * x) + 
  #geom_line(aes(y = y, x = x), size=0.75, data=data.frame(x=10**phylo.alpha.lag.x.line, y=phylo.alpha.lag.y.line)) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=23), 
        axis.title.y = element_text(color="black", size=23), 
        axis.text.x=element_text(size = 14),
        axis.text.y=element_text(size = 14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())




phylo.alpha.lag <- phylolm(alpha ~ log10(Lag), data= traits.merge, phy=ml.rooted.um, model = 'lambda')
summary(phylo.alpha.lag)
lag.shape.summary <- summary(lm(alpha~log10(Lag), traits.merge))
r2 <-lag.shape.summary$r.squared
p.value <- lag.shape.summary$coefficients[8]
phy.alpha.lag.plot.main <- ggplot(data = traits.merge, aes(x = Lag, y = alpha)) +
  geom_point(color='blue', alpha = 0.6, size=6) +
  ylab(TeX("Mean shape paramater, $\\bar{k}$")) +
  xlab(TeX("Lag time (hours)")) +
  #stat_function(fun = function(x) fit.trait.alpha.select$coefficients[1] + fit.trait.alpha.select$coefficients[4] * x) + 
  #geom_line(aes(y = y, x = x), size=0.75, data=data.frame(x=10**phylo.alpha.lag.x.line, y=phylo.alpha.lag.y.line)) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  stat_smooth(method = "lm", color='black', size=2, linetype = "dashed") +
  #scale_y_continuous(limits = c(0, 1)) +
  annotate("text", x=0.82, y=0.9, label=TeX(sprintf("$r^{2} = %g$", round(r2,2))), size = 6) +
  annotate("text", x=0.99, y=0.82, label=TeX(sprintf("$p = %g$", round(p.value,4))), size = 6) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=23), 
        axis.title.y = element_text(color="black", size=23), 
        axis.text.x=element_text(size = 14),
        axis.text.y=element_text(size = 14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


g <- ggarrange(phy.beta.lag.plot.main, phy.alpha.lag.plot.main,                                               # First row with scatter plot
               nrow = 1, ncol=2, 
               labels = "auto")
ggsave(file="figs/shape_scale_lag.png", g, width=10,height=5, units='in', dpi=600)

