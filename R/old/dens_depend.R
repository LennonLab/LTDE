#### trait regression
rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('phylolm')
library('ggplot2')
library("ape")
library('latex2exp')

df <- read.table("data/demography/weibull_results_clean_species.csv", 
                 header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
df<-df[!(df$Species=="KBS0727"),]# | df$Species=="KBS0812"),]
#df<-df[!(df$Species=="KBS0727"),]
rownames(df) <- df$Species
df$N_diff <- log10(df$N_0 - df$N_final)

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


phylo.beta.log10 <- phylolm(beta.log10 ~ N_diff, data = df, 
                       ml.rooted.um, model = 'lambda', boot = 0)
phylo.alpha <- phylolm(alpha ~ N_diff, data = df, 
                      ml.rooted.um, model = 'lambda', boot = 0)

summary(phylo.beta.log10)
summary(phylo.alpha)


beta.log10 <- ggplot(data = df, aes(x = 10**N_diff, y = 10**beta.log10)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Mean scale paramater, $\\bar{\\lambda}$") ) + 
  xlab(TeX("Mean change in cell number, $\\bar{\\Delta N}$")) +
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


alpha <- ggplot(data = df, aes(x = 10**N_diff, y = alpha)) +
  geom_point(color='blue', alpha = 0.6, size=4) +
  ylab(TeX("Mean shape paramater, $\\bar{k}$")) +
  xlab(TeX("Mean change in cell number, $\\bar{\\Delta N}$")) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  #geom_line(aes(y = y, x = x), size=0.75, data=data.frame(x=10**x.line, y=y.line)) +
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


g <- ggarrange(beta.log10, alpha,
                 # First row with scatter plot
                 # Second row with box and dot plots
                 ncol = 2, nrow = 1,
                 labels = "auto")#, label.y = c(1, 0.5, 0.25)    


ggsave(file="figs/dens_depend.png", g, width=10,height=5, units='in', dpi=600)


