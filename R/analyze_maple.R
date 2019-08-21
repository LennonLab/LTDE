rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

set.seed(123456)

library('ggplot2')
library('latex2exp')
library('ggpubr')
library('vegan')
library('ggrepel')

#source("http://www.phytools.org/utilities/v4.6/utilities.R")
#source("http://www.phytools.org/phyl.pca/v0.5/phyl.pca.R")
#require("corpcor")

df.met <- read.table("data/genomes/genomes_ncbi_maple.txt", 
                     header = TRUE, sep = "\t", row.names = 1)
df.met <- subset(df.met, select = -c( KBS0727))#, KBS0710, KBS0721))
# remove rows with all ones
df.met<- t(df.met[apply(df.met[,-1], 1, function(x) !all(x==1)),])

df.met.no0 <-  df.met[, colSums(df.met != 0) > 0]
df.met.no0or1 <-  df.met.no0[, colSums(df.met.no0 != 1) < dim(df.met)[1]]

# run regression
df.species <- read.table("data/traits/traits_weibull.txt", 
                         header = TRUE, sep = "\t", 
                         row.names = 1, stringsAsFactors = FALSE)

rownames(df.species) <- df.species$Species
rownames(df.species) <- lapply(rownames(df.species), as.character)
df.species.subset <- df.species[,c("alpha", "beta.log10","A","umax", "Lag")]
df.species.subset$Lag <- log10(df.species.subset$Lag)

# rda
m1.rda <- vegan::rda(df.met.no0or1 ~ ., as.data.frame(df.species.subset))
m0.rda <- vegan::rda(df.met.no0or1 ~ 1, as.data.frame(df.species.subset))
m.rda <- step(m0.rda, scope=formula(m1.rda), test="perm")
# significance of final model
anova.cca(m1.rda, by="term")



# rda w/out bacillus
df.met.no0or1.noKBS0812 <- df.met.no0or1[!rownames(df.met.no0or1) %in% c('KBS0812'), ]
df.species.subset.noKBS0812 <- df.species.subset[!rownames(df.species.subset) %in% c('KBS0812'), ]

m1.rda.noKBS0812 <- vegan::rda(df.met.no0or1.noKBS0812 ~ ., as.data.frame(df.species.subset.noKBS0812))
m0.rda.noKBS0812 <- vegan::rda(df.met.no0or1.noKBS0812 ~ 1, as.data.frame(df.species.subset.noKBS0812))
m.rda.noKBS0812 <- step(m0.rda.noKBS0812, scope=formula(m1.rda.noKBS0812), test="perm")
# reduced model
m.rda.reduced.noKBS0812 <- vegan::rda(df.met.no0or1.noKBS0812 ~ umax, as.data.frame(df.species.subset.noKBS0812))

anova.cca(m.rda.reduced.noKBS0812, by="term", permutations = how(nperm=9999))



# make plot
# the code for this plot was adapted from the GitHub repo HJA-streams
# under a GNU General Public License v3.0
# https://github.com/nwisnoski/HJA-streams
rda.plot <- cbind.data.frame(scores(m1.rda.noKBS0812)$sites)
rda.var1 <- round(eigenvals(m1.rda.noKBS0812)[1] / sum(eigenvals(m1.rda.noKBS0812)) * 100, 2)
rda.var2 <- round(eigenvals(m1.rda.noKBS0812)[2] / sum(eigenvals(m1.rda.noKBS0812)) * 100, 2)
rda.vecs <- as.data.frame(m1.rda.noKBS0812$CCA$biplot)
rda.vecs$predictor <- c("alpha", "beta.log10", "A", "umax", "Lag")
rda.vecs$origin <- 0
scale.arrows = 4

ggplot(data = rda.plot, aes(RDA1, RDA2)) +
  geom_hline(aes(yintercept = 0), color = 'black', linetype = "dashed", alpha = 1, size = 0.25) +
  geom_vline(aes(xintercept = 0),color = 'black',  linetype = "dashed", alpha = 1, size = 0.25) +
  geom_point(color='blue', alpha = 0.8, size=5)+
  #geom_point(data = subset(rda.plot, group == "Sediment"), shape = 1, color = "black", size = 2) +
  #geom_point(data = subset(rda.plot, group == "Water"), shape = 2, color = "black", size = 2) +
  # geom_path(data = df_ell, 
  #           aes(x = Dim1, y = Dim2, color = group),
  #           size = .8, alpha = 1, linetype = 2) +
  labs(x = paste0("RDA1 (", rda.var1, "%)"),
       y = paste0("RDA2 (", rda.var2, "%)"),
       color = "Habitat", shape = "Habitat") +
  #coord_fixed() +
  #coord_fixed()
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_segment(data = rda.vecs, size = .5,
               aes(x = origin, y = origin, 
                   xend = scale.arrows*RDA1, 
                   yend = scale.arrows*RDA2),
               alpha = 1, color = "black",
               arrow = arrow(angle = 20,
                             length = unit(.1, "inches"),
                             type = "open")) +
  annotate("text", x=(rda.vecs[4,1]*scale.arrows) + 0.3, y=(rda.vecs[4,2]*scale.arrows-0.5), label=TeX("$\\mu_{max}^{*}$"), size = 4) +
  annotate("text", x=(rda.vecs[3,1]*scale.arrows) - 0.2, y=(rda.vecs[3,2]*scale.arrows)+0.6, label=TeX("$Yield$"), size = 4) +
  annotate("text", x=(rda.vecs[2,1]*scale.arrows) + 0.1, y=(rda.vecs[2,2]*scale.arrows)-0.3, label=TeX("$\\log_{10} (\\lambda) $"), size = 4) +
  annotate("text", x=(rda.vecs[1,1]*scale.arrows) + 0.2, y=rda.vecs[1,2]*scale.arrows, label=TeX("$k$"), size = 4) +
  annotate("text", x=(rda.vecs[5,1]*scale.arrows) - 0.3, y=rda.vecs[5,2]*scale.arrows, label=TeX("$\\log_{10} (Lag \\, time) $"), size = 4) +
  ggsave("figs/RDA.png", width = 84, height = 84, units = "mm", dpi = 500)


