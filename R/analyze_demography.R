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

set.seed(123456789)

df <- read.csv("data/demography/weibull_results_clean.csv", 
               header = TRUE, stringsAsFactors = FALSE)
df.species <- read.table("data/demography/weibull_results_clean_species.csv", 
                 header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
rownames(df.species) <- df.species$Species
df.species<-df.species[!(df.species$Species=="KBS0727"),]


###### ggplot KDE
bw <- bw.CV(df.species$mttf.log10, method="LCV", lower=0, upper=100)
# mean of means
10**mean(df.species$mttf.log10)
# anova analysis
aov <- aov(log10(mttf) ~ strain, data = df)
summary(aov)
df_numerator <- nrow(df.species) - 1
df_denominator <- nrow(df) - nrow(df.species)
# 20, 73

kde.plot <- ggplot(df.species, aes(10** mttf.log10)) +
            xlab(TeX("Mean time to death, $\\bar{T}_{d}$ (days)") ) + 
            ylab('Density') +
            geom_density(fill = "blue", alpha = 0.2) +
            theme_bw() +
            scale_x_log10(
              limits = c(0.1, 5000),
              breaks = scales::trans_breaks("log10", function(x) 10^x),
              labels = scales::trans_format("log10", scales::math_format(10^.x))
            ) +
            geom_vline(xintercept= 10**mean(df.species$mttf.log10), linetype = "longdash") 
  

kde.plot <- kde.plot + theme(axis.title.x = element_text(color="black", size=10), 
                            axis.title.y = element_text(color="black", size=14), 
                            panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank())


# reorder rows so that beta is in increasing order
df.species.order <- df.species[rev(order(df.species$beta.log10)),]
df.species.order.beta <- df.species.order$beta.log10
names(df.species.order.beta) <- df.species.order$Species
df.species.order.alpha <- df.species.order$alpha
names(df.species.order.alpha) <- df.species.order$Species

chull.idx <- chull(df.species.order.beta, df.species.order.alpha)
# remove second to last object in list
chull.idx <- chull.idx[-c(6)]
chull.beta <- df.species.order.beta[chull.idx]
chull.alpha <- df.species.order.alpha[chull.idx]

# complete null test
Fxn <- approxfun(x=chull.beta, chull.alpha)
null.pareto.test <- function(betas, alphas){
  num.above.line.sims <- c()
  for (i in seq(1, 10000)){
    beta.rndm <- sample(betas)
    alpha.rndm <- sample(alphas)
    num.above.line <- 0
    for (j in seq(1, length(beta.rndm))){

      # the lower the alpha, the higher the mean time to death
      # so to test for the Pareto front, we want to see if the 
      # simlated alpha is less than the alpha we get from Fxn
      if (alpha.rndm[j] < Fxn(beta.rndm[j]) ){
        num.above.line <- num.above.line + 1
      }
    }
    num.above.line.sims <- c(num.above.line.sims, num.above.line)
  }
  return(length(num.above.line.sims[num.above.line.sims==0]) / length(num.above.line.sims))
}

p.value <- null.pareto.test(df.species.order.beta, df.species.order.alpha)



phylo.params <- ggplot(data = df.species, aes(x = 10**(beta.log10 ), y = alpha)) +
                geom_point(color='blue', alpha = 0.6, size=4) +
                xlab(TeX("Scale paramater, $\\lambda$") ) + 
                ylab(TeX("Shape paramater, $\\k$")) +
                scale_y_continuous(limits = c(0, 1.05)) +
                scale_x_log10(
                  limits = c(0.000001, 200),
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))
                ) +

                geom_hline(yintercept= 1, linetype = "longdash") +
                theme_bw() +
                geom_line(aes(y = y, x = x), size=0.75, data=data.frame(x=10**chull.beta, y=chull.alpha))
                #geom_segment(aes(x=min(chull.beta), xend=max(chull.beta),y=min(chull.alpha),yend=min(chull.alpha)))


phylo.params <- phylo.params + 
                theme(axis.title.x = element_text(color="black", size=10), 
                            axis.title.y = element_text(color="black", size=10), 
                            panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank()) + 
                scale_y_reverse() +
                geom_line(aes(y = y, x = x), size=0.75, 
                          data=data.frame(x=10**min(chull.beta):10**max(chull.beta), 
                                          y=min(chull.alpha)), linetype = 'dotted') +
                geom_line(aes(y = y, x = x), size=0.75, 
                          data=data.frame(y=seq(from = min(chull.alpha), to = 1, length.out = 10), 
                                          x=10**max(chull.beta)), linetype = 'dotted') +
                geom_point(aes(x=10**max(chull.beta), y=min(chull.alpha)), colour="black")





#### make boxplot ggplot
boxplot <- ggplot(data = df.species) +
  geom_point(aes(x = reorder(Species, -mttf.log10), y = 10**mttf.log10), color='blue', alpha = 0.6, size=2.2) +
  geom_point(aes(x = reorder(Species, -mttf.log10), y = 10**(mttf.log10-mttf.log10.se)), shape=124,size=2.5 ) +
  geom_point(aes(x = reorder(Species, -mttf.log10), y = 10**(mttf.log10+mttf.log10.se)), shape=124,size=2.5 ) +
  geom_segment(aes(x = Species, y = 10**mttf.log10, xend = Species, yend = 10**(mttf.log10-mttf.log10.se)), size = 0.5) +
  geom_segment(aes(x = Species, y = 10**mttf.log10, xend = Species, yend = 10**(mttf.log10+mttf.log10.se)), size = 0.5) +
  ylab(TeX("Mean time to death, $\\bar{T}_{d}$ (days)") ) + 
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  coord_flip() +
  theme_bw() + 
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
  


g <- ggarrange(boxplot,                                                 # First row with scatter plot
          ggarrange(kde.plot, phylo.params, ncol = 2, labels = c("b", "c")), # Second row with box and dot plots
          nrow = 2, 
          labels = "auto")


ggsave(file="figs/demography.png", g, units='in', dpi=600)

