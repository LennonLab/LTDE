rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('ggplot2')
library('phylolm')
library('latex2exp')


df.metab <- read.table("data/metab_stain/staining.all.new.txt", sep="\t", 
                                 header=TRUE, stringsAsFactors = FALSE)

df.metab.noNA <- df.metab[complete.cases(df.metab), ]
df.metab.noNA.anc <- df.metab.noNA[df.metab.noNA$hist == "anc", ] 
df.metab.noNA.anc$rel.met.anc <- df.metab.noNA.anc$active / (df.metab.noNA.anc$active + df.metab.noNA.anc$dormant)
# derived
df.metab.noNA.der <- df.metab.noNA[df.metab.noNA$hist == "der", ] 
df.metab.noNA.der$active.prop.der <- df.metab.noNA.der$active / (df.metab.noNA.der$active + df.metab.noNA.der$dormant)
df.metab.noNA.der.mean <- aggregate(df.metab.noNA.der[, 7], list(df.metab.noNA.der$strain), mean)
colnames(df.metab.noNA.der.mean) <- c("strain", "rel.met.der")
# merge
df.metab.mean.merge <- merge(df.metab.noNA.anc, df.metab.noNA.der.mean, by=c("strain"))
df.metab.mean.merge$rel.met.diff <- df.metab.mean.merge$rel.met.der - df.metab.mean.merge$rel.met.anc
rownames(df.metab.mean.merge) <- df.metab.mean.merge$strain

#plot(df.metab.mean.merge$)
# Death curve

df.species <- read.table("data/demography/weibull_results_clean_species.csv", 
                         header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
rownames(df.species) <- df.species$Species
df.species<-df.species[!(df.species$Species=="KBS0727") ,]
df.species<-df.species[!(df.species$Species=="KBS0812") ,]

x <- merge(df.metab.mean.merge,df.species,by="row.names",all.x=TRUE)
rownames(x) <- x$Row.names
x <- x[!(x$strain=="KBS0711W" ),]
x <- x[!(x$strain=="KBS0812" ),]


plot(x$rel.met.diff, x$mttf.log10, xlab = "Change in proportion of active cells", ylab = "mean time to death, log10")
abline(lm(x$mttf.log10 ~ x$rel.met.diff))

plot(x$rel.met.diff, x$alpha, xlab = "Change in proportion of active cells", ylab = "shape param")
abline(lm(x$alpha ~ x$rel.met.diff))
summary(lm(x$alpha ~ x$rel.met.diff))


plot(x$rel.met.diff, x$beta.log10, xlab = "Change in proportion of active cells", ylab = "scale param")
abline(lm(x$beta.log10 ~ x$rel.met.diff))
summary(lm(x$beta.log10 ~ x$rel.met.diff))



