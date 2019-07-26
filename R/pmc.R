rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

# unblock to install package
#library("devtools")
#install_github("cboettig/pmc")

library("pmc")
library("ape")


df <- read.table("data/demography/weibull_results_clean_species.csv", 
                 header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
df<-df[!(df$Species=="KBS0727" ),]
df.no812 <- df[!(df$Species=="KBS0812" ),]
rownames(df) <- df$Species
rownames(df.no812) <- df.no812$Species


traits <-  read.table("data/traits/traits.txt", 
                      header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
traits$Species <- rownames(traits)
traits.merge <- merge(df, traits, by="row.names")
rownames(traits.merge) <- traits.merge$Species

traits.no812 <- traits[!(traits$Species=="KBS0812" ),]
traits.no812.merge <- merge(df.no812, traits.no812, by="row.names")
rownames(traits.no812.merge) <- traits.no812.merge$Species

# Load ML tree
ml.tree <- read.tree("data/tree/RAxML_bipartitionsBranchLabels.ltde")
# Define the outgroup
outgroup <- match("NC_005042.1.353331-354795", ml.tree$tip.label)
# Create a rooted tree {ape} using the outgroup
ml.rooted <- root(ml.tree, outgroup, resolve.root = TRUE)
# Keep rooted but drop outgroup branch
ml.rooted <- drop.tip(ml.rooted, c("NC_005042.1.353331-354795"))
ml.rooted.no812 <- drop.tip(ml.rooted, c("KBS0812"))

is.ultrametric(ml.rooted)
ml.rooted.um  <- chronos(ml.rooted)
is.ultrametric(ml.rooted.um)

is.ultrametric(ml.rooted.no812)
ml.rooted.no812.um  <- chronos(ml.rooted.no812)
is.ultrametric(ml.rooted.no812.um)

beta.log10 <- df$beta.log10
names(beta.log10) <- df$Species

alpha <- df$alpha
names(alpha) <- df$Species

umax <- traits.merge$umax
names(umax) <- traits.merge$Species

lag <- log10(traits.merge$Lag)
names(lag) <- traits.merge$Species

yield <- traits.merge$A
names(yield) <- traits.merge$Species

# no KBS0812
beta.log10.no812 <- df.no812$beta.log10
names(beta.log10.no812) <- df.no812$Species

alpha.no812 <- df.no812$alpha
names(alpha.no812) <- df.no812$Species

umax.no812 <- traits.no812.merge$umax
names(umax.no812) <- traits.no812.merge$Species.x

lag.no812 <- log10(traits.no812.merge$Lag)
names(lag.no812) <- traits.no812.merge$Species.x

yield.no812 <- traits.no812.merge$A
names(yield.no812) <- traits.no812.merge$Species.x


iter <- 1000
BM.PL.beta.log10 <- pmc(ml.rooted.um, beta.log10, "BM", "lambda", nboot = iter)
BM.PL.alpha <- pmc(ml.rooted.um, alpha, "BM", "lambda", nboot = iter )
BM.PL.umax <- pmc(ml.rooted.um, umax, "BM", "lambda", nboot = iter)
BM.PL.lag <- pmc(ml.rooted.um, lag, "BM", "lambda", nboot = iter)
BM.PL.yield <- pmc(ml.rooted.um, yield, "BM", "lambda", nboot = iter)

BM.PL.beta.log10.no812 <- pmc(ml.rooted.no812.um, beta.log10.no812, "BM", "lambda", nboot = iter)
BM.PL.alpha.no812 <- pmc(ml.rooted.no812.um, alpha.no812, "BM", "lambda", nboot = iter)
BM.PL.umax.no812 <- pmc(ml.rooted.no812.um, umax.no812, "BM", "lambda", nboot = iter)
BM.PL.lag.no812 <- pmc(ml.rooted.no812.um, lag.no812, "BM", "lambda", nboot = iter)
BM.PL.yield.no812 <- pmc(ml.rooted.no812.um, yield.no812, "BM", "lambda", nboot = iter)

# get p values
p_value.BM.PL.beta.log10 <- length(BM.PL.beta.log10$null[BM.PL.beta.log10$null > BM.PL.beta.log10$lr] ) / iter
p_value.BM.PL.alpha <- length(BM.PL.alpha$null[BM.PL.alpha$null > BM.PL.alpha$lr] ) / iter
p_value.BM.PL.umax <- length(BM.PL.umax$null[BM.PL.umax$null > BM.PL.umax$lr] ) / iter
p_value.BM.PL.lag <- length(BM.PL.lag$null[BM.PL.lag$null > BM.PL.lag$lr] ) / iter
p_value.BM.PL.yield <- length(BM.PL.yield$null[BM.PL.yield$null > BM.PL.yield$lr] ) / iter


p_value.BM.PL.beta.log10.no812 <- length(BM.PL.beta.log10.no812$null[BM.PL.beta.log10.no812$null > BM.PL.beta.log10.no812$lr] ) / iter
p_value.BM.PL.alpha.no812 <- length(BM.PL.alpha.no812$null[BM.PL.alpha.no812$null > BM.PL.alpha.no812$lr] ) / iter
p_value.BM.PL.umax.no812 <- length(BM.PL.umax.no812$null[BM.PL.umax.no812$null > BM.PL.umax.no812$lr] ) / iter
p_value.BM.PL.lag.no812 <- length(BM.PL.lag.no812$null[BM.PL.lag.no812$null > BM.PL.lag.no812$lr] ) / iter
p_value.BM.PL.yield.no812 <- length(BM.PL.yield.no812$null[BM.PL.yield.no812$null > BM.PL.yield.no812$lr] ) / iter


df.summary <- cbind(parameter=c("beta.log10","alpha"), 
                    test=c("BM_PL","BM_PL"), 
                    llr=c(BM.PL.beta.log10$lr, BM.PL.alpha$lr),
                    p_value=c(p_value.BM.PL.beta.log10, p_value.BM.PL.alpha),
                    sigsq_BM=c(BM.PL.beta.log10$A$opt$sigsq, BM.PL.alpha$A$opt$sigsq),
                    sigsq_PL=c(BM.PL.beta.log10$B$opt$sigsq, BM.PL.alpha$B$opt$sigsq),
                    lambda=c(BM.PL.beta.log10$B$opt$lambda, BM.PL.alpha$B$opt$lambda)
)

df.summary <- as.data.frame(df.summary)
df.summary$llr <- as.numeric(as.character(df.summary$llr))
df.summary$p_value <- as.numeric(as.character(df.summary$p_value))
write.csv(df.summary, file = "data/pmc/pmc_death_BM_PL_summary.csv")


df.summary.no812 <- cbind(parameter=c("beta.log10","alpha", "umax", "lag", "yield"), 
                    test=c("BM_PL","BM_PL","BM_PL","BM_PL","BM_PL"), 
                    llr=c(BM.PL.beta.log10.no812$lr, BM.PL.alpha.no812$lr, BM.PL.umax.no812$lr, BM.PL.lag.no812$lr, BM.PL.yield.no812$lr),
                    p_value=c(p_value.BM.PL.beta.log10.no812, p_value.BM.PL.alpha.no812, p_value.BM.PL.umax.no812, p_value.BM.PL.lag.no812, p_value.BM.PL.yield.no812),
                    sigsq_BM=c(BM.PL.beta.log10.no812$A$opt$sigsq, BM.PL.alpha.no812$A$opt$sigsq, BM.PL.umax.no812$A$opt$sigsq, BM.PL.lag.no812$A$opt$sigsq, BM.PL.yield.no812$A$opt$sigsq),
                    sigsq_PL=c(BM.PL.beta.log10.no812$B$opt$sigsq, BM.PL.alpha.no812$B$opt$sigsq, BM.PL.umax.no812$B$opt$sigsq, BM.PL.lag.no812$B$opt$sigsq, BM.PL.yield.no812$B$opt$sigsq),
                    lambda=c(BM.PL.beta.log10.no812$B$opt$lambda, BM.PL.alpha.no812$B$opt$lambda, BM.PL.umax.no812$B$opt$lambda, BM.PL.lag.no812$B$opt$lambda, BM.PL.yield.no812$B$opt$lambda)
)


df.summary.no812 <- as.data.frame(df.summary.no812)
df.summary.no812$llr <- as.numeric(as.character(df.summary.no812$llr))
df.summary.no812$p_value <- as.numeric(as.character(df.summary.no812$p_value))
write.csv(df.summary.no812, file = "data/pmc/pmc_death_BM_PL_no812_summary.csv")



# do the same analysis for the three growth curve parameters




















# w/out bacillus
test.beta.log10.lambda <- pmc(ml.rooted.um, beta.log10, "BM", "lambda", nboot = iter)
test.alpha.lambda <- pmc(ml.rooted.um, alpha, "BM", "lambda", nboot = iter)

test.mttf.lambda <- pmc(ml.rooted.um, mttf.log10, "BM", "lambda", nboot = 2)



iter <- 1000
BM.OU.beta.log10 <- pmc(ml.rooted.um, beta.log10, "BM", "OU", nboot = iter)
BM.OU.alpha <- pmc(ml.rooted.um, alpha, "BM", "OU", nboot = iter)
#BM.OU.mttf.log10 <- pmc(ml.rooted.um, mttf.log10, "BM", "OU", nboot = 2)

BM.GBM.beta.log10 <- pmc(ml.rooted.um, beta.log10, "BM", "trend", nboot = iter)
BM.GBM.alpha <- pmc(ml.rooted.um, alpha, "BM", "trend", nboot = iter)
#BM.GBM.mttf.log10 <- pmc(ml.rooted.um, mttf.log10, "BM", "trend", nboot = iter)

# run after this
GBM.OU.beta.log10 <- pmc(ml.rooted.um, beta.log10, "trend", "OU", nboot = iter)
GBM.OU.alpha <- pmc(ml.rooted.um, alpha, "trend", "OU", nboot = iter)
#GBM.OU.mttf.log10 <- pmc(ml.rooted.um, mttf.log10, "trend", "OU", nboot = iter)

p_value.BM.OU.beta.log10 <- length(BM.OU.beta.log10$null[BM.OU.beta.log10$null > BM.OU.beta.log10$lr] ) / iter
p_value.BM.GBM.beta.log10 <- length(BM.GBM.beta.log10$null[BM.GBM.beta.log10$null > BM.GBM.beta.log10$lr] ) / iter
p_value.GBM.OU.beta.log10 <- length(GBM.OU.beta.log10$null[GBM.OU.beta.log10$null > GBM.OU.beta.log10$lr] ) / iter

p_value.BM.OU.alpha <- length(BM.OU.alpha$null[BM.OU.alpha$null > BM.OU.alpha$lr] ) / iter
p_value.BM.GBM.alpha <- length(BM.GBM.alpha$null[BM.GBM.alpha$null > BM.GBM.alpha$lr] ) / iter
p_value.GBM.OU.alpha <- length(GBM.OU.alpha$null[GBM.OU.alpha$null > GBM.OU.alpha$lr] ) / iter

#p_value.BM.OU.mttf.log10 <- length(BM.OU.mttf.log10$null[BM.OU.mttf.log10$null > BM.OU.mttf.log10$lr] ) / iter
#p_value.BM.GBM.mttf.log10 <- length(BM.GBM.mttf.log10$null[BM.GBM.mttf.log10$null > BM.GBM.mttf.log10$lr] ) / iter
#p_value.GBM.OU.mttf.log10 <- length(GBM.OU.mttf.log10$null[GBM.OU.mttf.log10$null > GBM.OU.mttf.log10$lr] ) / iter



df.summary <- cbind(parameter=c("beta.log10","beta.log10","beta.log10","alpha","alpha","alpha"), 
                    test=c("BM_OU","BM_GBM", "GBM_OU", "BM_OU","BM_GBM", "GBM_OU"), 
                    llr=c(BM.OU.beta.log10$lr, BM.GBM.beta.log10$lr, GBM.OU.beta.log10$lr, BM.OU.alpha$lr, BM.GBM.alpha$lr, GBM.OU.alpha$lr),
                    p_value=c(p_value.BM.OU.beta.log10, p_value.BM.GBM.beta.log10, p_value.GBM.OU.beta.log10, p_value.BM.OU.alpha, p_value.BM.GBM.alpha, p_value.GBM.OU.alpha),
                    sigsq_BM=c(BM.OU.beta.log10$A$opt$sigsq, BM.GBM.beta.log10$A$opt$sigsq, NA, BM.OU.alpha$A$opt$sigsq, BM.GBM.alpha$A$opt$sigsq, NA),
                    sigsq_OU=c(BM.OU.beta.log10$B$opt$sigsq, NA, GBM.OU.beta.log10$B$opt$sigsq, BM.OU.alpha$B$opt$sigsq, NA,  GBM.OU.alpha$B$opt$sigsq),
                    alpha_OU=c(BM.OU.beta.log10$B$opt$alpha, NA, GBM.OU.beta.log10$B$opt$alpha, BM.OU.alpha$B$opt$alpha, NA,  GBM.OU.alpha$B$opt$alpha),
                    sigsq_GBM=c(NA, BM.GBM.beta.log10$B$opt$sigsq, GBM.OU.beta.log10$A$opt$sigsq, NA, BM.GBM.alpha$B$opt$sigsq, GBM.OU.alpha$A$opt$sigsq),
                    slope_GBM=c(NA, BM.GBM.beta.log10$B$opt$slope, GBM.OU.beta.log10$A$opt$slope, NA, BM.GBM.alpha$B$opt$slope, GBM.OU.alpha$A$opt$slope)
                    )


#df.summary <- cbind(parameter=c("beta.log10","beta.log10","beta.log10","alpha","alpha","alpha", "mttf.log10", "mttf.log10", "mttf.log10"), 
#                    test=c("BM_OU","BM_GBM", "GBM_OU", "BM_OU","BM_GBM", "GBM_OU", "BM_OU","BM_GBM", "GBM_OU"), 
#                    llr=c(BM.OU.beta.log10$lr, BM.GBM.beta.log10$lr, GBM.OU.beta.log10$lr, BM.OU.alpha$lr, BM.GBM.alpha$lr, GBM.OU.alpha$lr, BM.OU.mttf.log10$lr, BM.GBM.mttf.log10$lr, GBM.OU.mttf.log10$lr),
#                    p_value=c(p_value.BM.OU.beta.log10, p_value.BM.GBM.beta.log10, p_value.GBM.OU.beta.log10, p_value.BM.OU.alpha, p_value.BM.GBM.alpha, p_value.GBM.OU.alpha, p_value.BM.OU.mttf.log10, p_value.BM.GBM.mttf.log10, p_value.GBM.OU.mttf.log10),
#                    sigsq_BM=c(BM.OU.beta.log10$A$opt$sigsq, BM.GBM.beta.log10$A$opt$sigsq, NA, BM.OU.alpha$A$opt$sigsq, BM.GBM.alpha$A$opt$sigsq, NA, BM.OU.mttf.log10$A$opt$sigsq, BM.GBM.mttf.log10$A$opt$sigsq, NA),
#                    sigsq_OU=c(BM.OU.beta.log10$B$opt$sigsq, NA, GBM.OU.beta.log10$B$opt$sigsq, BM.OU.alpha$B$opt$sigsq, NA,  GBM.OU.alpha$B$opt$sigsq, BM.OU.mttf.log10$B$opt$sigsq, NA,  GBM.OU.mttf.log10$B$opt$sigsq),
#                    alpha_OU=c(BM.OU.beta.log10$B$opt$alpha, NA, GBM.OU.beta.log10$B$opt$alpha, BM.OU.alpha$B$opt$alpha, NA,  GBM.OU.alpha$B$opt$alpha, BM.OU.mttf.log10$B$opt$alpha, NA,  GBM.OU.mttf.log10$B$opt$alpha),
#                    sigsq_GBM=c(NA, BM.GBM.beta.log10$B$opt$sigsq, GBM.OU.beta.log10$A$opt$sigsq, NA, BM.GBM.alpha$B$opt$sigsq, GBM.OU.alpha$A$opt$sigsq, NA, BM.GBM.mttf.log10$B$opt$sigsq, GBM.OU.mttf.log10$A$opt$sigsq),
#                    slope_GBM=c(NA, BM.GBM.beta.log10$B$opt$slope, GBM.OU.beta.log10$A$opt$slope, NA, BM.GBM.alpha$B$opt$slope, GBM.OU.alpha$A$opt$slope, NA, BM.GBM.mttf.log10$B$opt$slope, GBM.OU.mttf.log10$A$opt$slope)
#)

df.summary <- as.data.frame(df.summary)
df.summary$llr <- as.numeric(as.character(df.summary$llr))
df.summary$p_value <- as.numeric(as.character(df.summary$p_value))

write.csv(df.summary, file = "data/pmc/pmc_summary.csv")



# write code for table here
#######




