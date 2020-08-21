rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('phylolm')
library('MuMIn')
library('viridis')
library('lme4')


df <- read.csv("data/demography/weibull_results_clean.csv", 
               header = TRUE, stringsAsFactors = FALSE)

#df.sub <- df[df$strain %in% c('ATCC13985', 'ATCC43928', 'KBS0702', 'KBS0703', 'KBS0705', 'KBS0706', 'KBS0711', 'KBS0712', 'KBS0713', 'KBS0721', 'KBS0725', 'KBS0801'),]


df.species <- read.table("data/demography/weibull_results_clean_species.csv", 
                         header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
rownames(df.species) <- df.species$Species
df.species <- df.species[!(df.species$Species=="KBS0727"),]
df.species <- df.species[order(df.species$mttf.log10),]
rownames(df.species) <- df.species$Species

#### examine relationship between number of dead cells and shape paraemter
df$N_0_beta <- df$N_0/df$beta

mixed.b0 <- lmer(log10(alpha) ~ log10(N_0_beta) + (1|strain), data = df)
#AIC(ols, mixed.inter)
mixed.b0_and_b1 <- lmer(log10(alpha) ~ log10(N_0_beta) + (log10(N_0_beta)|strain), data = df)
qqnorm(resid(mixed.b0_and_b1))
qqline(resid(mixed.b0_and_b1))

# test whether random slopes explain enough variation
anova(mixed.b0, mixed.b0_and_b1)
# test whether the fixed effect alpha is significant (slope test)
summary(mixed.b0_and_b1  )

print(VarCorr(mixed.b0_and_b1),comp="Variance")

anova(mixed.b0, mixed.b0_and_b1)$"Pr(>Chisq)"[2]

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


phylolm.fit <- phylolm(log10(alpha) ~ N_0_beta.log10, data=df.species, phy=ml.rooted.um, model = 'lambda', boot=10)



# print coefficients to file

r2 <- r.squaredGLMM(mixed.b0_and_b1)

merBoot <- bootMer(mixed.b0_and_b1, predict, nsim = 10000, re.form = NA)
pred <- predict(mixed.b0_and_b1, re.form = NA)

std.err <- apply(merBoot$t, 2, sd)
CI.lower <- pred - std.err*1.96
CI.upper <- pred + std.err*1.96

# output CI as txt
CI.df <- data.frame(merBoot$data$`log10(N_0_beta)`, CI.lower, CI.upper)
colnames(CI.df) <- c('x', 'CI_lower', 'CI_upper')
CI.df <- CI.df[order(CI.df$x),]
# output coefficients as txt
coef.df <- coef(mixed.b0_and_b1)$strain


labels <- c('phylom_intercept', 'phylom_slope', 'lmm_intercept', 'lmm_slope', 'r2_m', 'p_value')
model_items <- c(phylolm.fit$coefficients[1], phylolm.fit$coefficients[2], mixed.b0_and_b1@beta[1], mixed.b0_and_b1@beta[2], r2[1], anova(mixed.b0, mixed.b0_and_b1)$"Pr(>Chisq)"[2])


model.df <- data.frame(labels, model_items)

write.csv(model.df, "data/demography/model_features.csv", row.names = FALSE)
write.csv(CI.df, "data/demography/model_CIs.csv", row.names = FALSE)
write.csv(coef.df, "data/demography/slope_coefficients.csv", row.names = FALSE)




