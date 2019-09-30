rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('lme4')
library('nlme')
library('optimx')

df <- read.csv("data/demography/weibull_results_clean.csv", 
               header = TRUE, stringsAsFactors = FALSE)

df<-df[!(df$strain=="KBS0727"),]
#df<-df[!(df$strain=="KBS0812"),]


df$N_0_log10 <- log10(df$N_0)

me0.inter <- lmer(alpha ~ N_0_log10, data = df)
ols <- lm(alpha ~ N_0_log10, data = df)
mixed.inter <- lmer(alpha ~ N_0_log10 + (1|strain), data = df)

AIC(ols, mixed.inter)



summary(me.inter)

summary(lmer(alpha ~ (1|strain), data = df))

summary(lm(alpha ~ N_0_log10 , data = df))


lmList(alpha ~ N_0_log10 | (1|strain), data=df)



basic.lm <- lm(alpha ~ N_0_log10, data = df)
plot(basic.lm, which = 1)
plot(basic.lm, which = 2)


mixed.lmer <- lmer(alpha ~ N_0_log10+ (N_0_log10|strain), data = df)
                                                                                         
mixed.lmer <- lmer(alpha ~ N_0_log10+ (N_0_log10|strain), data = df, control=lmerControl(optimizer="optimx",
                                                                                         optCtrl=list(method='nlminb')))
#lmer(alpha ~ N_0_log10+ (N_0_log10|strain))
print(VarCorr(mixed.lmer),comp="Variance")

coef(mixed.lmer)$strain


plot(df$N_0_log10 ~ df$alpha)

df$N_0_beta <- df$N_0/df$beta
df$N_0_beta_log10 <- log10(df$N_0_beta)






plot(log10(df$N_0_beta_log10), df$alpha, xlab="N0 / lambda", ylab="k")
abline(lm(df$alpha ~ log10(df$N_0_beta_log10)))

df <- log10(df$N_0_beta_log10)
basic.lm <- lm(alpha ~ log10(N_0_beta_log10), data = df)
plot(basic.lm, which = 1)
plot(basic.lm, which = 2)





# fit non-linear model
mod <- nls(alpha ~ exp(a + b * df$N_0_beta_log10), data = df, start = list(a = 0, b = 0))

# add fitted curve
lines(df$N_0_beta_log10, predict(mod, list(x = df$N_0_beta_log10)))

