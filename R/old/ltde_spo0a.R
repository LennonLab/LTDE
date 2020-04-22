rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('bbmle')

obs <- read.csv("data/demography/spo0a_counts_pbs.csv", 
                header = TRUE, stringsAsFactors = FALSE)
obs$Abund <- (as.numeric(obs$nt) +1)* (1000 / as.numeric(100 )) * ( 10 ^  as.numeric(obs$nt_dil) )
obs$Days <- as.numeric(obs$time / 24)

reps <- c("A","B","C","D")

obs.A <- obs[(obs$rep == "A"), ] 
obs.B <- obs[(obs$rep == "B"), ] 
obs.C <- obs[(obs$rep == "C"), ] 
obs.D <- obs[(obs$rep == "D"), ] 

obs.A$log.prop <- as.numeric(log(obs.A$Abund/obs.A$Abund[1]))




fit <- mle2(minuslogl=log.prop.A ~ dnorm(mean =  -1 * ((Days / beta)^ alpha), sd = z), 
            start = c(beta=0.01, alpha=1, z=0.1), data = obs.A,
            control=list(parscale=pscale, maxit=1000), 
            method="Nelder-Mead", hessian = T)


fittt <- mle2(minuslogl=log.prop ~ dnorm(mean =  -1 * ((Days / beta)^ alpha), sd = z), 
              data=obs.A, start=c(beta=10, alpha=1, z=0.1))




repObs["prop"] <- log(abund / N_0)

plot(obs.A$Days, log10(obs.A$Abund), col="red")
points(obs.B$Days, log10(obs.B$Abund), col="green")
points(obs.C$Days, log10(obs.C$Abund), col="blue")
points(obs.D$Days, log10(obs.D$Abund), col="orange")


plot(obs.B$Days, log10(obs.B$Abund))


fit <- mle2(minuslogl=prop ~ dnorm(mean =  -1 * ((time / beta)^ alpha), sd = z), 
            start = new.start, data = repObs,
            control=list(parscale=pscale, maxit=1000), 
            method="Nelder-Mead", hessian = T)

