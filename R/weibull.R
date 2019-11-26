rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('bbmle')
library('devtools')
library('plotrix')
library('pracma')
#install_github("rmcelreath/rethinking")
#library('rethinking')

## Load Data
obs <- read.csv("data/demography/longtermdormancy_20190528_nocomments.csv", 
                header = TRUE, stringsAsFactors = FALSE)
## Adding 1 to deal with log(0) observations
obs$Abund <- (as.numeric(obs$Colonies) +1)* (1000 / as.numeric(obs$Inoculum )) * ( 10 ^  as.numeric(obs$Dilution) )
strains <- sort(unique(obs$Strain))
#strains <- strains[table(obs$Strain)>10]
#strains <- c('KBS0713')
obs <- obs[obs$Strain%in%strains,]
summ <- matrix(NA,length(strains)*max(obs$Rep),18)
pdf('figs/weibull_fits.pdf') # Uncomment to create pdf that will plot data and fits
counter <- 1
for(i in 1:length(strains)){
  strainObs=obs[obs$Strain==strains[i],]
  reps=unique(strainObs$Rep)
  print(strains[i])
  for(j in 1:length(reps)){
    repObs=strainObs[strainObs$Rep==reps[j],]
    # minimum of 10 data points
    if(nrow(repObs)>10){
      start=repObs[1,1]
      time<-(as.numeric(strptime(repObs$Firstread_date,format="%d-%b-%y",tz="EST"))-
              as.numeric(strptime(repObs$Dormstart_date,format="%d-%b-%y",tz="EST")))/(3600*24)
      #time=(as.numeric(strptime(repObs$Firstread_date,format="%d-%b-%y",tz="EST"))-
      #        as.numeric(strptime(start,format="%d-%b-%y",tz="EST")))/(3600*24)
      # sort the data
      index <- order(time)
      time <- time[index]
      abund <- repObs$Abund[index]
      # add one to each day to avoid multiplication by zero error
      repObs["time"] <- time + 1
      repObs["logabund"] <- log10(abund)
      if (repObs["logabund"][[1]][2] - repObs["logabund"][[1]][1] > 1){
        repObs <- repObs[-c(1), ]
      }

      N_0 <- abund[1]
      N_final <- abund[length(abund)]
      repObs["prop"] <- log(abund / N_0)
      # Initial parameters
      #beta = Initial death (larger = slower) 
      #alpha = Bend (upper = 1 = first-order decay)
      #Z = Error
      grids<-list(beta=c(1,10,50,100,200),alpha=c(0.05,0.1,0.5,1,1.1,1.5),z=c(0.1,1,10))
      start<-list(beta=NA,alpha=NA,z=NA)
      grid.starts<-as.matrix(expand.grid(grids))
      ncombos<-dim(grid.starts)[[1]]
      # cycle through each combo
      res.mat<-matrix(NA,nrow=ncombos,ncol=I(length(start)+1))
      res.mod<-list()
      for(k in 1:dim(grid.starts)[[1]]){
        #some how need to match grid parameters to start lists.
        mod.start<-as.list(grid.starts[k,])	
        new.start<-start
        new.start[names(start) %in% names(mod.start)]<-mod.start
        pscale<-as.numeric(new.start)
        names(pscale)<-names(new.start)
        fit <- mle2(minuslogl=prop ~ dnorm(mean =  -1 * ((time / beta)^ alpha), sd = z), 
                    start = new.start, data = repObs,
                    control=list(parscale=pscale, maxit=1000), 
                    method="Nelder-Mead", hessian = T)
        res.mat[k,]<-c(coef(fit),AIC(fit))	
        res.mod[[k]]<-fit
      }
      
      # same thing for null model
      grids.exp<-list(beta=c(1,10,50,100,200),z=c(0.1,1,10))
      start.exp<-list(beta=NA,z=NA)
      grid.starts.exp<-as.matrix(expand.grid(grids.exp))
      ncombos.exp<-dim(grid.starts.exp)[[1]]
      # cycle through each combo
      res.mat.exp<-matrix(NA,nrow=ncombos.exp,ncol=I(length(start.exp)+1))
      res.mod.exp<-list()
      for(k in 1:dim(grid.starts.exp)[[1]]){
        #some how need to match grid parameters to start lists.
        mod.start.exp<-as.list(grid.starts.exp[k,])	
        new.start.exp<-start.exp
        new.start.exp[names(start.exp) %in% names(mod.start.exp)]<-mod.start.exp
        pscale.exp<-as.numeric(new.start.exp)
        names(pscale.exp)<-names(new.start.exp)
        fit.exp <- mle2(minuslogl=prop ~ dnorm(mean =  -1 * ((time / beta)), sd = z), 
                    start = new.start.exp, data = repObs,
                    control=list(parscale=pscale.exp, maxit=1000), 
                    method="Nelder-Mead", hessian = T)
        res.mat.exp[k,]<-c(coef(fit.exp),AIC(fit.exp))	
        res.mod.exp[[k]]<-fit.exp
      }
      
      colnames(res.mat)<-c(names(coef(fit)),"AIC")
      best.fit<-res.mod[[which(res.mat[,'AIC']==min(res.mat[,'AIC']))[1]]]
      
      colnames(res.mat.exp)<-c(names(coef(fit.exp)),"AIC")
      best.fit.exp<-res.mod.exp[[which(res.mat.exp[,'AIC']==min(res.mat.exp[,'AIC']))[1]]]
      LLR <- -2* (as.numeric(logLik(best.fit.exp) - logLik(best.fit)))
      p.value <- pchisq(LLR,df=1,lower.tail=FALSE)
      
      beta <- coef(best.fit)[1]
      alpha <- coef(best.fit)[2]
  
      summ[counter,1]=strains[i]
      summ[counter,2]=reps[j]
      # beta
      summ[counter,3]=beta
      # alpha
      summ[counter,4]=alpha
      # z
      summ[counter,5]=coef(best.fit)[3]
      summ[counter,6]=AIC(best.fit)
      summ[counter,7]=length(repObs$time)
      #CIs <- confint(profile(best.fit))
      best.fit.summary <- summary(best.fit)
      # S.E. beta
      summ[counter,8]=best.fit.summary@coef[1,2]
      # S.E. alpha
      summ[counter,9]=best.fit.summary@coef[2,2]
      # S.E. z
      summ[counter,10]=best.fit.summary@coef[3,2]
      # MTTF 
      summ[counter,11] <- coef(best.fit)[1] * gamma(1 + (1/coef(best.fit)[2]))
      
      
      dT_dBeta <- gamma(1 + (1/alpha))
      dT_dAlpha <- -1* (beta/ (alpha**2)) * gamma(1 + (1/alpha)) * digamma(1 + (1/alpha))
      dT_vector <- c(dT_dBeta, dT_dAlpha)
      summ[counter,12] <- sqrt(t(dT_vector) %*% best.fit@vcov[1:2,1:2] %*% dT_vector)

      dlog10T_dBeta <-log10(exp(1)) * (beta**-1)
      dlog10T_dAlpha <- -1*(alpha**-2)*log10(exp(1)) * digamma(1 + (1/alpha))
      dlog10T_vector <- c(dlog10T_dBeta, dlog10T_dAlpha)
      summ[counter,13] <- sqrt(t(dlog10T_vector) %*% best.fit@vcov[1:2,1:2] %*% dlog10T_vector)
      
      summ[counter,14] <- N_0
      summ[counter,15] <- N_final
      summ[counter,16] <- LLR
      summ[counter,17] <- p.value
      summ[counter,18] <- max(repObs["time"]-1)

      ### *** Comment/Uncomment following code to make pdf figs*** ###
      title=paste(strains[i],"  rep ",reps[j])
      plot(repObs$time,repObs$prop,main=title,ylim=c(min(repObs$prop),0), 
           xlab = 'Time (days)', ylab = 'Proportion surviving, log' )
      predTime=seq(0,max(repObs$time))
      lines(repObs$time, (-1 * ((repObs$time /beta )^ alpha )), 
            lwd=4, lty=2, col = "red")
      counter=counter+1
    }
  }
}

dev.off()
summ=summ[!is.na(summ[,1]),]
colnames(summ)=c('strain','rep','beta','alpha','std_dev','AIC', 'N.obs', 'beta.sd', 'alpha.sd', 'z.sd', 'mttf', 'mttf.sd', 'log10.mttf.sd',  "N_0", "N_final", "LR", "p.value", "Last_date")
write.csv(summ,"data/demography/weibull_results.csv")



# clean the results file
df <- read.table("data/demography/weibull_results.csv", 
                 header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
df <- df[!(df$strain=="KBS0711W"),]
# remove KBS0711 replicates 10, 11, and 12
# these samples were only sampled starting at day 100
# We don't have time to re-do them, so just remove them
df <- df[!(df$strain == "KBS0711" & df$rep == 10 ),] 
df <- df[!(df$strain == "KBS0711" & df$rep == 11 ),] 
df <- df[!(df$strain == "KBS0711" & df$rep == 12 ),] 
df <- df[!(df$strain=="KBS0727"),]
df$N_0_beta <- df$N_0/df$beta

# multiple testing correction 
df$p.value.BH <- p.adjust(df$p.value, method = "BH", n = length(df$p.value))
write.csv(df, file = "data/demography/weibull_results_clean.csv")


# get mean time to failure and CIs
df$beta.log10 <- log10(df$beta)
df$alpha.log10 <- log10(df$alpha)
df$mttf.log10 <- log10(df$mttf)
df$N_0.log10 <- log10(df$N_0)
df$N_0_beta.log10 <- log10(df$N_0_beta)

df$N_final.log10 <- log10(df$N_final)
df$delta_N.log10 <- log10(df$N_0 - df$N_final) 
df.species.mean <- aggregate(df[, c('beta', 'alpha', 'mttf', 'N_0', 'N_final', 'beta.log10','mttf.log10', 'delta_N.log10', 'Last_date', 'LR', 'p.value.BH', 'N_0.log10', 'N_final.log10', 'N_0_beta.log10', 'alpha.log10')], list(df$strain), mean)
colnames(df.species.mean)[1] <- "Species"

df.species.log10.se <- aggregate(df[, c('beta.log10', 'mttf.log10')], list(df$strain), std.error)
colnames(df.species.log10.se) <-  c("Species", "beta.log10.se", "mttf.log10.se")
df.species <- merge(df.species.mean, df.species.log10.se,by="Species")

# function to calculate pooled standard error
get.pooled.se <- function(strains){
  df.strain.new <- data.frame()
  for (strain in strains)
  {
    df.strain <- df[ which(df$strain==strain), ]
    N.reps <- nrow(df.strain)
    log10.N_0.se <- sd(log10(df.strain$N_0)) / sqrt(N.reps)
    log10.N_0 <- mean(log10(df.strain$N_0))
    # remove rows with NAs
    if(strain == "KBS0812"){
      pooled.log10.mttf.se <- sd(df.strain$mttf.log10) / sqrt(N.reps)
      pooled.alpha.se <- sd(df.strain$alpha) / sqrt(N.reps)
    }
    else {
      df.strain <- df.strain[complete.cases(df.strain), ]
      
      pooled.log10.mttf.var <- sum((df.strain$N.obs-1) * (df.strain$log10.mttf.sd ** 2)) / sum(df.strain$N.obs-1)
      pooled.alpha.var <- sum((df.strain$N.obs-1) * (df.strain$alpha.sd ** 2)) / sum(df.strain$N.obs-1)
      
      pooled.log10.mttf.se <- sqrt(pooled.log10.mttf.var) / sqrt(N.reps)
      pooled.alpha.se <- sqrt(pooled.alpha.var) / sqrt(N.reps)
    }
    
    df.strain.new.row <- data.frame(strain, pooled.log10.mttf.se, pooled.alpha.se, log10.N_0, log10.N_0.se)
    df.strain.new <- rbind(df.strain.new, df.strain.new.row)
  }
  return(df.strain.new)
}

pooled.se <- get.pooled.se(unique(df$strain))
colnames(pooled.se)[1] <- "Species"
df.species.pool <- merge(df.species, pooled.se,by="Species")
write.csv(df.species.pool, file = "data/demography/weibull_results_clean_species.csv")

