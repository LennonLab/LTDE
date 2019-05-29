rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('bbmle')
library('devtools')
#install_github("rmcelreath/rethinking")
#library('rethinking')

#library('rethinking')
## Load Data
obs <- read.csv("data/demography/longtermdormancy_20170620_nocomments.csv", 
                header = TRUE, stringsAsFactors = FALSE)
## Adding 1 to deal with log(0) observations
#obs$Abund <- as.numeric(obs$Colonies) * 10 ^ as.numeric(obs$Dilution) + 1
obs$Abund <- (as.numeric(obs$Colonies) +1)* (1000 / as.numeric(obs$Inoculum )) * ( 10 ^  as.numeric(obs$Dilution) )
strains <- sort(unique(obs$Strain))
strains <- strains[table(obs$Strain)>10]
#strains <- c('KBS0721')
#print(strains[-c('KBS0714', 'KBS0715')])


obs <- obs[obs$Strain%in%strains,]
summ <- matrix(NA,length(strains)*max(obs$Rep),15)
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
      time=(as.numeric(strptime(repObs$Firstread_date,format="%d-%b-%y",tz="EST"))-
              as.numeric(strptime(start,format="%d-%b-%y",tz="EST")))/(3600*24)
      repObs["time"] <- time + 1
      repObs["logabund"] <- log10(repObs$Abund)
      if (repObs["logabund"][[1]][2] - repObs["logabund"][[1]][1] > 1){
        repObs <- repObs[-c(1), ]
      }
      repObs["prop"] <- log(repObs$Abund / repObs$Abund[1])
      N_0 <- repObs$Abund[1]
      # Initial parameters
      #beta = Initial death (larger = slower) 
      #alpha = 1 # Bend (upper = 1 = first-order decay)
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
      colnames(res.mat)<-c(names(coef(fit)),"AIC")
      best.fit<-res.mod[[which(res.mat[,'AIC']==min(res.mat[,'AIC']))[1]]]
      
      
      summ[counter,1]=strains[i]
      summ[counter,2]=reps[j]
      # beta
      beta <- coef(best.fit)[1]
      summ[counter,3]=beta
      # alpha
      alpha <- coef(best.fit)[2]
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
      summ[counter,13] <- sqrt( (beta**2) * (gamma(1 + (2/alpha)) - (gamma(1 + (1/alpha)) **2)) )
      dV_dBeta <- 2 * beta * (gamma(1 + (2/alpha)) - (gamma(1 + (1/alpha)) **2))
      dV_dAlpha <- (beta**2) * (gamma(1 + (2/alpha)) * digamma(1 + (2/alpha)) - (2*gamma(1 + (1/alpha)) * gamma(1 + (1/alpha)) * digamma(1 + (1/alpha)))  )
      dV_vector <- c(dV_dBeta, dV_dAlpha)
      summ[counter,14] <- sqrt(t(dV_vector) %*% best.fit@vcov[1:2,1:2] %*% dV_vector)
      summ[counter,15] <- N_0
      
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
colnames(summ)=c('strain','rep','beta','alpha','std_dev','AIC', 'N.obs', 'beta.sd', 'alpha.sd', 'z.sd', 'mttf', 'mttf.sd', "sd.ttf", "sd.ttf.sd", "N_0")
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

write.csv(df, file = "data/demography/weibull_results_clean.csv")


# get mean time to failure and CIs
df.species.mean <- aggregate(df[, c('beta', 'alpha', 'mttf', 'N_0')], list(df$strain), mean)
colnames(df.species.mean)[1] <- "Species"

# function to calculate pooled standard error
get.pooled.se <- function(strains){
  df.strain.new <- data.frame()
  for (strain in strains)
  {
    df.strain <- df[ which(df$strain==strain), ]
    # remove rows with NAs
    df.strain <- df.strain[complete.cases(df.strain), ]
    N.reps <- nrow(df.strain)
    N_0.se <- sd(df.strain$N_0) / sqrt(N.reps)
    
    pooled.mttf.var <- sum((df.strain$N.obs-1) * (df.strain$mttf.sd ** 2)) / sum(df.strain$N.obs-1)
    pooled.beta.var <- sum((df.strain$N.obs-1) * (df.strain$beta.sd ** 2)) / sum(df.strain$N.obs-1)
    pooled.alpha.var <- sum((df.strain$N.obs-1) * (df.strain$alpha.sd ** 2)) / sum(df.strain$N.obs-1)
    
    pooled.mttf.se <- sqrt(pooled.mttf.var) / sqrt(N.reps)
    pooled.beta.se <- sqrt(pooled.beta.var) / sqrt(N.reps)
    pooled.alpha.se <- sqrt(pooled.alpha.var) / sqrt(N.reps)
    df.strain.new.row <- data.frame(strain, pooled.mttf.se, pooled.beta.se, pooled.alpha.se, N_0.se)
    df.strain.new <- rbind(df.strain.new, df.strain.new.row)
  }
  return(df.strain.new)
}

pooled.se <- get.pooled.se(unique(df$strain))
colnames(pooled.se)[1] <- "Species"

df.species <- merge(df.species.mean, pooled.se,by="Species")


write.csv(df.species, file = "data/demography/weibull_results_clean_species.csv")

