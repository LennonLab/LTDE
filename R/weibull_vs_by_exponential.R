rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('bbmle')
library('devtools')
library('pracma')
library('plotrix')
#install_github("rmcelreath/rethinking")
#library('rethinking')

## Load Data
obs <- read.csv("data/demography/longtermdormancy_20190528_nocomments.csv", 
                header = TRUE, stringsAsFactors = FALSE)
## Adding 1 to deal with log(0) observations
obs$Abund <- (as.numeric(obs$Colonies) +1)* (1000 / as.numeric(obs$Inoculum )) * ( 10 ^  as.numeric(obs$Dilution) )
strains <- sort(unique(obs$Strain))
strains <- strains[table(obs$Strain)>10]
#strains <- c('KBS0711')
obs <- obs[obs$Strain%in%strains,]
summ <- matrix(NA,length(strains)*max(obs$Rep),20)
counter <- 1

strainObs=obs[obs$Strain==strains[1],]
reps=unique(strainObs$Rep)
repObs=strainObs[strainObs$Rep==reps[1],]

time<-(as.numeric(strptime(repObs$Firstread_date,format="%d-%b-%y",tz="EST"))-
         as.numeric(strptime(repObs$Dormstart_date,format="%d-%b-%y",tz="EST")))/(3600*24)

index <- order(time)
time <- time[index]
abund <- repObs$Abund[index]
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
      print(length(time))
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
      repObs["prop"] <- abund / N_0
      repObs["log_prop"] <- log(abund / N_0)
      # Initial parameters
      #beta = Initial death (larger = slower) 
      #alpha = Bend (upper = 1 = first-order decay)
      #Z = Error
      grids<-list(beta=c(1,10,50,100,200),alpha=c(0.05,0.1,0.5,1),z=c(0.1,1,10))
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
        fit <- mle2(minuslogl=log_prop ~ dnorm(mean =  -1 * ((time / beta)^ alpha), sd = z), 
                    start = new.start, data = repObs,
                    control=list(parscale=pscale, maxit=1000), 
                    method="Nelder-Mead", hessian = T)
        
        #fit <- mle2(minuslogl=prop ~ dnorm(mean =  -1 * ((time / beta)^ alpha), sd = z), 
        #            start = new.start, data = repObs, lower = c(beta = 2, alpha =0.001), upper = c(beta=1000, alpha =1),
        #            control=list(parscale=pscale, maxit=1000), 
        #            method="L-BFGS-B", hessian = T)
        
        res.mat[k,]<-c(coef(fit),AIC(fit))	
        res.mod[[k]]<-fit
      }
      
      # same thing for null model
      grids.exp<-list(delta=c(0.0001,0.1,0.99999), lambda_1=c(0.1, 1,10), lambda_2 = c(0.1,1,10), z=c(0.1,1,10))
      start.exp<-list(delta=NA,lambda_1=NA,lambda_2=NA,z=NA)
      grid.starts.exp<-as.matrix(expand.grid(grids.exp))
      ncombos.exp<-dim(grid.starts.exp)[[1]]
      # cycle through each combo
      res.mat.exp<-matrix(NA,nrow=ncombos.exp,ncol=I(length(start.exp)+1))
      res.mod.exp<-list()
      
      vector <- c()
      for(k in 1:dim(grid.starts.exp)[[1]]){
        #some how need to match grid parameters to start lists.
        mod.start.exp<-as.list(grid.starts.exp[k,])	
        new.start.exp<-start.exp
        new.start.exp[names(start.exp) %in% names(mod.start.exp)]<-mod.start.exp
        #print(new.start.exp)
        pscale.exp<-as.numeric(new.start.exp)
        names(pscale.exp)<-names(new.start.exp)
        #print(pscale.exp)

        #fit.exp <- mle2(minuslogl=log_prop ~ dnorm(mean =  log(delta* exp(-1* time * lambda_1 ) + (1-delta) * exp(-1* time * lambda_2 ) ), sd = z), 
        #            start = new.start.exp, data = repObs, lower = c(delta = 0, lambda_1 = 0.001, lambda_2 =0.001), upper = c(delta=1, lambda_1=100, lambda_2=100),
        #            control=list(parscale=pscale.exp, maxit=1000), 
        #            method="L-BFGS-B", hessian = F)
        
        
        
        f1 = function(par, x){
          log(par[0]* exp(-1* time * par[1] ) + (1-par[0]) * exp(-1* time * par[2] ) )
        }
        result <- optim(log(abund / N_0), f1, x=time, 
                        method = "L-BFGS-B", 
                        lower = c(0, 0.001, 0.001), upper = c(1, 100, 100),
                        control = list(trace = 5,fnscale=-1)) 
        
        #skip_to_next <- FALSE
        
        print(result)
        
        # Note that print(b) fails since b doesn't exist
        
        #tryCatch(print(b), error = function(e) { skip_to_next <<- TRUE})
        
        #if(skip_to_next) { next }   
        
        #print(fit.exp)
        #print(AIC(fit.exp))
        #res.mat.exp[k,]<-c(coef(fit.exp),AIC(fit.exp))	
        #res.mod.exp[[k]]<-fit.exp
       # if(inherits(fit.exp, "try-error"))
        #{
        #  #error handling code, maybe just skip this iteration using
        #  next
        #}
        
        print( AIC(fit.exp))
        
        vector <- c(vector, AIC(fit.exp))
      }
      
      print(vector)
      
      colnames(res.mat)<-c(names(coef(fit)),"AIC")
      best.fit<-res.mod[[which(res.mat[,'AIC']==min(res.mat[,'AIC']))[1]]]
      
      print(best.fit)
      
      #colnames(res.mat.exp)<-c(names(coef(fit.exp)),"AIC")
      best.fit.exp<-res.mod.exp[[which(res.mat.exp[,'AIC']==min(res.mat.exp[,'AIC']))[1]]]
     
      print(AIC(best.fit))
      print(AIC(best.fit.exp))
      #LLR <- -2* (as.numeric(logLik(best.fit.exp) - logLik(best.fit)))
      #p.value <- pchisq(LLR,df=1,lower.tail=FALSE)
      

      #beta <- coef(best.fit)[1]
      #alpha <- coef(best.fit)[2]
      
      #summ[counter,1]=strains[i]
      #summ[counter,2]=reps[j]
      # beta
      #summ[counter,3]=beta
      # alpha
     # summ[counter,4]=alpha
      # z

    
    }
  }
}

#dev.off()
#summ=summ[!is.na(summ[,1]),]
#colnames(summ)=c('strain','rep','beta','alpha','std_dev','AIC', 'N.obs', 'beta.sd', 'alpha.sd', 'z.sd', 'mttf', 'mttf.sd', 'log10.mttf.sd',  "N_0", "N_final", "LR", "p.value", "Last_date", "T_ext", "log10.T_ext.sd")
#write.csv(summ,"data/demography/weibull_results.csv")


