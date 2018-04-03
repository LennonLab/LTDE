rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('bbmle')
## Load Data
obs <- read.csv("data/demography/longtermdormancy_20170620_nocomments.csv", 
                header = TRUE, stringsAsFactors = FALSE)
## Adding 1 to deal with log(0) observations
obs$Abund <- as.numeric(obs$Colonies) * 10 ^ as.numeric(obs$Dilution) + 1
strains <- sort(unique(obs$Strain))
strains <- strains[table(obs$Strain)>10]
obs <- obs[obs$Strain%in%strains,]
summ <- matrix(NA,length(strains)*max(obs$Rep),7)
pdf('figs/weibull_fits.pdf') # Uncomment to create pdf that will plot data and fits
counter <- 1

for(i in 1:length(strains)){
  strainObs=obs[obs$Strain==strains[i],]
  
  reps=unique(strainObs$Rep)
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
      # Initial parameters
      #A = 200 # Initial death (larger = slower) 
      #B = 1 # Bend (upper = 1 = first-order decay)
      #C = round(max(repObs$logabund),1) # intercept
      #Z = 6 # Error
      grids<-list(a=c(1,10,50,100,200),b=c(0.05,0.1,0.5,1,1.1,1.5),z=c(0.1,1,10))
      #start<-list(a=NA,b=NA,c=round(max(repObs$logabund),1),z=NA)
      start<-list(a=NA,b=NA,z=NA)
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
        fit <- mle2(minuslogl=prop ~ dnorm(mean =  -1 * ((time / a)^ b), sd = z), 
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
      #CIs <- confint( profile(best.fit))
      print(coef(best.fit))
      # a
      summ[counter,3]=coef(best.fit)[1]
      # b
      summ[counter,4]=coef(best.fit)[2]
      # c
      summ[counter,5]=coef(best.fit)[3]
      #summ[counter,5]= round(max(repObs$logabund),1)
      # z
      summ[counter,6]=AIC(best.fit)
      #summ[counter,8]=CIs[1,1]
      #summ[counter,9]=CIs[1,2]
      #summ[counter,10]=CIs[2,1]
      #summ[counter,11]=CIs[2,2]
      #summ[counter,12]=CIs[3,1]
      #summ[counter,13]=CIs[3,2]
      #summ[counter,14]=CIs[4,1]
      #summ[counter,15]=CIs[4,2]
      summ[counter,7]=length(repObs$time)

      ### *** Comment/Uncomment following code to make pdf figs*** ###
      title=paste(strains[i],"  rep ",reps[j])
      plot(repObs$time,repObs$prop,main=title,ylim=c(min(repObs$prop),0), 
           xlab = 'Time (days)', ylab = 'Proportion surviving, log' )
      predTime=seq(0,max(repObs$time))
      print(strains[i])
      print(reps[j])
      #coef(best.fit)[3] 
      #lines(repObs$time, exp( -1 * ((repObs$time / coef(best.fit)[1] )^ coef(best.fit)[2])), 
      #        lwd=4, lty=2, col = "red")
      lines(repObs$time, (-1 * ((repObs$time / coef(best.fit)[1] )^ coef(best.fit)[2])), 
              lwd=4, lty=2, col = "red")
      counter=counter+1
    }
  }
}
  
dev.off() 
summ=summ[!is.na(summ[,1]),]
#colnames(summ)=c('strain','rep','a','b','c','z','AIC', 'a.CI.2.5', 'a.CI.97.5', 'b.CI.2.5', 'b.CI.97.5', 'c.CI.2.5', 'c.CI.97.5', 'z.CI.2.5', 'z.CI.97.5')
colnames(summ)=c('strain','rep','beta','alpha','std_dev','AIC', 'N.obs')
write.csv(summ,"data/demography/weibull_results.csv")