rm(list = ls())
getwd()
setwd("~/GitHub/LTDE/")

library('bbmle')
library('devtools')
library('plotrix')

## Load Data
obs <- read.csv("data/demography/longtermdormancy_20190528_nocomments.csv", 
                header = TRUE, stringsAsFactors = FALSE)
## Adding 1 to deal with log(0) observations
obs$Abund <- (as.numeric(obs$Colonies) +1)* (1000 / as.numeric(obs$Inoculum )) * ( 10 ^  as.numeric(obs$Dilution) )
strains <- sort(unique(obs$Strain))
#strains <- c('KBS0711')



obs <- obs[obs$Strain%in%strains,]
summ <- matrix(NA,length(strains)*max(obs$Rep),7)
pdf('figs/gompertz_fits.pdf') # Uncomment to create pdf that will plot data and fits
counter <- 1
for(i in 1:length(strains)){
  strainObs=obs[obs$Strain==strains[i],]
  reps=unique(strainObs$Rep)
  print(strains[i])
  for(j in 1:length(reps)){
    print(reps[j])
    repObs=strainObs[strainObs$Rep==reps[j],]
    # minimum of 10 data points
    if(nrow(repObs)>10){
      start=repObs[1,1]
      time<-(as.numeric(strptime(repObs$Firstread_date,format="%d-%b-%y",tz="EST"))-
               as.numeric(strptime(repObs$Dormstart_date,format="%d-%b-%y",tz="EST")))/(3600*24)
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
      repObs["prop"] <- log(abund / N_0)
      #repObs <- repObs[repObs$time < 50,]
      # Initial parameters
      grids<-list(b=c(0.1, 0.01, 0.001),eta=c(1, 1.5, 2),z=c(0.1,1,10))
      start<-list(b=NA,eta=NA,z=NA)
      grid.starts<-as.matrix(expand.grid(grids))
      ncombos<-dim(grid.starts)[[1]]
      # cycle through each combo
      res.mat<-matrix(NA,nrow=ncombos,ncol=I(length(start)+1))
      res.mod<-list()
      
      if(strains[i] == "ATCC13985"){
        b_lower <- 0.004
      } else if ( strains[i] == "ATCC43928") {
        b_lower <- 0.005
      } else if ( strains[i] == "KBS0701") {
        b_lower <- 0.01
      } else if ( strains[i] == "KBS0702") {
        b_lower <- 0.005
      } else if ( strains[i] == "KBS0703") {
        b_lower <- 0.006
      } else if ( strains[i] == "KBS0705") {
        b_lower <- 0.005
      } else if ( strains[i] == "KBS0706") {
        b_lower <- 0.006
      } else if ( strains[i] == "KBS0707") {
        b_lower <- 0.003
      } else if ( strains[i] == "KBS0710") {
        b_lower <- 0.005
      } else if ( strains[i] == "KBS0711") {
        b_lower <- 0.03
      } else if ( strains[i] == "KBS0712") {
        b_lower <- 0.005
      } else if ( strains[i] == "KBS0713") {
        b_lower <- 0.003
      } else if ( strains[i] == "KBS0714") {
        b_lower <- 0.02
      } else if ( strains[i] == "KBS0715") {
        b_lower <- 0.003
      } else if ( strains[i] == "KBS0721") {
        b_lower <- 0.003
      } else if ( strains[i] == "KBS0722") {
        b_lower <- 0.006
      } else if ( strains[i] == "KBS0724") {
        b_lower <- 0.004
      } else if ( strains[i] == "KBS0725") {
        b_lower <- 0.005
      } else if ( strains[i] == "KBS0727") {
        b_lower <- 0.006
      } else if ( strains[i] == "KBS0801") {
        b_lower <- 0.004
      } else if ( strains[i] == "KBS0802") {
        b_lower <- 0.003
      } else if ( strains[i] == "KBS0812") {
        b_lower <- 0.003
      } else {
        b_lower <- 0.005
      }
       
      for(k in 1:dim(grid.starts)[[1]]){
        #some how need to match grid parameters to start lists.
        mod.start<-as.list(grid.starts[k,])	
        new.start<-start
        new.start[names(start) %in% names(mod.start)]<-mod.start
        pscale<-as.numeric(new.start)
        names(pscale)<-names(new.start)
        
        fit <- mle2(minuslogl=prop ~ dnorm(mean = -1 * eta * (exp(b*time)-1), sd = z), 
                         start = new.start, data = repObs,# optimizer ="constrOptim",
                         control=list(parscale=pscale, maxit=1000), lower=c(b=b_lower,eta=0.01,z=0.01), upper=c(b=100,eta=1000,z=100),
                         method="L-BFGS-B", hessian = T)
       
        res.mat[k,]<-c(coef(fit), AIC(fit))		
        res.mod[[k]]<-fit
      }
      
      colnames(res.mat)<-c(names(coef(fit)),"AIC")
      best.fit<-res.mod[[which(res.mat[,'AIC']==min(res.mat[,'AIC']))[1]]]

      summ[counter,1]=strains[i]
      summ[counter,2]=reps[j]
      # b
      b <- coef(best.fit)[1]
      summ[counter,3]=b
      # eta
      eta <- coef(best.fit)[2]
      summ[counter,4]=eta
      # z
      summ[counter,5]=coef(best.fit)[3]
      summ[counter,6]=AIC(best.fit)
      summ[counter,7]=length(repObs$time)
      #CIs <- confint(profile(best.fit))

      ### *** Comment/Uncomment following code to make pdf figs*** ###
      title=paste(strains[i],"  rep ",reps[j])
      plot(repObs$time,repObs$prop,main=title,ylim=c(min(repObs$prop),0), 
           xlab = 'Time (days)', ylab = 'Proportion surviving, log' )
      predTime=seq(0,max(repObs$time))
      lines(repObs$time, -1 * eta * (exp(b*time)-1), 
            lwd=4, lty=2, col = "red")
      counter=counter+1
    }
  }
}
  
dev.off() 
print(summ)
summ=summ[!is.na(summ[,1]),]
colnames(summ)=c('strain','rep','b','eta','std_dev','AIC', 'N.obs')
write.csv(summ,"data/demography/gompertz_results.csv")
  
  
# clean the results file
df <- read.table("data/demography/gompertz_results.csv", 
                 header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
df <- df[!(df$strain=="KBS0711W"),]
# remove KBS0711 replicates 10, 11, and 12
# these samples were only sampled starting at day 100
# We don't have time to re-do them, so just remove them
df <- df[!(df$strain == "KBS0711" & df$rep == 10 ),] 
df <- df[!(df$strain == "KBS0711" & df$rep == 11 ),] 
df <- df[!(df$strain == "KBS0711" & df$rep == 12 ),] 

write.csv(df, file = "data/demography/gompertz_results_clean.csv")
  
  
