---
title: "Microbial death model"
author: "Stuart E. Jones and Jay T. Lennon"
date: "`r format(Sys.time(), '%d %B, %Y')`"
header-includes:
  - \usepackage{array}
output: pdf_document
geometry: margin=2.54cm
---

## OVERVIEW
  In January 2013, we started an experiment using ~24 strains of bacteria -- some "lab" strains others enviornmental strains -- to estimate survivorship during long-term starvation. Each strain was grown un into logarithmic phase and harvested. The cells were pelleted and washed 5x before being put in replicate 50 mL Falcon tubes. Subsamples were taken on overtime and plated onto R2A plates for enumeration as CFUs. After looking at plots of the data, in became apparent that CFUs were declining at a non-constant rate on semi-log plots. Here, we test whether or not the Weibull function can help explain variation in death curves.  




## 1) SET WORKING DIRECTORY, LOAD DATA, LOAD PACKAGE



```{r}  
  
  
  
  
  
  ## OVERVIEW
  In January 2013,

``{r}


In January 2013
# single model that is capable of simulating three scenarios
#	1. constant death with no scavenging and no dormancy
#	2. constant death with no dormancy
# 3. dormancy and scavenging

rm(list=ls())

# Install package
require("deSolve")
````
persist<-function(t,x,parms){
	with(as.list(c(parms,x)),{
		dAdt=A*(C*Va)/(C+Ka)*Ea-A*da-A*t2d
		dDdt=A*t2d-D*dd
		dCdt=A*da*m-A*(C*Va)/(C+Ka)+D*dd*m
		
		res=c(dAdt,dDdt,dCdt)
		list(res)
	})
}

times=1:200

##########################################################
# Scenario 1: simulate with no dormancy and no recycling #
##########################################################

# parameters and initial conditions
fracD=0
parmsNoDormNoScav = c(Va = 12, Ka = 0.0011, Ea = 0.75, da = 0.1, t2d = 0, dd=0.01, m=0)
initNoDormNoScav = c(A = 0.02-0.02*fracD, D=0.02*fracD, C = 0)

outNoDormNoScav=ode(y=initNoDormNoScav,times=times,func=persist,parms=parmsNoDormNoScav)

# summary plots of model with no dormancy and no recycling
dev.new()
par(mfrow=c(2,3))
plot(outNoDormNoScav[,1],log10(rowSums(outNoDormNoScav[,2:3])/(20e-15)),xlab="time",ylab="log10 cells",type='l')
plot(outNoDormNoScav[,1],log10(outNoDormNoScav[,2]/(20e-15)),xlab="time",ylab="log10 active cells",type='l')
if(any(outNoDormNoScav[,3]>0)){
	plot(outNoDormNoScav[,1],log10(outNoDormNoScav[,3]/(20e-15)),xlab="time",ylab="log10 dormant cells",type='l')
}else{
	plot(outNoDormNoScav[,1],rep(0,nrow(outNoDormNoScav)),xlab="time",ylab="dormant cells",type='l')
}
plot(outNoDormNoScav[,1],outNoDormNoScav[,2]/rowSums(outNoDormNoScav[,2:3]),xlab="time",ylab="prop. active",type='l')
plot(outNoDormNoScav[,1],outNoDormNoScav[,4],type='l',xlab='time',ylab='carbon (g)')


#######################################################
# Scenario 2: simulate with no dormancy but recycling #
#######################################################

# parameters and initial conditions
fracD=0
parmsNoDorm = c(Va = 12, Ka = 0.0011, Ea = 0.75, da = 0.1, t2d = 0, dd=0.01, m=0.5)
initNoDorm = c(A = 0.02-0.02*fracD, D=0.02*fracD, C = 0)

outNoDorm=ode(y=initNoDorm,times=times,func=persist,parms=parmsNoDorm)

# summary plots of model with no dormancy but with recycling
dev.new()
par(mfrow=c(2,3))
plot(outNoDorm[,1],log10(rowSums(outNoDorm[,2:3])/(20e-15)),xlab="time",ylab="log10 cells",type='l')
plot(outNoDorm[,1],log10(outNoDorm[,2]/(20e-15)),xlab="time",ylab="log10 active cells",type='l')
if(any(outNoDorm[,3]>0)){
	plot(outNoDorm[,1],log10(outNoDorm[,3]/(20e-15)),xlab="time",ylab="log10 dormant cells",type='l')
}else{
	plot(outNoDorm[,1],rep(0,nrow(outNoDorm)),xlab="time",ylab="dormant cells",type='l')
}
plot(outNoDorm[,1],outNoDorm[,2]/rowSums(outNoDorm[,2:3]),xlab="time",ylab="prop. active",type='l')
plot(outNoDorm[,1],outNoDorm[,4],type='l',xlab='time',ylab='carbon (g)')


####################################################
# Scenario 3: simulate with dormancy and recycling #
####################################################

# parameters and initial conditions
fracD=1e-5
parms = c(Va = 12, Ka = 0.0011, Ea = 0.75, da = 0.1, t2d = 0.0001, dd=0.01, m=0.5)
init = c(A = 0.02-0.02*fracD, D=0.02*fracD, C = 0)

out=ode(y=init,times=times,func=persist,parms=parms)

# summary plots of model with dormancy and recycling
dev.new()
par(mfrow=c(2,3))
plot(out[,1],log10(rowSums(out[,2:3])/(20e-15)),xlab="time",ylab="log10 cells",type='l')
plot(out[,1],log10(out[,2]/(20e-15)),xlab="time",ylab="log10 active cells",type='l')
if(any(out[,3]>0)){
	plot(out[,1],log10(out[,3]/(20e-15)),xlab="time",ylab="log10 dormant cells",type='l')
}else{
	plot(out[,1],rep(0,nrow(out)),xlab="time",ylab="dormant cells",type='l')
}
plot(out[,1],out[,2]/rowSums(out[,2:3]),xlab="time",ylab="prop. active",type='l')
plot(out[,1],out[,4],type='l',xlab='time',ylab='carbon (g)')


# plot with all three dynamics on same panel
dev.new()
plot(outNoDormNoScav[,1],log10(rowSums(outNoDormNoScav[,2:3])/(20e-15)),xlab="time",ylab="log10 cells",type='l',lwd=3)
lines(outNoDorm[,1],log10(rowSums(outNoDorm[,2:3])/(20e-15)),lty=2,lwd=3)
lines(out[,1],log10(rowSums(out[,2:3])/(20e-15)),col='grey',lwd=3)
legend('topright',c('no dormancy, no scavenging','no dormancy, scavenging','dormancy, scavenging'),lty=c(1,2,1),col=c('black','black','grey'),box.lty=0)