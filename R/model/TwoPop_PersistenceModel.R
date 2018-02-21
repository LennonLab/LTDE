# Install package
require("deSolve")

# straight decay
persist<-function(t,x,parms){
	with(as.list(c(parms,x)),{
		dBadt=-Ba*da		#g (L*day)^-1
		
		res=c(dBadt)
		list(res)
	})
}

#parameters
# da= mortality rate [day^-1]

# assumptions
# 20e-15 g (femtograms of C per cell)
# start with 1e12 cells per L (0.02 g C in Ba)

parms = c(da = 0.1)
times = seq(0,1000,0.1)
init = c(Ba = 0.02)

outP=ode(y=init,times=times,func=persist,parms=parms)

dev.new()
plot(outP[,1],log10(outP[,2]/(20e-15)),xlab="time",ylab="log10 cells",type='l')


# population decay without evolution, but C recycling
persistCannibalism<-function(t,x,parms){
	with(as.list(c(parms,x)),{
		dBadt=Ba*(C*Va)/(C+Ka)*Ea-Ba*da		#g (L*day)^-1
		dCdt=Ba*da*m-Ba*(C*Va)/(C+Ka)		#g (L*day)^-1
		
		res=c(dBadt,dCdt)
		list(res)
	})
}

# parameters
# Va= mass-specific maximum uptake rate [g (g*day)^-1]
# Ka= carbon half-saturation constant [g L^-1]
# Ea= growth efficiency
# da= mortality rate [day^-1]
# m= fraction of dead bacteria available for consumption

# assumptions
# 20e-15 g (femtograms of C per cell)
# start with 1e12 cells per L (0.02 g C in Ba)

# Apsergillus niger??? Ka=15 uM glucose, Va=1 umol (gram dry weight)^-1 second^-1 Jorgensen et al. 2007 Microbiology, 153: 1963-1973
					#  Ka=0.0011 g C L^-1; Va=12.4 g C (g bacterial C)^-1 day^-1

#parms = c(Va = 10^-10, Ka = 10^-10, Ea = 0.75, da = 0.1, m=0.5) # these are from Jay's old script
parms = c(Va = 12, Ka = 0.0011, Ea = 0.75, da = 0.1, m=0.5)
times = seq(0,1000,0.1)
init = c(Ba = 0.02, C = 0)

outPC=ode(y=init,times=times,func=persistCannibalism,parms=parms)

dev.new()
par(mfrow=c(2,1))
plot(outPC[,1],log10(outPC[,2]/(20e-15)),xlab="time",ylab="log10 cells",type='l')
plot(outPC[,1],outPC[,3],type='l',xlab='time',ylab='carbon (g)')
## interesting prameter = time until extinction! (i.e., x-intercept)

# add mutant population
persistCannibalismMut<-function(t,x,parms){
	with(as.list(c(parms,x)),{
		dBadt=Ba*(C*Va)/(C+Ka)*Ea-Ba*da
		dBmdt=Bm*(C*Vm)/(C+Km)*Em-Bm*dm
		dCdt=Ba*da*m-Ba*(C*Va)/(C+Ka)+Bm*dm*m-Bm*(C*Vm)/(C+Km)
		
		res=c(dBadt,dBmdt,dCdt)
		list(res)
	})
}


# simulate with lower death rate in mutant; no recycling

# fraction of initial population that are mutant
fracMut=1e-5

parmsMut = c(Va = 0, Ka = 0.0011, Ea = 0.75, da = 0.1, Vm=0, Km=0.0011, Em=0.75, dm=0.05, m=0.5)
initMut = c(Ba = 0.02-0.02*fracMut, Bm=0.02*fracMut, C = 0)

outPCmut=ode(y=initMut,times=times,func=persistCannibalismMut,parms=parmsMut)

# summary plots of model with mutant dynamics
dev.new()
par(mfrow=c(2,3))
plot(outPCmut[,1],log10(rowSums(outPCmut[,2:3])/(20e-15)),xlab="time",ylab="log10 cells",type='l')
plot(outPCmut[,1],log10(outPCmut[,2]/(20e-15)),xlab="time",ylab="log10 ancestral cells",type='l')
plot(outPCmut[,1],log10(outPCmut[,3]/(20e-15)),xlab="time",ylab="log10 mutant cells",type='l')
plot(outPCmut[,1],outPCmut[,3]/rowSums(outPCmut[,2:3]),xlab="time",ylab="prop. cells",type='l')
plot(outPCmut[,1],outPCmut[,4],type='l',xlab='time',ylab='carbon (g)')



# simulate with lower death rate in mutant and recycling

# fraction of initial population that are mutant
fracMut=1e-5

parmsMut = c(Va = 12, Ka = 0.0011, Ea = 0.75, da = 0.1, Vm=12, Km=0.0011, Em=0.75, dm=0.05, m=0.5)
initMut = c(Ba = 0.02-0.02*fracMut, Bm=0.02*fracMut, C = 0)

outPCmut=ode(y=initMut,times=times,func=persistCannibalismMut,parms=parmsMut)

# summary plots of model with mutant dynamics
dev.new()
par(mfrow=c(2,3))
plot(outPCmut[,1],log10(rowSums(outPCmut[,2:3])/(20e-15)),xlab="time",ylab="log10 cells",type='l')
plot(outPCmut[,1],log10(outPCmut[,2]/(20e-15)),xlab="time",ylab="log10 ancestral cells",type='l')
plot(outPCmut[,1],log10(outPCmut[,3]/(20e-15)),xlab="time",ylab="log10 mutant cells",type='l')
plot(outPCmut[,1],outPCmut[,3]/rowSums(outPCmut[,2:3]),xlab="time",ylab="prop. cells",type='l')
plot(outPCmut[,1],outPCmut[,4],type='l',xlab='time',ylab='carbon (g)')


# simulate with higher efficiency in mutant and recycling

# fraction of initial population that are mutant
fracMut=1e-5

parmsMut = c(Va = 12, Ka = 0.0011, Ea = 0.5, da = 0.1, Vm=12, Km=0.0011, Em=0.9, dm=0.1, m=0.5)
initMut = c(Ba = 0.02-0.02*fracMut, Bm=0.02*fracMut, C = 0)

outPCmut=ode(y=initMut,times=times,func=persistCannibalismMut,parms=parmsMut)

# summary plots of model with mutant dynamics
dev.new()
par(mfrow=c(2,3))
plot(outPCmut[,1],log10(rowSums(outPCmut[,2:3])/(20e-15)),xlab="time",ylab="log10 cells",type='l')
plot(outPCmut[,1],log10(outPCmut[,2]/(20e-15)),xlab="time",ylab="log10 ancestral cells",type='l')
plot(outPCmut[,1],log10(outPCmut[,3]/(20e-15)),xlab="time",ylab="log10 mutant cells",type='l')
plot(outPCmut[,1],outPCmut[,3]/rowSums(outPCmut[,2:3]),xlab="time",ylab="prop. cells",type='l')
plot(outPCmut[,1],outPCmut[,4],type='l',xlab='time',ylab='carbon (g)')



# simulate with lower death rate & higher efficiency in mutant and recycling

# fraction of initial population that are mutant
fracMut=1e-5

parmsMut = c(Va = 12, Ka = 0.0011, Ea = 0.5, da = 0.1, Vm=12, Km=0.0011, Em=0.9, dm=0.05, m=0.5)
initMut = c(Ba = 0.02-0.02*fracMut, Bm=0.02*fracMut, C = 0)

outPCmut=ode(y=initMut,times=times,func=persistCannibalismMut,parms=parmsMut)

# summary plots of model with mutant dynamics
dev.new()
par(mfrow=c(2,3))
plot(outPCmut[,1],log10(rowSums(outPCmut[,2:3])/(20e-15)),xlab="time",ylab="log10 cells",type='l')
plot(outPCmut[,1],log10(outPCmut[,2]/(20e-15)),xlab="time",ylab="log10 ancestral cells",type='l')
plot(outPCmut[,1],log10(outPCmut[,3]/(20e-15)),xlab="time",ylab="log10 mutant cells",type='l')
plot(outPCmut[,1],outPCmut[,3]/rowSums(outPCmut[,2:3]),xlab="time",ylab="prop. cells",type='l')
plot(outPCmut[,1],outPCmut[,4],type='l',xlab='time',ylab='carbon (g)')


## what are key things we want to see with sims?  
# Change in frequency of genotypes with associted parameters
# Carbon should go down initially and then rate should slow down, but dynamics will depend on uptake parms
# Revisit my box model. 
