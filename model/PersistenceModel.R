# biological models for bacterial persistence during starvation
# 10/08/2014
# SEJ, JLT, KL

# clear workspace
rm(list=ls())

# load library for simulating
library(deSolve)

#population decay with constant death rate and no cannibalism --> only allows log-linear decline
persist<-function(t,x,parms){
	with(as.list(c(parms,x)),{
		dBdt=-B*d

		res=c(dBdt)
		list(res)
	})
}

# population decay with constant death rate and cannibalism --> only allows log-linear decline
persistCannibalism<-function(t,x,parms){
	with(as.list(c(parms,x)),{
		dBdt=B*umax*C/(C+k)*Yb-B*d
		dCdt=B*d*Yc-B*umax*C/(C+k)
		
		res=c(dBdt,dCdt)
		list(res)
	})
}

# population decay with constant death rate and death rate evolution, no cannibalism --> only allows log-non-linear, but monotonic, decline
persistEvolveD<-function(t,x,parms){
	with(as.list(c(parms,x)),{
		d=d0+dm*t
		dBdt=-B*d
		
		res=c(dBdt)
		list(res)
	})
}

# population decay with death rate evolution and cannibalism
persistEvolveDcannibalism<-function(t,x,parms){
	with(as.list(c(parms,x)),{
		d=d0+dm*t
		dBdt=B*umax*C/(C+k)*Yb-B*d
		dCdt=B*d*Yc-B*umax*C/(C+k)
		
		res=c(dBdt,dCdt)
		list(res)
	})
}

# population decay with constant death rate and cannibalism with evolving yield --> can show increase in population after initial decline
persistCannibalismEvolveYb<-function(t,x,parms){
	with(as.list(c(parms,x)),{
		Yb=Yb0+Ybm*t
		dBdt=B*umax*C/(C+k)*Yb-B*d
		dCdt=B*d*Yc-B*umax*C/(C+k)
		
		res=c(dBdt,dCdt)
		list(res)
	})
}

# population decay with evolving death rate and cannibalism with evolving yield --> increase after initial decline possible, but depends on parameterization
persistCannibalismEvolveYbEvolveD<-function(t,x,parms){
	with(as.list(c(parms,x)),{
		Yb=Yb0+Ybm*t
		d=d0+dm*t
		dBdt=B*umax*C/(C+k)*Yb-B*d
		dCdt=B*d*Yc-B*umax*C/(C+k)
		
		res=c(dBdt,dCdt)
		list(res)
	})
}


# parameter values
tsteps=150		# number of time steps to simulate
BGE=0.5			# bacterial growth efficiency
Yc=20e-9		# ug C per cell (currently assumes all goes to C at death); could make some slowly decay
Yb=1/Yc*BGE		# cells per ug C consumed
Yb0=1/Yc*BGE		# intercept cells per ug C consumed
Ybm=Yb0/100		# change of yield due to evolution yields a doubling over 100 days
#POTENTIAL ISSUE HERE: really this might change if cells get smaller (lower Yc) or if BGE improves. Right now Yc is independent of Yb...
umax=1.2e-8		# per capita maximum rate of uptake of C (ugC cell-1 day-1); E.coli 1.2e-8
k=1				# half saturation of carbon ugC L-1; E. coli 2.4
d=0.05			# death rate of cells
d0=0.05			# death rate intercept
dm=-0.0003		# death rate slope; NOTE- slope needs to be adjusted for longer simulations because you can generate positive death rates...


# simulate constant death rate and no cannibalism
parms=c(d=d)
xstart=c(B=1e9)
times=seq(0,tsteps,0.1)
outPersist=ode(y=xstart,times=times,func=persist,parms=parms)

dev.new()
plot(outPersist[,1],log10(outPersist[,2]),type='l',xlab="time",ylab="log Bacteria (cells)",ylim=c(5,9))

print(outPersist[nrow(outPersist),])


# simulate constant death rate and cannibalism
parms=c(Yc=Yc,Yb=Yb,umax=umax,k=k,d=d)
xstart=c(B=1e9,C=0)
times=seq(0,tsteps,0.1)
outPC=ode(y=xstart,times=times,func=persistCannibalism,parms=parms)

dev.new()
par(mfrow=c(2,1))
plot(outPC[,1],log10(outPC[,2]),type='l',xlab="time",ylab="log Bacteria (cells)",ylim=c(5,9))
plot(outPC[,1],outPC[,3],type='l',xlab="time",ylab="Carbon (ug L-1)")

print(outPC[nrow(outPC),])


# simulate constant death rate and death rate evolution, no cannibalism
parms=c(d0=d0,dm=dm)
xstart=c(B=1e9)
times=seq(0,tsteps,0.1)
outED=ode(y=xstart,times=times,func=persistEvolveD,parms=parms)

dev.new()
plot(outED[,1],log10(outED[,2]),type='l',xlab="time",ylab="log Bacteria (cells)",ylim=c(0,9))

print(outED[nrow(outED),])


# simulate evolving death rate and cannibalism
parms=c(d0=d0,dm=dm,Yc=Yc,Yb=Yb,umax=umax,k=k)
xstart=c(B=1e9,C=0)
times=seq(0,tsteps,0.1)
outEDC=ode(y=xstart,times=times,func=persistEvolveDcannibalism,parms=parms)

dev.new()
par(mfrow=c(2,1))
plot(outEDC[,1],log10(outEDC[,2]),type='l',xlab="time",ylab="log Bacteria (cells)",ylim=c(0,9))
plot(outEDC[,1],outEDC[,3],type='l',xlab='time',ylab="Carbon (ug L-1)")

print(outED[nrow(outED),])


# simulate constant death rate and cannibalism with evolving yield
parms=c(Yc=Yc,Yb0=Yb0,Ybm=Ybm,umax=umax,k=k,d=d)
xstart=c(B=1e9,C=0)
times=seq(0,tsteps,0.1)
outPCEY=ode(y=xstart,times=times,func=persistCannibalismEvolveYb,parms=parms)

dev.new()
par(mfrow=c(2,1))
plot(outPCEY[,1],log10(outPCEY[,2]),type='l',xlab="time",ylab="log Bacteria (cells)",ylim=c(0,9))
plot(outPCEY[,1],outPCEY[,3],type='l',xlab="time",ylab="Carbon (ug L-1)")

print(outPCEY[nrow(outPCEY),])


# simulate decay with evolving death rate and cannibalism with evolving yield
parms=c(Yc=Yc,Yb0=Yb0,Ybm=Ybm,umax=umax,k=k,d0=d0,dm=dm)
xstart=c(B=1e9,C=0)
times=seq(0,tsteps,0.1)
outPCEYED=ode(y=xstart,times=times,func=persistCannibalismEvolveYbEvolveD,parms=parms)

dev.new()
par(mfrow=c(2,1))
plot(outPCEYED[,1],log10(outPCEYED[,2]),type='l',xlab="time",ylab="log Bacteria (cells)",ylim=c(0,15))
plot(outPCEYED[,1],outPCEYED[,3],type='l',xlab="time",ylab="Carbon (ug L-1)")

print(outPCEYED[nrow(outPCEYED),])



# plot four simulation outcomes on one axis
dev.new()
plot(outPersist[,1],log10(outPersist[,2]),type='l',xlab="time",ylab="log Bacteria (cells)",ylim=c(5,9))
lines(outPC[,1],log10(outPC[,2]),col='red',lwd=2)
lines(outED[,1],log10(outED[,2]),col='green',lwd=2)
lines(outPCEY[,1],log10(outPCEY[,2]),col='blue',lwd=2)
lines(outPCEYED[,1],log10(outPCEYED[,2]),col='orange',lwd=2)
legend('bottomleft',c('simple decay','cannibalism','evolving death rate','cannibalism + evolving yield','cannibalism+evolving yield + evolving death rate'),lty=1,col=c('black','red','green','blue','orange'),box.lty=0)


# plot four simulation outcomes on one axis -- modifications by JTL for talk
png(filename="/Users/lennonj/Desktop/Persistence_Plot.png",width=1200, height=1200, res=96*2) 
dev.new()
par(bg="white",mar=c(5,5,1,1),lwd=3)
# plot four simulation outcomes on one axis
dev.new()
plot(outPersist[,1],log10(outPersist[,2]),type='l',lty=2, xlab="time",ylab="Cells / mL",ylim=c(5,9),xlim=c(-10,175),lwd=4,bty="n")
box(lwd=4)

lines(outPC[,1],log10(outPC[,2]),col='black',lwd=4) # scavenging
lines(outPCEYED[,1],log10(outPCEYED[,2]),col='black',lwd=4) # evolving


lines(outED[,1],log10(outED[,2]),col='green',lwd=2)
lines(outPCEY[,1],log10(outPCEY[,2]),col='blue',lwd=2)
lines(outPCEYED[,1],log10(outPCEYED[,2]),col='orange',lwd=2)
legend('bottomleft',c('simple decay','cannibalism','evolving death rate','cannibalism + evolving yield','cannibalism+evolving yield + evolving death rate'),lty=1,col=c('black','red','green','blue','orange'),box.lty=0)
