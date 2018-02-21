## Retrieve and Set Your Working Directory
rm(list = ls())
getwd()
setwd("~/GitHub/Dimensions/Aim1/DeathCurves/")

## Load Data
obs <- read.csv("longtermdormancy_20150526b.csv", header = TRUE, stringsAsFactors = FALSE)

## Estimating CFUs
## Adding 1 to deal with log(0) observations --> should we just remove instead?

obs$Abund <- as.numeric(obs$Colonies) * 10 ^ as.numeric(obs$Dilution) + 1

# likelihood function for linear fit to log transformed data
fitLogLinearDecay<-function(p,N,time){
	k=p[1]
	N0=p[2]
	sd=exp(p[3])
	
	
	Nt=N0-k*time
	
	-sum(dnorm(N,Nt,sd,log=TRUE))
}

# likelihood function for quadratic fit to log transformed data
fitLogQuadDecay<-function(p,N,time){
	k2=p[1]
	k=p[2]
	N0=p[3]
	sd=exp(p[4])

	Nt=N0-k2*time^2-k*time
	
	-sum(dnorm(N,Nt,sd,log=TRUE))
}


strains=sort(unique(obs$Strain))
#only fit strains with more than 10 observations
strains=strains[table(obs$Strain)>10]
obs=obs[obs$Strain%in%strains,]

#pdf('decayFits.pdf')

# matrix for storing model fit output
summ=matrix(NA,length(strains)*max(obs$Rep),10)
counter=1

for(i in 1:length(strains)){
	strainObs=obs[obs$Strain==strains[i],]
	
	reps=unique(strainObs$Rep)
	for(j in 1:length(reps)){
		
		repObs=strainObs[strainObs$Rep==reps[j],]

		if(nrow(repObs)>3){
			
			start=repObs[1,1]
		
			time=(as.numeric(strptime(repObs$Firstread_date,format="%d-%b-%y",tz="EST"))-as.numeric(strptime(start,format="%d-%b-%y",tz="EST")))/(3600*24)
		
			logLinCur=optim(c((log10(repObs$Abund[1])-log10(repObs$Abund[nrow(repObs)]))/(time[length(time)]-time[1]),log10(max(repObs$Abund)),1),fitLogLinearDecay,N=log10(repObs$Abund),time=time)
			tempQuad=optim(c(((log10(repObs$Abund[1])-log10(repObs$Abund[nrow(repObs)]))/(time[length(time)]-time[1]))/10,(log10(repObs$Abund[1])-log10(repObs$Abund[nrow(repObs)]))/(time[length(time)]-time[1]),log10(max(repObs$Abund)),1),fitLogQuadDecay,N=log10(repObs$Abund),time=time)
			tempQuad2=optim(c(-5e-6,0.005,5,log(0.1)),fitLogQuadDecay,N=log10(repObs$Abund),time=time)
		
			if(tempQuad2$value<=tempQuad$value){
				logQuadCur=tempQuad2
			}else{
				logQuadCur=tempQuad	
			}
		
			summ[counter,1]=strains[i]
			summ[counter,2]=reps[j]
			summ[counter,3:4]=logLinCur$par[1:2]
			summ[counter,5]=round(2*logLinCur$value+2*length(logLinCur$par),2)
			summ[counter,6:8]=logQuadCur$par[1:3]
			summ[counter,9]=round(2*logQuadCur$value+2*length(logQuadCur$par),2)
		
			# run likelihood ratio test to compare the linear and quadratic model; store p-value in summ
			summ[counter,10]=pchisq((2*logLinCur$value-2*logQuadCur$value),df=1,lower.tail=FALSE)
		
			#### add indicator of non-linearity (****) to plot title if quadratic model is better
#			if(!is.na(summ[counter,10])){
#				if(as.numeric(summ[counter,10])<0.05){
#					title=paste(strains[i]," rep ",reps[j],"****")
#				}else{
#					title=paste(strains[i]," rep ",reps[j])
#				}
#			}else{
#				title=paste(strains[i]," rep ",reps[j])
#			}
			
#			plot(time,log10(repObs$Abund),main=title,ylim=c(0,9))
#			predTime=seq(0,max(time))
#			lines(predTime,logLinCur$par[2]-logLinCur$par[1]*predTime,lwd=2,lty=2)
#			lines(predTime,logQuadCur$par[3]-logQuadCur$par[1]*predTime^2-logQuadCur$par[2]*predTime,col='red',lwd=2,lty=2)
		
			counter=counter+1
		}
	}
}

#dev.off()

summ=summ[!is.na(summ[,1]),]
colnames(summ)=c('strain','rep','linearfit_K','linearfit_N0','linearfit_AIC','quadfit_k2','quadfit_k','quadfit_N0','quadfit_AIC','LRT_pvalue')


#percent of tubes that are better fit by quadratic (AIC_linear > AIC_quadratic)
sum(as.numeric(summ[,5])>as.numeric(summ[,9]))/nrow(summ)*100	# 51.3%

# percent of tubes that are better based on likelihood ratio test
sum(as.numeric(summ[,10])<0.05)/nrow(summ)*100	# 46%

# add bonferroni correction
0.05/nrow(summ)	#  0.00044
sum(as.numeric(summ[,10])<0.0044)/nrow(summ)*100	# 31.9%

## histogram of likelihood ratio p-values
hist(log10(as.numeric(summ[,10])),breaks=seq(-20,0,1))
abline(v=log10(0.05),lwd=2,lty=2,col='red')
abline(v=log10(0.00044),lwd=2,lty=2,col='green')
legend('topleft',c('alpha=0.05','alpha=0.00044 (bonferroni)'),lty=2,col=c('red','green'),box.lty=0,lwd=2)

### confirm that likelihood ratio test is showing what deltaAIC between linear and quadratic would suggest
plot(log10(as.numeric(summ[,10])),as.numeric(summ[,9])-as.numeric(summ[,5]),xlab="log10(LRT pvalue)",ylab="deltaAIC (quadratic - linear)")
abline(v=log10(0.05),lwd=2,lty=2,col='red')
abline(v=log10(0.00044),lwd=2,lty=2,col='green')

### look at agreement amongst reps
repAgreement=matrix(NA,length(strains),4)
repAgreement[,1]=strains
for(i in 1:length(strains)){
	cur=summ[summ[,1]==strains[i],]
	
	repAgreement[i,2]=nrow(cur)
	repAgreement[i,3]=sum(as.numeric(cur[,10])>0.05)
	repAgreement[i,4]=nrow(cur)-sum(as.numeric(cur[,10])>0.05)
}

repAgree=data.frame(strain=repAgreement[,1],Nreps=as.numeric(repAgreement[,2]),Nlinear=as.numeric(repAgreement[,3]),Nquad=as.numeric(repAgreement[,4]),stringsAsFactors=FALSE)


### generate trait output table
# when quadratic is not significantly better by likelihood ratio test report "evolvability" of 0
traitsReport=matrix(NA,nrow(summ),5)
traitsReport[,1]=summ[,1]
traitsReport[,2]=summ[,2]
for(i in 1:nrow(summ)){
	if(as.numeric(summ[i,10])<(0.05/nrow(summ))){
		traitsReport[i,3]=summ[i,7]
		traitsReport[i,4]=summ[i,6]
		traitsReport[i,5]=summ[i,8]
	}else{
		traitsReport[i,3]=summ[i,3]
		traitsReport[i,4]=0
		traitsReport[i,5]=summ[i,4]
	}
}
colnames(traitsReport)=c('strain','rep','decay','evolvability','N0')

# remove the extremely bad fit
traitsReport=traitsReport[-27,]

strains=sort(unique(traitsReport[,1]))
meanDecay=tapply(as.numeric(traitsReport[,3]),traitsReport[,1],FUN=mean)
meanEvolve=tapply(as.numeric(traitsReport[,4]),traitsReport[,1],FUN=mean)
meanN0=tapply(as.numeric(traitsReport[,5]),traitsReport[,1],FUN=mean)

summaryTraitsReport=data.frame(strains=strains,decay=meanDecay,evolvability=meanEvolve,N0=meanN0,stringsAsFactors=FALSE)

write.table(traitsReport,"perRepDeathCurveTraits.txt",row.names=FALSE,quote=FALSE)
write.table(summaryTraitsReport,"DeathCurveTraits.txt",row.names=FALSE,quote=FALSE)