
# Gomp fits
# Initial parameters
#Z = Error
grids.gomp<-list(b=c(0.001, 0.005),eta=c(0.001),z=c(0.1,1,10))
start.gomp<-list(b=NA,eta=NA,z=NA)
grid.starts.gomp<-as.matrix(expand.grid(grids.gomp))
ncombos.gomp<-dim(grid.starts.gomp)[[1]]
# cycle through each combo
res.mat.gomp<-matrix(NA,nrow=ncombos.gomp,ncol=I(length(start.gomp)+1))
res.mod.gomp<-list()
for(k.gomp in 1:dim(grid.starts.gomp)[[1]]){
  #some how need to match grid parameters to start lists.
  mod.start.gomp<-as.list(grid.starts.gomp[k.gomp,])	
  new.start.gomp<-start.gomp
  new.start.gomp[names(start.gomp) %in% names(mod.start.gomp)]<-mod.start.gomp
  print(new.start.gomp)
  pscale.gomp<-as.numeric(new.start.gomp)
  names(pscale.gomp)<-names(new.start.gomp)
  print(pscale.gomp)
  
  fit.gomp <- mle2(minuslogl=prop ~ dnorm(mean = -1 * (eta/b) * (exp(b*time)-1), sd = z), 
                   start = new.start.gomp, data = repObs,# optimizer ="constrOptim",
                   control=list(parscale=pscale.gomp, maxit=1000), lower=c(b=0,eta=0), 
                   method="L-BFGS-B", hessian = T)
  print(summary(fit.gomp) )
  res.mat.gomp[k.gomp,]<-c(coef(fit.gomp), AIC(fit.gomp))		
  res.mat.gomp[[k.gomp]]<-fit.gomp
}


colnames(res.mat.gomp)<-c(names(coef(fit.gomp)),"AIC")
best.fit.gomp<-res.mod.gomp[[which(res.mat.gomp[,'AIC']==min(res.mat.gomp[,'AIC']))[1]]]      
print(best.fit.gomp)







perm.number <- 10000
null.dist <- replicate(perm.number, mean(sample(abs(ltde.met.ppca.PC1.L), 50, FALSE)))
ltde.met.ppca.PC1.L.null <- as.matrix(replicate(perm.number, sample(ltde.met.ppca.PC1.L, length(ltde.met.ppca.PC1.L), FALSE)))
ltde.met.ppca.PC1.L.null.sort <- as.matrix(t(apply(ltde.met.ppca.PC1.L.null, 1, sort)))
ltde.met.ppca.PC1.L.null.sort.merge <- cbind(ltde.met.ppca.PC1.L.null.sort, ltde.met.ppca.PC1.L)

perm.results <- matrix(ncol=2, nrow=length(ltde.met.ppca.PC1.L))
rownames(perm.results) <- names(ltde.met.ppca.PC1.L)
colnames(perm.results) <- c("corr", "p.value")
for(i in rownames(ltde.met.ppca.PC1.L.null.sort.merge)){
  row <- ltde.met.ppca.PC1.L.null.sort.merge[rownames(ltde.met.ppca.PC1.L.null.sort.merge) == i,] 
  #print(length(row[1:10000]))
  row.perm <- row[1:perm.number]
  row.slope <- row[perm.number + 1]
  if (row.slope > 0) {
    row.p <- (length(row.perm[row.perm > row.slope]) + 1 ) / (perm.number + 1)
  } else {
    row.p <- (length(row.perm[row.perm < row.slope]) + 1 ) / (perm.number + 1)
  }
  perm.results[i,] <- c(as.numeric(row.slope), as.numeric(row.p))
}

perm.results <- matrix(ncol=2, nrow=length(ltde.met.ppca.PC1.L))
rownames(perm.results) <- names(ltde.met.ppca.PC1.L)
colnames(perm.results) <- c("corr", "p.value")
#for(corr.i in ltde.met.ppca.PC1.L){
#if (corr.i > 0) {
#  corr.i.p <- length(null.dist[null.dist > corr.i]) / perm.number
#  print(corr.i.p)
#} else {
#  corr.i.p <- length(null.dist[null.dist < corr.i]) / perm.number
#}
#  print(abs(corr.i))
#  corr.i.p <- length(null.dist[null.dist > abs(corr.i)]) / perm.number
#  perm.results[i,] <- c(as.numeric(corr.i), as.numeric(corr.i.p))
#}

## ask nathan about this.....

perm.results <- data.frame(perm.results)
perm.results$p.value.BH <- p.adjust(perm.results$p.value, method="BH")

ltde.met.ppca.PC1.L > lower.CI
#lower.CI <- sort(boot)[ 10000 * 0.025 ]
#upper.CI <- sort(boot)[ 10000 * 0.975 ]
