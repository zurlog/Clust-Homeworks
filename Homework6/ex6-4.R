### EXERCISES 6-4 - MS&BDA 
### ZURLO GIOVANNI 21/11/2021

library(mclust)
library(cluster)
library(mclust)
library(smacof)
library(factoextra)
library(dbscan)
library(knitr)
library(ggdendro)
library(fpc)
library(EMMIXskew)

## EXERCISE 1 ####
setwd("C:/Users/giova/Dropbox/Magistrale/[NEW] MS & BD Analysis/Datasets")
stars<- read.table("stars5000.dat", header=TRUE)
str(stars)
library(clusterSim)
# positional standardization ((x-median)/mad)
sstars=data.Normalization(stars,type="n2",normalization="column")
pairs(sstars)

bicvals.sn<-bicvals.st<-bicvals.td<-bicvals.nd<-numeric(0)

## GAUSSIAN MIX ####
mclust=Mclust(sstars,G=1:20,verbose = T)
summary(mclust) # For sure a VVV model
pairs(sstars, pch=20,cex=0.1)
mclust$BIC[,14] # VVV BIC sequence
plot(mclust$BIC[-(1:6),14],type='b') # Best is 17 - BIC -51820.40



set.seed(1234)
ndist <- list()

for (i in 1:20){
  print(i)
  tryattempts <- 3
  trycounter <- 1
  tst <- try(ndist[[i]] <- EmSkew(sstars,g=i,distr="mvn",ncov=3))
  while((is.null(tst) | class(tst)=="try-error") & trycounter<tryattempts+1){
    print("Error, try again")
    tst <- try(ndist[[i]] <- EmSkew(sstars,g=1,distr="mvn",ncov=3))
    trycounter <- trycounter+1
  }
  trycounter <- 1
  while((is.null(tst) | class(tst)=="try-error") & trycounter<tryattempts+1){
    print("Error, try again")
    tst <- try(ndist[[i]] <- EmSkew(sstars,g=i,distr="mvn",ncov=4))
    trycounter <- trycounter+1
  }
  trycounter <- 1
  while((is.null(tst) | class(tst)=="try-error") & trycounter<tryattempts+1){
    print("Error, try again")
    tst <- try(ndist[[i]] <- EmSkew(sstars,g=i,distr="mvn",ncov=2))
    trycounter <- trycounter+1
  }
  bicvals.nd[i] <- ndist[[i]]$bic
  rm(tst)
}

bicvals.nd1 = c(110475.07,69536.17, 66641.99,  60541.57, 58248.76,  54870.09,
               56888.23,  54822.53, 54314.79, 56158.90,  53335.61,  53038.59,
               53247.56,  52475.27,  52430.73,  52087.31, 51933.18,
               51700.16,  52134.60,  52324.85)
plot(bicvals.nd,type='b') # Best is 18 - BIC 51700.16
which.min(bicvals.nd)

## T DISTRIBUTION ####
set.seed(1234)
tdist <- list()

for (i in 1:20){
  print(i)
  tryattempts <- 3
  trycounter <- 1
  tst <- try(tdist[[i]] <- EmSkew(sstars,g=i,distr="mvt",ncov=3))
  while((is.null(tst) | class(tst)=="try-error") & trycounter<tryattempts+1){
    print("Error, try again")
    tst <- try(tdist[[i]] <- EmSkew(sstars,g=1,distr="mvt",ncov=3))
    trycounter <- trycounter+1
  }
  trycounter <- 1
  while((is.null(tst) | class(tst)=="try-error") & trycounter<tryattempts+1){
    print("Error, try again")
    tst <- try(tdist[[i]] <- EmSkew(sstars,g=i,distr="mvt",ncov=4))
    trycounter <- trycounter+1
  }
  trycounter <- 1
  while((is.null(tst) | class(tst)=="try-error") & trycounter<tryattempts+1){
    print("Error, try again")
    tst <- try(tdist[[i]] <- EmSkew(sstars,g=i,distr="mvt",ncov=2))
    trycounter <- trycounter+1
  }
  bicvals.td[i] <- tdist[[i]]$bic
  rm(tst)
}

# MINIMUM VALUE FOR K=11 - BIC 52293.52
bicvals.td1=c(72171.13,64608.79, 60915.28, 58730.46, 56550.67, 54525.43, 53648.07,
              52898.45, 53096.22, 53264.08, 53195.38, 52721.67, 52612.79, 52645.55,
              52510.79, 52565.42, 52437.13, 51993.30, 52785.95, 53001.69)
plot(bicvals.td,type='b') # Best is 18 - BIC 51700.16
which.min(bicvals.td)


## SKEW NORMALS ####
set.seed(1234)
skew.norm <- list()

for (i in 1:15){
  print(i)
  tryattempts <- 3
  trycounter <- 1
  tst <- try(skew.norm[[i]] <- EmSkew(sstars,g=i,distr="msn",ncov=3,itmax = 5000))
  while((is.null(tst) | class(tst)=="try-error") & trycounter<tryattempts+1){
    print("Error, try again")
    tst <- try(skew.norm[[i]] <- EmSkew(sstars,g=1,distr="msn",ncov=3,itmax = 5000))
    trycounter <- trycounter+1
  }
  trycounter <- 1
  while((is.null(tst) | class(tst)=="try-error") & trycounter<tryattempts+1){
    print("Error, try again")
    tst <- try(skew.norm[[i]] <- EmSkew(sstars,g=i,distr="msn",ncov=4,itmax = 5000))
    trycounter <- trycounter+1
  }
  trycounter <- 1
  while((is.null(tst) | class(tst)=="try-error") & trycounter<tryattempts+1){
    print("Error, try again")
    tst <- try(skew.norm[[i]] <- EmSkew(sstars,g=i,distr="msn",ncov=2,itmax = 5000))
    trycounter <- trycounter+1
  }
  bicvals.sn[i] <- skew.norm[[i]]$bic
  rm(tst)
}

bicvals.sn1=c(104041.15, 67236.72,  62314.28, 59576.50, 58645.87, 57313.79,
              56669.74,  56590.62, 53148.74, 104041.15, 52611.90, 104041.15,
              52661.83, 104041.15, 104041.15)
plot(bicvals.sn,type='b')
which.min(bicvals.sn)
# MINIMUM VALUE FOR K=11 - BIC 52611.9

## SKEW T DISTR. ####
set.seed(1234)
skew.t <- list()

for (i in 1:15){
  print(i)
  tryattempts <- 3
  trycounter <- 1
  tst <- try(skew.t[[i]] <- EmSkew(sstars,g=i,distr="mst",ncov=3,itmax = 5000))
  while((is.null(tst) | class(tst)=="try-error") & trycounter<tryattempts+1){
    print("Error, try again")
    tst <- try(skew.t[[i]] <- EmSkew(sstars,g=1,distr="mst",ncov=3,itmax = 5000))
    trycounter <- trycounter+1
  }
  trycounter <- 1
  while((is.null(tst) | class(tst)=="try-error") & trycounter<tryattempts+1){
    print("Error, try again")
    tst <- try(skew.t[[i]] <- EmSkew(sstars,g=i,distr="mst",ncov=4,itmax = 5000))
    trycounter <- trycounter+1
  }
  trycounter <- 1
  while((is.null(tst) | class(tst)=="try-error") & trycounter<tryattempts+1){
    print("Error, try again")
    tst <- try(skew.t[[i]] <- EmSkew(sstars,g=i,distr="mst",ncov=2,itmax = 5000))
    trycounter <- trycounter+1
  }
  bicvals.st[i] <- skew.t[[i]]$bic
  rm(tst)
}

bicvals.st1=c(69241.32, 62566.32, 59418.47, 55750.76, 54844.15, 53892.91, 52848.98,
              53202.09, 53108.56, 52845.28, 52584.94, 52086.00, 52745.76, 52838.98,
              69241.32)
plot(bicvals.st,type='b') #convergence error for K=15
which.min(bicvals.st)
# MINIMUM VALUE FOR K=12 - BIC 52086


## COMPARISON ####
BestBICs=c(min(bicvals.nd),min(bicvals.td),min(bicvals.sn),min(bicvals.st))
BestCompN=c(which.min(bicvals.nd),which.min(bicvals.td),which.min(bicvals.sn),which.min(bicvals.st))


plot(1:20,bicvals.nd,typ="l",xlab="Number of Mixture Components",ylab="BIC",lty=2,col='steelblue1', lwd=2)
lines(1:20,bicvals.td,typ="l",col='firebrick1', lwd=2.3)
lines(1:15,bicvals.sn,typ="l",col='royalblue4',lty=2, lwd=2)
lines(1:15,bicvals.st,typ="l",col='tan2', lwd=2.3)
legend(x=4,y=90000,c("Normal","t","Skew-Normal","Skew-t"),col=c('steelblue1','firebrick1','royalblue4','tan2'),
       lty = c(2,1,2,1), cex=0.8, lwd=3)

plot(5:20,bicvals.nd[-(1:4)],typ="l",xlab="Number of Mixture Components",ylab="BIC",lty=2,col='steelblue1', lwd=2)
lines(5:20,bicvals.td[-(1:4)],typ="l",col='firebrick1', lwd=2.3)
lines(5:15,bicvals.sn[-(1:4)],typ="l",col='royalblue4',lty=2, lwd=2)
lines(5:15,bicvals.st[-(1:4)],typ="l",col='tan2', lwd=2.3)
legend(x=16,y=57000,c("Normal","t","Skew-Normal","Skew-t"),col=c('steelblue1','firebrick1','royalblue4','tan2'),
       lty = c(2,1,2,1), cex=0.8, lwd=3)

# Further Compare Skew-norm (13), Normal(18), t(18)
library(fpc)
pairs(sstars, col=skew.t[[12]]$clust, pch=20, cex=0.4,oma=c(2,2,2,2))
pairs(sstars, col=skew.t[[12]]$clust, pch=1, cex=0.7,oma=c(2,2,2,2))
pairs(sstars[,1:3], col=skew.norm[[13]]$clust, pch=clusym[skew.norm[[13]]$clust], cex=0.8,oma=c(2,2,2,2))
pairs(sstars[,4:6], col=skew.norm[[13]]$clust, pch=clusym[skew.norm[[13]]$clust], cex=0.8,oma=c(2,2,2,2))

pairs(sstars, col=tdist[[18]]$clust, pch=20, cex=0.4,oma=c(2,2,2,2))
pairs(sstars, col=tdist[[18]]$clust, pch=1, cex=0.7,oma=c(2,2,2,2))
pairs(sstars[,1:3], col=tdist[[18]]$clust, pch=clusym[tdist[[18]]$clust], cex=0.8,oma=c(2,2,2,2)) #hard to see
pairs(sstars[,4:6], col=tdist[[18]]$clust, pch=clusym[tdist[[18]]$clust], cex=0.8,oma=c(2,2,2,2))

pca=princomp(sstars)
summary(pca)

par(mfrow=c(1,3))
plot(x=pca$scores[,1],y=pca$scores[,2], cex=1, lwd=1.7,col=skew.norm[[13]]$clust,
     pch=clusym[skew.norm[[13]]$clust], xlab="PC1",ylab="", main = "Skew-t")
plot(x=pca$scores[,1],y=pca$scores[,2], cex=1, lwd=1.7,col=tdist[[18]]$clust,
     pch=clusym[tdist[[18]]$clust], xlab="PC1",ylab="", main = "T-Dist")
# By means of tdist mixture, extreme values all fall in cluster 1 (39)
pairs(sstars[,2:6], col=tdist[[18]]$clust, pch=clusym[tdist[[18]]$clust], cex=0.8,oma=c(2,2,2,2))
# Does not happen with mixture of skew-ts (splitted between 2,9,6, smallest. inside not overclustering)
pairs(sstars[,2:6], col=skew.norm[[13]]$clust, pch=clusym[skew.norm[[13]]$clust], cex=0.8,oma=c(2,2,2,2))

pairs(sstars[,2:6], col=ndist[[18]]$clust, pch=clusym[ndist[[18]]$clust], cex=0.8,oma=c(2,2,2,2))

for (i in 1:15){
errori[i] <- skew.t[[i]]$error
}
errori
for (i in 1:20){
errori[i] <- ndist[[i]]$error
}
errori
for (i in 1:20){
errori[i] <- tdist[[i]]$error
}
errori




## EXERCISE 2 ####
setwd("C:/Users/giova/Dropbox/Magistrale/[NEW] MS & BD Analytics/Datasets")
wdbc <- read.csv("wdbc.dat",header=FALSE)
data <- wdbc[,3:12]
dim(data)
p = 10
K = 4

VVV =Mclust(data, G=K, modelNames = 'VVV')
VVV$df
# (K-1) + K(p + p(p+1)/2)
((K-1) + K*(p + p*(p+1)/2)) # 263

VII =Mclust(data, G=K, modelNames = 'VII')
VII$df
# (K-1) + K(p + p(p+1)/2)
((K-1) + K*(p + 1)) #47

EEE =Mclust(data, G=K, modelNames = 'EEE')
EEE$df
# (K-1) + K(p + p(p+1)/2)
((K-1) + K*p + p*(p+1)/2) # 98

SkNorm=EmSkew(data, g=K, distr = 'msn', itmax = 5000, nhclust = T)
((K-1) + K*(p + p*(p+1)/2)) + K*p  # 303

Tdist=EmSkew(data, g=K, distr = 'mvt', itmax = 5000)
((K-1) + K*(p + p*(p+1)/2)) + K  # 267

SkewT=EmSkew(data, g=2, distr = 'mst')
((K-1) + K*(p + p*(p+1)/2)) + K*p + K  # 307

SkewT=EmSkew(data, g=2, distr = 'mst')
((K-1) + K*p + p*(p+1)/2 + p + 1)  # 109







## EXERCISE 4 ####
library(mclust)
setwd("C:/Users/giova/Dropbox/Magistrale/[NEW] MS & BD Analysis/Datasets")
stars<- read.table("stars5000.dat", header=TRUE)
str(stars)
library(clusterSim)
# positional standardization ((x-median)/mad)
sstars=data.Normalization(stars,type="n2",normalization="column")


library(parallel)
options(mc.cores = parallel::detectCores())
set.seed(1234)
index=sample(1:5000,1000)
rsubs=sstars[index,]
extsubs=sstars[-index,]
str(rsubs)
fit=Mclust(rsubs, G=1:15, verbose = T)
summary(fit)
# Some clusters leaved out
full=Mclust(sstars, G=1:15, verbose = T)
adjustedRandIndex(fit$classification,full$classification)


MclustBD=function(x,ns=2000,seed=1234, showProgress=T,...){
  # The ... arguments are passed on to the Mclust function
  # Setting a seed is required for reproducibility
  # Some default values are specified for seed and subset size
  start_time <- Sys.time()
  set.seed(seed)
  index=sample(1:nrow(x), size=ns)
  rsubs=x[index,]
  extsubs=x[-index,]
  fit=Mclust(rsubs,verbose = showProgress,...)
  pred=predict.Mclust(fit,extsubs)
  end_time <- Sys.time()
  print(paste("Running Time: ", round(end_time - start_time,3)),quote = F)
  print(summary(fit))
  
  # Returning a list with comprehensive results 
  out=list()
  out$call=fit$call
  out$summary=summary(fit)
  out$modelName=fit$modelName
  out$train_n=fit$n
  out$d=fit$d
  out$G=fit$G
  out$train_data=fit$data
  out$ext_data=as.matrix(extsubs)
  out$BIC=fit$BIC
  out$loglik=fit$loglik
  out$df=fit$df
  out$parameters=fit$parameters
  out$uncertainty=fit$uncertainty
  out$classification=integer(nrow(x))
  out$classification[index]=fit$classification
  out$classification[-index]=pred$classification
  out$probability=matrix(nrow = nrow(x),ncol = fit$G)
  out$probability[index,]=fit$z
  out$probability[-index,]=pred$z
  return(out)
}

system.time(Mclust(sstars, verbose = T)) # 28.71
start_time <- Sys.time()
Mclust(sstars, verbose = T)
end_time <- Sys.time()
print(paste("Running Time: ", round(end_time - start_time,3))) #28.22

## Optimal model ####
subs.fit=MclustBD(sstars,G=1:15,ns=1000)
full=Mclust(sstars, G=1:15, verbose = T)
adjustedRandIndex(subs.fit$classification,full$classification) # 0.40

library(knitr)
library(kableExtra)
table(subs.fit$classification,full$classification) %>%
  kable("html", caption = 'Optimal Mixture   -   Subset vs. Full', row.names = T) %>%
  kable_styling(full_width = F, position = "center")

#####


## VVV (7) model ####
subs.fit7=MclustBD(sstars,G=7,ns=1000,modelNames = 'VVV')
full7=Mclust(sstars, G=7, verbose = T,modelNames = 'VVV')
adjustedRandIndex(subs.fit7$classification,full7$classification) # 0.427

library(knitr)
library(kableExtra)
table(subs.fit7$classification,full7$classification) %>%
  kable("html", caption = 'VVV (7) Model   -   Subset vs. Full', row.names = T) %>%
  kable_styling(full_width = F, position = "float_left")

plot(full7, what = "density")
plot(subs.full7, what = "density")

#####

MclustBD <-function(x, ns = 2000, seed = 1234, showProgress = F, Summary = T,...){
  # The ... arguments are passed on to the Mclust function
  # Setting a seed is required for reproducibility
  # Some default values are specified for seed and subset size
  require(mclust)
  start_time <- Sys.time()
  set.seed(seed)
  # Draw a random subset of ns observations
  index <- sample(1:nrow(x), size=ns)
  rsubs <- x[index,]
  extsubs <- x[-index,]
  # Fitting Mclust to the subset
  fit=Mclust(rsubs,verbose = showProgress,...)
  # Extending the fitted model to all observations
  pred <- predict.Mclust(fit,extsubs)
  end_time <- Sys.time()
  # Printing Running Time & Summary
  if(showProgress == T){
    print(paste("Running Time: ", round(end_time - start_time,3)),quote = F)}
  
  # Returning a list with comprehensive results 
  out <- list()
  out$classification <- integer(nrow(x))
  out$classification[index] <- fit$classification
  out$classification[-index] <- pred$classification
  out$probability <- matrix(nrow = nrow(x),ncol = fit$G)
  out$probability[index,] <- fit$z
  out$probability[-index,] <- pred$z
  out$fit <- fit
  
  if(Summary == T){
    print(summary(fit))
    print(paste("* Overall Predictions (Out-of-Sample = ",nrow(extsubs), ") *", quote=F))
    t(table(out$classification))}
  
  return(out)}
