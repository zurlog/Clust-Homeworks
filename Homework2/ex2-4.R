### EXERCISES 2-4 - MS&BDA 
### ZURLO GIOVANNI 11/10/2021

##  EXERCISE 1 #### 

library(pdfCluster)
data("oliveoil")
diag(var(oliveoil))
solive=scale(oliveoil[,3:10])

library(cluster)
# CACHE THIS!
set.seed(1234)
cg.solive=clusGap(solive,kmeans, K.max=10, B=150, d.power=2,
                  spaceH0 = 'scaledPCA', nstart=100, iter.max=100)
plot(cg.solive, main="Gap Statistic Plot for std Oliveoil Data")
print(cg.solive,method = 'globalSEmax',SE.factor = 2)
print(cg.solive,method = 'firstSEmax', SE.factor = 2)
print(cg.solive,method = 'Tibs2001SEmax', SE.factor = 3)
maxSE(cg.solive$Tab[,3],cg.solive$Tab[,4],"globalSEmax",SE.factor=2 )
# Scree plot shows the elbow at k=5, but the gap follows his growth at constant rate
# Every criterion with se.factor = 2 points towards k=9 with the expection of tibs
# I would stick with k=9 since, as we saw in the previous exercises, this choice was
# also suggested by the precise matching with regions


artdata=read.table("clusterdata2.dat")
sartdata=scale(artdata)
library(cluster)
set.seed(1234)
cg.sart=clusGap(sartdata,kmeans, K.max=10, B=130, d.power=2,
                  spaceH0 = 'scaledPCA', nstart=100, iter.max=100)
plot(cg.sart, main="Gap Statistic Plot for std Artificial Data 2")
print(cg.sart,method = 'globalSEmax',SE.factor = 2)
# Global GapS optimiz criterias point tow k=9 or 10
print(cg.sart,method = 'firstSEmax', SE.factor = 2)
print(cg.sart,method = 'firstmax')
# First max settle with the peak in the graph at k=3
print(cg.sart,method = 'Tibs2001SEmax', SE.factor = 2)
# This chooses k = 1 since all stdev are larger than the differences f(k+1)-f(k).
# Best sol appears t be 3 since we know the DGP of such dataset - local maximum

###

set.seed(1234)
U1=runif(140, max = max(artdata$V1),min = min(artdata$V1))
U2=runif(140, max = max(artdata$V2),min = min(artdata$V2))
U=data.frame(U1,U2)
plot(U)

km.out=list()
set.seed(1234)
for (i in 1:10) {
  km.out[[i]]=kmeans(U,iter.max = 100,nstart =100, centers = i)
}

plot(U,col=(km.out[[2]]$cluster+2),cex=1.2,pch=0,lwd=2,
     main='Artificial Data Clustering,  K=2')
plot(U,col=(km.out[[4]]$cluster),cex=1.4,pch=0,lwd=2,
     main='Artificial Data Clustering,  K=4')
plot(U,col=(km.out[[6]]$cluster+2),cex=1.2,pch=0,lwd=2,
     main='Artificial Data Clustering,  K=6')

U.logSk=rep(NA,10)
for (i in 1:10) {
U.logSk[i]=log(km.out[[i]]$tot.withinss)  
}

set.seed(1234)
cg.sart=clusGap(artdata,kmeans, K.max=10, B=300, d.power=2,
                spaceH0 = 'scaledPCA', nstart=100, iter.max=100)
AD2.logSk=cg.sart$Tab[,1]

plot(1:10,U.logSk,type='l',col='grey',lty=2,lwd=2,ylim=c(min(AD2.logSk),11),
     ylab='log S_k')
lines(1:10,AD2.logSk,lwd=2)
legend(7,11,legend=c('log(Sk) dataset','log(SK) Uniform'), lty=1:2,lwd=3,col=c(1,'grey'))




## EXERCISE 2 ####
# GENERATE FURTHER DATA SETS as ArtificialData1
data=list()
simruns=100
kmax=10
library(sn)
set.seed(1234)
for (i in 1:simruns) {
  v1 <- c(rnorm(50,0,1), rsn(70,5,1,8), rnorm(30,6,1))
  v2 <- c(rnorm(50,0,1), rsn(70,0,1,8), 8+rt(30,5))
  data[[i]] <- cbind(v1,v2)
}

require(cluster)

gapnc <- function(data,FUNcluster=kmeans,
                  K.max=10, B = 100, d.power = 2,
                  spaceH0 ="scaledPCA",
                  method ="globalSEmax", SE.factor = 2,...){
  # As in original clusGap function the ... arguments are passed on
  # to the clustering method FUNcluster (kmeans).
  # Run clusGap
  gap1 <- clusGap(data,kmeans,K.max, B, d.power,spaceH0,...)
  # Find optimal number of clusters; note that the method for
  # finding the optimum and the SE.factor q need to be specified here.
  nc <- maxSE(gap1$Tab[,3],gap1$Tab[,4],method, SE.factor)
  # Re-run kmeans with optimal nc.
  kmopt <- kmeans(data,nc,...)
  out <- list()
  out$gapout <- gap1
  out$nc <- nc
  out$kmopt <- kmopt
  out
}

library(parallel)

out.orig1=list()
set.seed(1234)
for (i in 1:100) {
  out.orig1[[i]]=gapnc(data[[i]], K.max = kmax ,spaceH0='original',SE.factor=1,
                       iter.max=100, nstart=2)}
clu.orig1=rep(NA,100)
for (i in 1:100) {clu.orig1[i]=out.orig1[[i]]$nc}


#####
# ORIGINAL SPACE, SE.FACTOR=2
out.orig2=list()
set.seed(1234)
for (i in 1:100) {
  out.orig2[[i]]=gapnc(data[[i]],K.max = kmax,spaceH0='original',SE.factor=2,
                       iter.max=100, nstart=20)}

out.pca1=list()
set.seed(1234)
for (i in 1:100) {
  out.pca1[[i]]=gapnc(data[[i]],K.max = kmax, spaceH0='scaledPCA',SE.factor=1,
                      iter.max=100, nstart=20)}

out.pca2=list()
set.seed(1234)
for (i in 1:100) {
  out.pca2[[i]]=gapnc(data[[i]],K.max = kmax,spaceH0='scaledPCA',SE.factor=2,
                      iter.max=100, nstart=20)}
#####

set.seed(77665544)
x1 <- rnorm(20)
y1 <- rnorm(20)
x2 <- rnorm(20,mean=10)
y2 <- rnorm(20)
x3 <- runif(100,-20,30)
y3 <- runif(100,20,40)
clusterdata2 <- cbind(c(x1,x2,x3),c(y1,y2,y3))
plot(clusterdata2)












