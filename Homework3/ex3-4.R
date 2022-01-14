### EXERCISES 3-4 - MS&BDA 
### ZURLO GIOVANNI 19/10/2021

##  EXERCISE 1 #### 
## a) ####
library(cluster)
x1 = c("blue" ,  T, T,  F,  12) #4th is F
x2 = c("red"  ,  F, F, NA,  NA)
x3 = c("red"  ,  T, F, NA,  17)
x4 = c("green",  T, F,  F,  21)

Gow.d12=(1+1+1)/3
Gow.d13=(1+0+1+(17-12)/9)/4
Gow.d14=(1+0+1+0+1)/4
Gow.d23=(0+1+0)/2
Gow.d24=(1+1)/2
Gow.d34=(1+0+0+(21-17)/9)/3

dist=c(Gow.d12,Gow.d13,Gow.d23,Gow.d14,Gow.d24,Gow.d34)
dist.mat <- matrix(0, nrow = 3, ncol = 3)
rownames(dist.mat)=c("x1","x2","x3")
colnames(dist.mat)=c("x2","x3","x4")
dist.mat[upper.tri(dist.mat, diag = TRUE)] <- dist
t(dist.mat)

## b) ####
j.d12=(1-0)
j.d13=(1-1/2)
j.d14=(1-1/2)
j.d23=(1)
j.d24=(1)
j.d34=(1-1)

bGow.d12=(1+3*j.d12)/4
bGow.d13=(1+3*j.d13+(17-12)/9)/5
bGow.d14=(1+3*j.d14+1)/5
bGow.d23=(0+3*j.d23)/4
bGow.d24=(1+3*j.d24)/4
bGow.d34=(1+3*j.d34+(21-17)/9)/5

bdist=c(bGow.d12,bGow.d13,bGow.d23, bGow.d14, bGow.d24, bGow.d34)
bdist.mat <- matrix(0, nrow = 3, ncol = 3)
rownames(bdist.mat)=c("x1","x2","x3")
colnames(bdist.mat)=c("x2","x3","x4")
bdist.mat[upper.tri(bdist.mat, diag = TRUE)] <- bdist 
t(bdist.mat)


## c) ####
x=as.data.frame(rbind(x1,x2,x3,x4))
is.na(x[2,4:5])
str(x)
x[,1]=as.factor(x[,1])
x[,2]=as.logical(x[,2])
x[,3]=as.logical(x[,3])
x[,4]=as.logical(x[,4])
x[,5]=as.numeric(x[,5])
str(x)
daisy(x,'gower')

## EXERCISE 2 ####
library(MASS)
set.seed(123)
corr.d <- function(a,b) {0.5*(1-cor(a,b))}

#(x=rnorm(3))
#(y=round(rnorm(1,6,1))*x)
#cor(x,y)
#(z=-x)
#cor(x,z)
#cor(y,z)

set.seed(1234)
S<- matrix(c(3,0.5,1,0.5), nrow=2) # covariance matrix of X and Z
rmatrix <- mvrnorm(n=4, mu=c(0,0),Sigma=S, empirical=TRUE)
X <- rmatrix[,1] # mean 0, variance 3
Z <- rmatrix[,2] # mean 0, variance 1
cor(X,Z)
Y  <- X + Z
corr.d(X,Z)
corr.d(X,Y) 
corr.d(Y,Z)

# PROOF
cor(X,Y)+cor(Y,Z)<=1+cor(X,Z)
# OR
corr.d(X,Z)<=corr.d(X,Y) + corr.d(Y,Z)

## GOWER TRIANGLE ####
X = c("blue"  ,  T,NA, T,  1,"1")
Y = c("blue"  ,  NA,T, F,  4,"4")
Z = c("green" ,  F,T, F, 1,"2")

x=as.data.frame(rbind(X,Y,Z))
x[,1]=as.factor(x[,1])
x[,2]=as.logical(x[,2])
x[,3]=as.logical(x[,3])
x[,4]=as.logical(x[,4])
x[,5]=as.numeric(x[,5])
x[,6]=as.ordered(x[,6])
str(x)
library(cluster)
d=daisy(x,'gower', type = list(symm=2))
(diss=as.matrix(d))

diss[3,1]<=diss[2,1] + diss[3,2]


## EXERCISE 4####
library(factoextra)
library(cluster)
library(prabclus)
library(smacof)
library(ggdendro)
library(mclust)

covid2021 <- read.table("covid2021.dat")
x <- covid2021[,5:559] # This selects the variables for clustering
dim(x)

plot(1:555,x[1,],type="l",ylim=c(0,25),
     ylab="New cases over one week per 1000 inhabitants",
     xlab="Day (1 April 2020-7 October 2021)")
for (i in 2:179) {points(1:555,x[i,],type="l")}


## Hierarchical Ward Clustering + Euclidean ####
d.eucl=dist(x)
out.eucl=hclust(d.eucl,method = "ward.D2")
plot(out.eucl, hang=-1, xlab="",sub="" ,cex=0.6, cex.axis=1)

b=ggdendrogram(out.eucl, rotate = F, theme_dendro = F,labels = F)
print(b + ggtitle("Ward"))


par(mar=c(4,4,2,2))
# store the dedrogram in an object
dhc=as.dendrogram(out.eucl)
plot(dhc[[2]] , main= "zoom on a part of the dendrogram") #2nd group

tail(cbind(out.eucl$merge,out.eucl$height),13) # Last aggregations

plot(3:20,rev(out.eucl$height)[3:20],type="b",
     xlab="K",ylab="Height",
     cex.lab=0.8,cex.axis=0.7,
     main="Scree Plot", cex.main=0.9)
grid()
# Increse in the agglomerative distance at k=12, h=58.35897 OR 8!

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:20){
  tclusk[[k]] <- cutree(out.eucl,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=d.eucl)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=d.eucl))$avg.width
}
plot(1:20,tasw,type="b",xlab="Number of clusters",ylab="ASW")
grid() # 8 is suggested ASW=0.355, 3 clusters the best (western + all the rest)

plot(tsil[[8]], col=terrain.colors(8),
     main="Silhouette plot - Ward eucl - 8 cluster")

groups.eucl=cutree(out.eucl, k=3)
rect.hclust(out.eucl, k=3, border="red")
table(groups.eucl) # Many groups with very reduced size (outlier filtering)
# euclidean tend to isolate outliers

library(mclust)
adjustedRandIndex(covid2021$continent,groups.eucl)


## Hierarchical Ward Clustering + Correlation ####
r.mat=cor(t(x))
d.corr=as.dist(0.5*(1-r.mat)) # Why this version?
out.corr=hclust(d.corr,method = "ward.D2")
plot(out.corr, hang=-1, xlab="",sub="" ,cex=0.6, cex.axis=1)
rect.hclust(out.corr,k=7)

par(mar=c(6,6,6,6))
# store the dedrogram in an object
dhc=as.dendrogram(out.corr)
plot(dhc[[1]] , main= "Zoom on 1st Cluster - Western Countries",
     xlab="",sub="" ,cex=0.6, cex.axis=1)

tail(cbind(out.corr$merge,out.corr$height),15) # Last aggregations

plot(2:20,rev(out.corr$height)[2:20],type="b",
     xlab="K",ylab="Height",
     cex.lab=0.8,cex.axis=0.7,
     main="Scree Plot", cex.main=0.9)
grid()
# Increse in the agglomerative distance at k=7, h=1.0136557 OR 16(too much)!
# 14 but 7 may be also a good compromise

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:20){
  tclusk[[k]] <- cutree(out.corr,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=d.corr)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=d.corr))$avg.width
}
plot(1:20,tasw,type="b",xlab="Number of clusters",ylab="ASW")
grid() # 5 (corresp. with continent) is suggested with ASW=0.24, 14!

plot(tsil[[4]], col=terrain.colors(4),
     main="Silhouette plot - Ward corr - 5 cluster")

groups.corr=cutree(out.corr, k=5)
plot(out.corr, hang=-1, xlab="",sub="" ,cex=0.6, cex.axis=1)
rect.hclust(out.corr, k=5, border="red")
table(groups.corr) # Minimum to avoid small clusters.. good homogeneity

library(mclust)
adjustedRandIndex(groups.corr,groups.eucl) # Very very different results

## REPRESENTATION ####
mds.eucl=mds(d.eucl, type = 'ratio')
mds.corr=mds(d.corr, type = 'ratio')
c(mds.corr$stress,mds.eucl$stress) 
# 16%, not much information missing for such dataset - 0.265 a bit high

plot(mds.eucl$conf, main="MDS Plot for Euclidean DistMat - Ward",
     col=(groups.eucl+6), cex=1.2,pch=0,lwd=2, xlim=c(-2,0.6),ylim=c(-1.3,2.5),
     xlab="Stress  =  16.5%",ylab='')

plot(mds.corr$conf, asp=1,main="MDS Plot for Correlation DissMat - Ward",
     col=(groups.eucl+6), cex=1.2,pch=0,lwd=2,
     xlab="Stress  =  26.5%",ylab='')


## Hierarchical AVG Clustering + Euclidean ####
d.eucl=dist(x)
out.euclAVG=hclust(d.eucl,method = "average")
plot(out.euclAVG, hang=-1, xlab="",sub="" ,cex=0.6, cex.axis=1)

b=ggdendrogram(out.eucl, rotate = F, theme_dendro = F,labels = F)
print(b + ggtitle("Ward"))


par(mar=c(4,4,2,2))
# store the dedrogram in an object
dhc=as.dendrogram(out.euclAVG)
plot(dhc[[2]] , main= "zoom on a part of the dendrogram") #2nd group

tail(cbind(out.euclAVG$merge,out.euclAVG$height),13) # Last aggregations

plot(2:20,rev(out.euclAVG$height)[2:20],type="b",
     xlab="K",ylab="Height",
     cex.lab=0.8,cex.axis=0.7,
     main="Scree Plot", cex.main=0.9)
grid()
# CLEAR ELBOW in the agglomerative distance at k=9, h=41.33749 

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:20){
  tclusk[[k]] <- cutree(out.euclAVG,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=d.eucl)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=d.eucl))$avg.width
}
plot(1:20,tasw,type="b",xlab="Number of clusters",ylab="ASW")
grid() # not clear... cliff at k=9 with ASW=0.46

plot(tsil[[9]], col=terrain.colors(256),
     main="Silhouette plot - euclAVG - 9 cluster")
summary(tsil[[9]]) # Several individ. clusters, 9 high

groups.euclAVG=cutree(out.euclAVG, k=9)
rect.hclust(out.euclAVG, k=9, border="red")
table(groups.euclAVG) # It basically leaves out several "outliers"; don't like

adjustedRandIndex(groups.eucl,groups.euclAVG) #0.125 very low


## Hierarchical AVG Clustering + Correlation ####
r.mat=cor(t(x))
d.corr=as.dist(0.5*(1-r.mat))
out.corrAVG=hclust(d.corr,method = "average")
plot(out.corrAVG, hang=-1, xlab="",sub="" ,cex=0.6, cex.axis=1)
rect.hclust(out.corrAVG,k=6)

par(mar=c(5,5,2,2))
# store the dedrogram in an object
dhc=as.dendrogram(out.corr)
plot(dhc[[1]] , main= "zoom on a part of the dendrogram")

tail(cbind(out.corrAVG$merge,out.corrAVG$height),15) # Last aggregations

plot(3:20,rev(out.corrAVG$height)[3:20],type="b",
     xlab="K",ylab="Height",
     cex.lab=0.8,cex.axis=0.7,
     main="Scree Plot", cex.main=0.9)
grid()
# No elbow... maybe 6 or 13

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:20){
  tclusk[[k]] <- cutree(out.corrAVG,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=d.corr)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=d.corr))$avg.width
}
plot(1:20,tasw,type="b",xlab="Number of clusters",ylab="ASW")
grid() # 4 or 6, about ASW= 0.26

plot(tsil[[4]], col=terrain.colors(4),
     main="Silhouette plot - avg corr - 4 cluster") #more homogeneous clust

groups.corrAVG=cutree(out.corrAVG, k=6)
rect.hclust(out.corr, k=6, border="red")
table(groups.corrAVG) 

adjustedRandIndex(groups.corr,groups.corrAVG) #0.29 little bit more similar
