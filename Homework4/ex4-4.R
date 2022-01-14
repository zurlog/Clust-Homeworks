### EXERCISES 4-4 - MS&BDA 
### ZURLO GIOVANNI 07/11/2021

##  EXERCISE 1 #### 
library(factoextra)
library(cluster)
library(prabclus)
library(smacof)

data(tetragonula)
# In this form, the alleles are coded by numbers (every locus has a six digit number)

ta <- alleleconvert(strmatrix=tetragonula)
# converts the allele codes to letters, so that each locus now has two letters

tai <- alleleinit(allelematrix=ta)
# produces a collection of ways to represent the data

dist=as.dist(tai$distmat)

## a)
mds=mds(dist, type = "ratio")
plot(mds$conf, main="MDS Plot for Tetragonula Genes DistMat",
     cex=1.2,pch=1,lwd=2.2,
      xlim=c(-1,0.8),ylim=c(-0.8,1),
     xlab="Stress  =  26.2%",ylab='')
# much information missing for such dataset - 0.262 a bit high

library(ggplot2)
library(ggdendro)
library(ggpubr)
library(dplyr)

## COMPLETE ####
complete <- hclust(dist,method="complete")
a=ggdendrogram(complete, rotate = F, theme_dendro = F,labels = F)
print(a + ggtitle("Complete Linkage") + xlab("") + ylab("Height"))

plot(complete, hang=-1, xlab="Dendrogram for Complete Linkage AAHC",sub="" ,cex=0.6, cex.axis=1)

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:30){
  tclusk[[k]] <- cutree(complete,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=dist)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=dist))$avg.width
}

plot(1:30,tasw,type="b",ylab="ASW",
     xlab=paste("Number of clusters   -   Max ASW =",as.character(round(max(tasw[-1]),3))))
k.complete=match(max(tasw[-1]),tasw[-1])+1
abline(v=k.complete,col="red",lty="dashed")

complete.out = cutree(complete,k.complete)
complete.f=as.factor(complete.out)

mds.df=as.data.frame(mds$conf)
mds.df= cbind(mds.df,complete.f)

## VISUALIZATION ####
library(fpc)
library()

plot(mds$conf, main="MDS Plot for Tetragonula Genes DistMat",
     cex=1.2,pch=clusym[complete.f],lwd=2.2,
     xlim=c(-1,0.8),ylim=c(-0.8,1),
     xlab="Stress  =  26.2%",ylab='',
     col=complete.f)

#####
## AVERAGE ####
average <- hclust(dist,method="average")
a=ggdendrogram(average, rotate = F, theme_dendro = F,labels = F)
print(a + ggtitle("Average Linkage") + xlab("") + ylab("Height"))

plot(average, hang=-1, xlab="Dendrogram for Average Linkage AAHC",sub="" ,cex=0.6, cex.axis=1)

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:30){
  tclusk[[k]] <- cutree(average,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=dist)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=dist))$avg.width
}
plot(1:30,tasw,type="b",ylab="ASW",
     xlab=paste("Number of clusters   -   Max ASW =",as.character(round(max(tasw[-1]),3))))
k.average=match(max(tasw[-1]),tasw[-1])+1
abline(v=k.average,col="red",lty="dashed")

average.out = cutree(average,k.average)
average.f=as.factor(average.out)

mds.df=as.data.frame(mds$conf)
mds.df= cbind(mds.df,average.f)

## VISUALIZATION ####
library(fpc)

plot(mds$conf, main="MDS Plot for Tetragonula Genes DistMat",
     cex=1.2,pch=clusym[average.f],lwd=2.2,
     xlim=c(-1,0.8),ylim=c(-0.8,1),
     xlab="Stress  =  26.2%",ylab='',
     col=average.f)


gg1=ggscatter(mds.df,
               x="D1", y="D2",
               color = "average.f", palette = "simpsons",
               shape = "average.f",
               ellipse = F, 
               mean.point = TRUE,
               star.plot = TRUE)
gg1+theme(legend.position = "none")





#####
## PAM ####
pasw <- NA
pclusk <- list()
psil <- list()
for (k in 2:30){
  pclusk[[k]] <- pam(dist,k)
  # Computation of silhouettes:
  psil[[k]] <- silhouette(pclusk[[k]],dist=dist)
  # ASW needs to be extracted:
  pasw[k] <- summary(psil[[k]])$avg.width
}


plot(1:30,pasw,type="b",ylab="ASW",
     xlab=paste("Number of clusters   -   Max ASW =",as.character(round(max(pasw[-1]),3))))
m=match(max(pasw[-1]),pasw[-1])+1
abline(v=m,col="red",lty="dashed")

pam.f=as.factor(pclusk[[m]]$cluster)
mds.df=as.data.frame(mds$conf)
mds.df= cbind(mds.df,pam.f)

## VISUALIZATION ####

plot(mds$conf, main="MDS Plot for Tetragonula Genes DistMat",
     cex=1.2,pch=clusym[pam.f],lwd=2.2,
     xlim=c(-1,0.8),ylim=c(-0.8,1),
     xlab="Stress  =  26.2%",ylab='',
     col=pam.f)

gg1=ggscatter(mds.df,
               x="D1", y="D2",
               color = "pam.f", palette = "simpsons",
               shape = "pam.f",
               ellipse = F, 
               mean.point = TRUE,
               star.plot = TRUE)
gg1+theme(legend.position = "none")

######

summary<- matrix(0, nrow = 1, ncol = 3)

colnames(summary)=c("Complete","Average","Pam")
rownames(summary)=c("ASW")
summary[1,]=c(0.472,0.486,0.468)
summary





## POINT B) ####
mds=mds(dist, type = "ratio")
head(mds$conf)
set.seed(1234)
cg1 <- clusGap(mds$conf,kmeans,12,B=200,d.power=2,spaceH0="scaledPCA",nstart=200)
plot(cg1,main="")
print(cg1,method="globalSEmax",SE.factor=2)
print(cg1,method="Tibs2001SEmax",SE.factor=1)
print(cg1,method="firstSEmax",SE.factor=2)

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:20){
  tclusk[[k]] <- kmeans(mds$conf,centers=k,iter.max = 100,nstart = 100)
  tsil[[k]] <- silhouette(tclusk[[k]]$cluster,dist=dist(mds$conf))
  tasw[k] <- summary(tsil[[k]])$avg.width
}
plot(1:20,tasw,type="b",ylab="ASW",
     xlab=paste("Number of clusters   -   Max ASW =",as.character(round(max(tasw[-1]),3))))
grid()
m=match(max(tasw[-1]),tasw[-1])+1
abline(v=m,col="red",lty="dashed")

# Choice 7
set.seed(1234)
kmeans.out=kmeans(mds$conf,7,iter.max = 300,nstart = 300)
kmeans.f=as.factor(kmeans.out$cluster)

mds.df=as.data.frame(mds$conf)
mds.df= cbind(mds.df,kmeans.f)

gg1=ggscatter(mds.df,
              x="D1", y="D2",
              color = "kmeans.f", palette = "simpsons",
              shape = "kmeans.f",
              ellipse = F, 
              mean.point = TRUE,
              star.plot = TRUE)
gg1+theme(legend.position = "none")+ggtitle("MDS Plot for the Genetic DistMat - KMeans (7)")

# All fantastic? not really
# Partition based on a dataset with loss of information - about 26%
# Looks good on paper, let's try it on actual data

library(mclust)
adjustedRandIndex(kmeans.f,average.f)
# Somehow Similar (94%) but...

real_ksil<- silhouette(kmeans.out$cluster,as.dist(tai$distmat))
real_avgwidth <- summary(real_ksil)$avg.width
plot(real_ksil)
# IT DROPS SIGNIFICANTLY - Major issues in clusters 2-4

## Try mds on 3dimensions...####
library(smacof)
mds3=mds(dist,ndim=3,type = "ratio")

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:20){
  tclusk[[k]] <- kmeans(mds3$conf,centers=k,iter.max = 100,nstart = 100)
  tsil[[k]] <- silhouette(tclusk[[k]]$cluster,dist=dist(mds3$conf))
  tasw[k] <- summary(tsil[[k]])$avg.width
}
plot(1:20,tasw,type="b",ylab="ASW",
     xlab=paste("Number of clusters   -   Max ASW =",as.character(round(max(tasw[-1]),3))))
grid()
m=match(max(tasw[-1]),tasw[-1])+1
abline(v=m,col="red",lty="dashed")

set.seed(1234)
kmeans11.out=kmeans(mds3$conf,11,iter.max = 300,nstart = 300)
kmeans11.f=as.factor(kmeans11.out$cluster)

mds.df=as.data.frame(mds3$conf)
mds.df= cbind(mds.df,kmeans11.f)

adjustedRandIndex(kmeans11.f,average.f)
# NOW ONLY 87%

real_k11sil<- silhouette(kmeans11.out$cluster,as.dist(tai$distmat))
real_avg11width <- summary(real_k11sil)$avg.width
# Small ASW improvement up to 0.44
plot(real_k11sil)

######


## EXERCISE 2 #####
library(pdfCluster)
library(cluster)
library(mclust)
data("oliveoil")
olivex=oliveoil[,3:10]
solive=scale(olivex)

## 1) ####
set.seed(1234)
library(parallel)
# Computing the Gap Statistic for each k-means solution up to k = 14
cg1 <- clusGap(solive,kmeans,15,B=100,d.power=2,spaceH0="scaledPCA",nstart=100)
plot(cg1,main="ClusGap Plot for K-means",xlab="Number of Clusters",ylab="Gap Statistic")
abline(v=18,col="red",lty="dashed")

print(cg1,method="globalSEmax",SE.factor=2)
print(cg1,method="Tibs2001SEmax",SE.factor=2)
print(cg1,method="firstSEmax",SE.factor=2)

kmeans18=kmeans(solive,18,iter.max = 2000,nstart = 2000)
adj.rand.index(kmeans18$cluster,oliveoil$macro.area) #0.19
adj.rand.index(kmeans18$cluster,oliveoil$region) #0.424

library(kableExtra)
rownomi=c("Macro Area","Region")

a=cbind(0.19,0.424) 
colnames(a)=rownomi  
rownames(a)="K-Means (18)"
a |>
  kable("html", caption = 'Clustering Similarity measured by ARI', row.names = T) |>
  kable_styling(full_width = F, position = "center")

## 2) ####

library(factoextra)
set.seed(1234)
fviz_nbclust(solive, FUNcluster = hcut, k.max = 20, hc_method = "ward.D2",
 method = "gap_stat", maxSE = list(method = "globalSEmax", SE.factor = 2))

ward16=hclust(dist(solive),method = "ward.D2") |>
cutree(k=16)
adj.rand.index(ward16,oliveoil$macro.area) #0.209
adj.rand.index(ward16,oliveoil$region) #0.472

a=cbind(0.209,0.472) 
colnames(a)=rownomi  
rownames(a)="Ward.D2 (16)"
a |>
  kable("html", caption = 'Clustering Similarity measured by ARI', row.names = T) |>
  kable_styling(full_width = F, position = "center")

## 3) ####

fviz_nbclust(solive, FUNcluster = hcut, k.max = 20, hc_method = "single",
             hc_metric="euclidean", method = "silhouette")
# prefers the two clust solution but this is a bias - look local max 5
# Similarity is approx. zero when selecting k=2

single.eucl5=hclust(dist(solive),method = "single") |>
  cutree(k=5)
adj.rand.index(single.eucl5,oliveoil$macro.area) #0.81
adj.rand.index(single.eucl5,oliveoil$region) #0.356

fviz_nbclust(solive, FUNcluster = hcut, k.max = 20, hc_method = "single",
             hc_metric="manhattan", method = "silhouette")
# look local max 7

single.manh7=hclust(dist(solive,'manhattan'),method = "single") |>
  cutree(k=7)
adj.rand.index(single.manh7,oliveoil$macro.area) #0.794
adj.rand.index(single.manh7,oliveoil$region) #0.346

## MAHALANOBIS DOESN'T WORK ####

mahaladist <- matrix(0,ncol=572,nrow=572)
olivecov <- cov(olivex)
for (i in 1:572) {
  mahaladist[i,] <- mahalanobis(olivex,as.numeric(olivex[i,]),olivecov)
}
mahaladist=as.dist(mahaladist)
mahalanobis <- hclust(mahaladist,method="single")

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:20){
  tclusk[[k]] <- cutree(mahalanobis,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=mahaladist)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=mahaladist))$avg.width
}

plot(1:20,tasw,type="b",ylab="ASW",
     xlab=paste("Number of clusters   -   Max ASW =",as.character(round(max(tasw[-1]),3))))
k.mahala=match(max(tasw[-1]),tasw[-1])+1
abline(v=k.mahala,col="red",lty="dashed")

mahala.out = cutree(mahalanobis,5)
mahala.f=as.factor(mahala.out)

adj.rand.index(mahala.f,oliveoil$macro.area) #0.008
adj.rand.index(mahala.f,oliveoil$region) #0.03


a=rbind(0.81,0.356)
b=rbind(0.794,0.346)
c=rbind(0.008,0.003)
table=cbind(a,b,c)
colnomi=c('Euclidean','Manhattan','Mahalanobis')
colnames(table)=colnomi  
rownames(table)=rownomi
table |>
  kable("html", caption = 'ARI - Single Linkage AAHC', row.names = T) |>
  kable_styling(full_width = F, position = "center")

## 4) ####

fviz_nbclust(solive, FUNcluster = hcut, k.max = 20, hc_method = "average",
             hc_metric="euclidean", method = "silhouette")
# global max 10
# Similarity is approx. zero when selecting k=2

avg.eucl10=hclust(dist(solive),method = "average") |>
  cutree(k=10)
adj.rand.index(avg.eucl10,oliveoil$macro.area) #0.559
adj.rand.index(avg.eucl10,oliveoil$region) #0.827

fviz_nbclust(solive, FUNcluster = hcut, k.max = 20, hc_method = "average",
             hc_metric="manhattan", method = "silhouette")
# global max 8

avg.manh8=hclust(dist(solive,'manhattan'),method = "average") |>
  cutree(k=8)
adj.rand.index(avg.manh8,oliveoil$macro.area) #0.481
adj.rand.index(avg.manh8,oliveoil$region) #0.777

## MAHALANOBIS DOESN'T WORK ####

mahaladist <- matrix(0,ncol=572,nrow=572)
olivecov <- cov(olivex)
for (i in 1:572) {
  mahaladist[i,] <- mahalanobis(olivex,as.numeric(olivex[i,]),olivecov)
}
mahaladist=as.dist(mahaladist)
mahalanobis <- hclust(mahaladist,method="average")

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:20){
  tclusk[[k]] <- cutree(mahalanobis,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=mahaladist)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=mahaladist))$avg.width
}

plot(1:20,tasw,type="b",ylab="ASW",
     xlab=paste("Number of clusters   -   Max ASW =",as.character(round(max(tasw[-1]),3))))
k.mahala=match(max(tasw[-1]),tasw[-1])+1
abline(v=k.mahala,col="red",lty="dashed")
# local max is 12

mahala.out = cutree(mahalanobis,12)
mahala.f=as.factor(mahala.out)

adj.rand.index(mahala.f,oliveoil$macro.area) #0.124
adj.rand.index(mahala.f,oliveoil$region) #0.171


a=rbind(0.559,0.827)
b=rbind(0.481,0.777)
c=rbind(0.124,0.171)
table=cbind(a,b,c)
colnomi=c('Euclidean','Manhattan','Mahalanobis')
colnames(table)=colnomi  
rownames(table)=rownomi
table |>
  kable("html", caption = 'ARI - Average Linkage AAHC', row.names = T) |>
  kable_styling(full_width = F, position = "center")


## 5) ####

fviz_nbclust(solive, FUNcluster = hcut, k.max = 20, hc_method = "complete",
             hc_metric="euclidean", method = "silhouette")
# k=2 bias... try 10
# Similarity is approx. zero when selecting k=2

compl.eucl10=hclust(dist(solive),method = "complete") |>
  cutree(k=10)
adj.rand.index(compl.eucl10,oliveoil$macro.area) #0.431
adj.rand.index(compl.eucl10,oliveoil$region) #0.670

fviz_nbclust(solive, FUNcluster = hcut, k.max = 20, hc_method = "complete",
             hc_metric="manhattan", method = "silhouette")
# global max 7

compl.manh7=hclust(dist(solive,'manhattan'),method = "complete") |>
  cutree(k=7)
adj.rand.index(compl.manh7,oliveoil$macro.area) #0.458
adj.rand.index(compl.manh7,oliveoil$region) #0.673

## MAHALANOBIS DOESN'T WORK ####

mahaladist <- matrix(0,ncol=572,nrow=572)
olivecov <- cov(olivex)
for (i in 1:572) {
  mahaladist[i,] <- mahalanobis(olivex,as.numeric(olivex[i,]),olivecov)
}
mahaladist=as.dist(mahaladist)
mahalanobis <- hclust(mahaladist,method="complete")

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:20){
  tclusk[[k]] <- cutree(mahalanobis,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=mahaladist)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=mahaladist))$avg.width
}

plot(1:20,tasw,type="b",ylab="ASW",
     xlab=paste("Number of clusters   -   Max ASW =",as.character(round(max(tasw[-1]),3))))
k.mahala=match(max(tasw[-1]),tasw[-1])+1
abline(v=k.mahala,col="red",lty="dashed")
# local max is 10 or 20

mahala.out = cutree(mahalanobis,20)
mahala.f=as.factor(mahala.out)

adj.rand.index(mahala.f,oliveoil$macro.area) #0.274
adj.rand.index(mahala.f,oliveoil$region) #0.520


a=rbind(0.431,0.670)
b=rbind(0.458,0.673)
c=rbind(0.274,0.520)
table=cbind(a,b,c)
colnomi=c('Euclidean','Manhattan','Mahalanobis')
colnames(table)=colnomi  
rownames(table)=rownomi
table |>
  kable("html", caption = 'ARI - Complete Linkage AAHC', row.names = T) |>
  kable_styling(full_width = F, position = "center")

## 6) ####

pasw <- NA
pclusk <- list()
psil <- list()
set.seed(1234)
for (k in 2:20){
  pclusk[[k]] <- pam(solive,k, metric = 'euclidean', keep.diss = T)
  # Computation of silhouettes:
  psil[[k]] <- silhouette(pclusk[[k]])
  # ASW needs to be extracted:
  pasw[k] <- summary(psil[[k]])$avg.width
}


plot(1:20,pasw,type="b",ylab="ASW",
     xlab=paste("Number of clusters   -   Max ASW =",as.character(round(max(pasw[-1]),3))))
m=match(max(pasw[-1]),pasw[-1])+1
abline(v=m,col="red",lty="dashed")
# Global max at *K = 5*

pam.eucl5=as.factor(pclusk[[m]]$cluster)
adj.rand.index(pam.eucl5,oliveoil$macro.area) #0.559
adj.rand.index(pam.eucl5,oliveoil$region) #0.755


pasw <- NA
pclusk <- list()
psil <- list()
for (k in 2:20){
  pclusk[[k]] <- pam(solive,k, metric = 'manhattan', keep.diss = T)
  # Computation of silhouettes:
  psil[[k]] <- silhouette(pclusk[[k]])
  # ASW needs to be extracted:
  pasw[k] <- summary(psil[[k]])$avg.width
}


plot(1:20,pasw,type="b",ylab="ASW",
     xlab=paste("Number of clusters   -   Max ASW =",as.character(round(max(pasw[-1]),3))))
m=match(max(pasw[-1]),pasw[-1])+1
abline(v=m,col="red",lty="dashed")
# Global max at *K = 5*

pam.manh5=as.factor(pclusk[[m]]$cluster)
adj.rand.index(pam.manh5,oliveoil$macro.area) #0.553
adj.rand.index(pam.manh5,oliveoil$region) #0.757


## MAHALANOBIS DOESN'T WORK ####

mahaladist <- matrix(0,ncol=572,nrow=572)
olivecov <- cov(olivex)
for (i in 1:572) {
  mahaladist[i,] <- mahalanobis(olivex,as.numeric(olivex[i,]),olivecov)
}
mahaladist=as.dist(mahaladist)


pasw <- NA
pclusk <- list()
psil <- list()
for (k in 2:20){
  pclusk[[k]] <- pam(mahaladist, k, keep.dis = T)
  # Computation of silhouettes:
  psil[[k]] <- silhouette(pclusk[[k]])
  # ASW needs to be extracted:
  pasw[k] <- summary(psil[[k]])$avg.width
}


plot(1:20,pasw,type="b",ylab="ASW",
     xlab=paste("Number of clusters   -   Max ASW =",as.character(round(max(pasw[-1]),3))))
m=match(max(pasw[-1]),pasw[-1])+1
abline(v=m,col="red",lty="dashed")
# Global max at *K = 6*

pam.mahala6=as.factor(pclusk[[m]]$cluster)
adj.rand.index(pam.mahala6,oliveoil$macro.area) #0.597
adj.rand.index(pam.mahala6,oliveoil$region) #0.712

a=rbind(0.559,0.755)
b=rbind(0.553,0.757)
c=rbind(0.597,0.712)
table=cbind(a,b,c)
colnomi=c('Euclidean (5)','Manhattan (5)','Mahalanobis (6)')
colnames(table)=colnomi
rownomi=c("Macro Area","Region")
rownames(table)=rownomi
table |>
  kable("html", caption = 'ARI - PAM Algorithm', row.names = T) |>
  kable_styling(full_width = F, position = "center")


## DBSCAN #### 

setwd("C:/Users/giova/Dropbox/Magistrale/[NEW] MS & BD Analytics/Datasets")
bundestag=read.table("bundestag.dat",header = T,sep="")
flexclust::bundestag(2005)
bundestag$ewb=as.factor(bundestag$ewb)
bundestag$state=as.factor(bundestag$state)
data=bundestag[,1:5]
library(dbscan)
library(cluster)
library(mclust)
library(fpc)

kNNdistplot(data, k=5)
abline(h=0.06, col=2, lty=2)
fpc::dbscan() 
print(db)
print.dbscan()

#One limitation of DBSCAN is that it is sensitive to the choice of epsi, global par
# density of the thinnest cluster, that is lowest density which is not considered noise
#in particular if clusters have different densities. 
#If epsi is too small, sparser clusters will be defined as noise.
#If epsi is too large, denser clusters may be merged together.

out.dbs5=dbscan(data, eps = 0.05, minPts = 5)
out.dbs6=dbscan(data, eps = 0.06, minPts = 5)
out.dbs7=dbscan(data, eps = 0.065, minPts = 5)
out.dbs8=dbscan(data, eps = 0.07, minPts = 5)
adjustedRandIndex(out.dbs5$cluster,out.dbs8$cluster)

print(out.dbs5)
out.dbs=dbscan::dbscan(data, eps = 0.068, minPts = 6)
out.dbs2=fpc::dbscan(data, eps = 0.068, MinPts = 6)
print(out.dbs2)
col2=c("gold2","mediumpurple","tomato2")
pairs(data,cex=0.8,col=(out.dbs$cluster+5))
pairs(data,cex=0.8,col=as.numeric(bundestag$ewb)+3)

table(out.dbs$cluster,bundestag$ewb)
