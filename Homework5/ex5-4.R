### EXERCISES 5-4 - MS&BDA 
### ZURLO GIOVANNI 15/11/2021

library(mclust)
library(cluster)
library(mclust)
library(smacof)
library(factoextra)
library(dbscan)
library(knitr)
library(ggdendro)
library(fpc)

##  EXERCISE 1 #### 
setwd("C:/Users/giova/Dropbox/Magistrale/[NEW] MS & BD Analytics/Datasets")
wdbc <- read.csv("wdbc.dat",header=FALSE)
wdbcc <- wdbc[,3:12]
swdbcc = scale(wdbcc) 
wdbcdiag <- as.integer(as.factor(wdbc[,2]))
table(wdbcdiag)
# 2 = M     /  1 = B

dist=dist(swdbcc)
pca<- prcomp(swdbcc)
summary(pca) # NOW 80%!!

## KMEANS ####
cg1 <- clusGap(swdbcc,kmeans,10,B=100,d.power=2,spaceH0="scaledPCA",nstart=100,iter.max = 100)
plot(cg1,main="") # 2
print(cg1,method="globalSEmax",SE.factor=2)
print(cg1,method="Tibs2001SEmax",SE.factor=1)
print(cg1,method="firstSEmax",SE.factor=2) 

tasw <- NA
tclusk <- list()
tsil <- list()
set.seed(1234)
for (k in 2:10){
  tclusk[[k]] <- kmeans(swdbcc,centers=k,iter.max = 100,nstart = 200)
  tsil[[k]] <- silhouette(tclusk[[k]]$cluster,dist)
  tasw[k] <- summary(tsil[[k]])$avg.width
}
plot(1:10,tasw,type="b",ylab="ASW",
     xlab=paste("Number of clusters   -   Max ASW =",as.character(round(max(tasw[-1]),3))))
grid()
m=match(max(tasw[-1]),tasw[-1])+1
abline(v=m,col="red",lty="dashed") # ASW = 0.395

kmeans.out=kmeans(swdbcc,centers=2,iter.max = 1000,nstart = 500)
table(kmeans.out$cluster) # balanced
kmeans.f=as.factor(kmeans.out$cluster)

plot(x=pca$x[,1],y=pca$x[,2], cex=1.3, lwd=1.8,
     col=kmeans.out$cluster, pch=1, xlab="",ylab="")
# Fine separation

## AVERAGE - Excluded ####
average <- hclust(dist,method="average")
a=ggdendrogram(average, rotate = F, theme_dendro = F,labels = F)
print(a + ggtitle("Average Linkage") + xlab("") + ylab("Height"))

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:10){
  tclusk[[k]] <- cutree(average,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=dist)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=dist))$avg.width
}
plot(1:10,tasw,type="b",ylab="ASW",
     xlab=paste("Number of clusters   -   Max ASW =",as.character(round(max(tasw[-1]),3))))
k.average=match(max(tasw[-1]),tasw[-1])+1
abline(v=k.average,col="red",lty="dashed") # ASW (2) = 0.483
# But k=2 just isolates 17 low density obs - tendency to select k=2
# Unable to separate those density clusters

average.out = cutree(average,7)
table(average.out)
average.f=as.factor(average.out)

plot(x=pca$x[,1],y=pca$x[,2], cex=1.3, lwd=1.8,
     col=average.out, pch=1, xlab="",ylab="")



## WARD ####
ward <- hclust(dist,method="ward.D2")
a=ggdendrogram(ward, rotate = F, theme_dendro = F,labels = F)
print(a + ggtitle("ward Linkage") + xlab("") + ylab("Height"))

tasw <- NA
tclusk <- list()
tsil <- list()
set.seed(1234)
for (k in 2:10){
  tclusk[[k]] <- cutree(ward,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=dist)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=dist))$avg.width
}
plot(1:10,tasw,type="b",ylab="ASW",
     xlab=paste("Number of clusters   -   Max ASW =",as.character(round(max(tasw[-1]),3))))
k.ward=match(max(tasw[-1]),tasw[-1])+1
abline(v=k.ward,col="red",lty="dashed") # ASW (2) = 0.33 worse than kmeans
# Would choose also k = 4

ward.out = cutree(ward,2)
table(ward.out) #unbalanced
adjustedRandIndex(ward.out,wdbcdiag)
ward.f=as.factor(ward.out)

plot(x=pca$x[,1],y=pca$x[,2], cex=1.3, lwd=1.8,
     col=ward.out, pch=1, xlab="",ylab="") 
# Correct separation of the two density modes

## PAM ####

pasw <- NA
pclusk <- list()
psil <- list()
set.seed(1234)
for (k in 2:10){
  pclusk[[k]] <- pam(dist,k, nstart = 100)
  # Computation of silhouettes:
  psil[[k]] <- silhouette(pclusk[[k]],dist=dist)
  # ASW needs to be extracted:
  pasw[k] <- summary(psil[[k]])$avg.width
}


plot(1:10,pasw,type="b",ylab="ASW",
     xlab=paste("Number of clusters   -   Max ASW =",as.character(round(max(pasw[-1]),3))))
m=match(max(pasw[-1]),pasw[-1])+1

abline(v=m,col="red",lty="dashed") #ASW = 0.383

pam.f=as.factor(pclusk[[m]]$cluster)
table(pam.f) # 2 leads to balanced classes 

plot(x=pca$x[,1],y=pca$x[,2], cex=1.3, lwd=1.8,
     col=pam.f, pch=1, xlab="",ylab="")
# Quite good


## MIXTURE MODELS ####

mixtures=Mclust(swdbcc, G=1:10, verbose = T) # Fully flexible are the best
summary(mixtures) # VVV 5 components
# Mixtures have the tendency to fit more and more gaussians - misleading
plot(mixtures, which=1, select=1)
mixtures=Mclust(swdbcc, G=2)

bic_vvv=mixtures$BIC[,14]
bic_vev=mixtures$BIC[,12]
bic_evv=mixtures$BIC[,13]

selected_mix=Mclust(swdbcc, G=1:10, verbose = T,
                    modelNames = c('VVV','VEV','EVV'))

library(factoextra)
fviz_mclust_bic(selected_mix, model.names = NULL, shape = 1, lwd=3,
                color = "model", palette = NULL, legend = NULL,cex=2,
                main = "Model selection", xlab = "Number of components",
                ylab = "BIC") ### TOPTOP

VVV=Mclust(swdbcc, G=2, 'VVV')
VEV=Mclust(swdbcc, G=2, 'VEV')
EVV=Mclust(swdbcc, G=2, 'EVV') # Most similar to PAM

# Way lower ASW with G=4 - tendency to fit more clusters 
# Also tendency of ASW to prefer k=2 solution -> relay on visual assess.
summary(silhouette(VVV$classification, dist=dist))$avg.width # .336
summary(silhouette(VEV$classification, dist=dist))$avg.width # .346
summary(silhouette(EVV$classification, dist=dist))$avg.width # .347 
adjustedRandIndex(VVV$classification,wdbcdiag) 
adjustedRandIndex(VEV$classification,wdbcdiag) 
adjustedRandIndex(EVV$classification,wdbcdiag) 

plot(x=pca$x[,1],y=pca$x[,2], cex=1.3, lwd=1.8,
     col=VVV$classification, pch=1, xlab="",ylab="")


tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:10){
  tclusk[[k]] <- Mclust(wdbcc, G=k, verbose = T)
  tsil[[k]] <- silhouette(tclusk[[k]]$classification, dist=dist)
  tasw[k] <- summary(tsil[[k]])$avg.width
}
plot(1:10,tasw,type="b",ylab="ASW",
     xlab=paste("Number of clusters   -   Max ASW =",as.character(round(max(tasw[-1]),3))))
k.mix=match(max(tasw[-1]),tasw[-1])+1
abline(v=k.mix,col="red",lty="dashed") # 0.46 BUT 2
# maybe not appropriate to use silhouette width


# 2D PCA representation contains whole information
plot(x=pca$x[,1],y=pca$x[,2], cex=1.3, lwd=1.8,
     col=VVV$classification, pch=1, xlab="",ylab="")

## ARI COMPARISON ####
ARI<- matrix(0, nrow = 1, ncol = 6)
colnames(ARI)=c("Kmeans (2)","Ward (4)", "PAM (2)","VVV (2)","VEV (2)","EVV (2)")
rownames(ARI)=c("ARI")

ARI[1,1]=adjustedRandIndex(kmeans.out$cluster,wdbcdiag)
#ARI[1,2]=adjustedRandIndex(average.out,wdbcdiag)
ARI[1,2]= adjustedRandIndex(ward.out,wdbcdiag) # with k=4 is 0.463
ARI[1,3]= adjustedRandIndex(pam.f,wdbcdiag) 
ARI[1,4]= adjustedRandIndex(VVV$classification,wdbcdiag) 
ARI[1,5]= adjustedRandIndex(VEV$classification,wdbcdiag) 
ARI[1,6]= adjustedRandIndex(EVV$classification,wdbcdiag) 

library(kableExtra)
ARI |>
  kable("html", caption = 'Clustering Similarity to Diagnoses wdbcdiag measured by ARI', row.names = T) |>
  kable_styling(full_width = T, position = "center")

plot(x=pca$x[,1],y=pca$x[,2], cex=1.3, lwd=1.8,
     col=wdbcdiag, pch=1, xlab="",ylab="") 

