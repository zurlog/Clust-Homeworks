## CHEATSHEET ####
## MODERN STATISTICS AND BIG DATA ANALYSIS ####

## LIBRARIES ####
library(pdfCluster)
library(fpc)
library(factoextra)
library(cluster)
library(smacof)
#####

set.seed(1234)
kmeans(x, centers=k, iter.max = 100, nstart = 100, trace=FALSE)
NbClust::NbClust(data = NULL, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 20, 
                 method = NULL)

## GAP STATISTIC ####
cluster::clusGap(x, FUNcluster, K.max, B = 100, d.power = 2,
                 spaceH0 = c("scaledPCA", "original"),
                 SE.factor = 2, method="globalSEmax",
                 nstart=100,...)

print(obj, method = "globalSEmax", SE.factor = 2, ...)
fviz_gap_stat(obj)

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
# The output of clusGap is in component gapout.
# The optimal number of clusters is in component nc.
# The optimal kmeans output is in component kmopt.

set.seed(1234)
factoextra::fviz_nbclust(x, FUNcluster = hcut, k.max = 20, hc_method = "ward.D2",
             method = "gap_stat", maxSE = list(method = "globalmax", SE.factor = 2))

# This can be used to compute S for any clustering:
fpc::cluster.stats(dist(p05),kmbundestag5$cluster)
kmb$within.cluster.ss
# This is the same as kmeans5$tot.withinss

#####

## DISSIMILARITIES ####

scale(x)
dist(x, method = c("euclidean","manhattan"))
as.matrix(dist)[1,2]  # Check a dissimilarity value

cluster::daisy(x, metric = "euclidean")  # Handles NAs

# The mahalanobis command can only compute a vector of Mahalanobis distances
mahalm <- matrix(0,ncol=572,nrow=572)
olivecov <- cov(olive)
for (i in 1:572){
  mahalm[i,] <- mahalanobis(olive,as.numeric(olive[i,]),olivecov)}
# Note that it doesn't make a difference whether the data set is scaled or not.

dist(veronica,method="binary")  # Jaccard (asymmetric)
dist(veronica,method="manhattan")/583  # SMC for dummies
nomclust::sm(x)

1-abs(cor(x)) #  Largest dissim. for r = 0
0.5-cor(x)/2  #  Largest dissim. for r = -1
as.dist(cordist)

daisy(housing, metric="gower", type=list(asymm=c(2,4), symm=c(3,6)))

#####

## AAHC ####
# Complete Linkage enforces within-cluster homogeneity ag. separation
# Single Linkage enforces between-cluster separation ag. homogeneity
hclust(diss, method = "average")
cutree(hclust, k)

plot(hclust, hang=-1, xlab="",sub="" ,cex=0.6, cex.axis=1)
b=ggdendrogram(out.eucl, rotate = F, theme_dendro = F,labels = F)
print(b + ggtitle("Ward"))

tail(cbind(hclust$merge, hclust$height),13)
plot(3:20,rev(hclust$height)[3:20],type="b",
     xlab="K",ylab="Height",
     cex.lab=0.8,cex.axis=0.7,
     main="Scree Plot", cex.main=0.9); grid()

plot(x=1:10,y=rev(tail(avg$height,10)),type="b", xlab="K",ylab="Height",cex.lab=0.8,cex.axis=0.7,
main="Scree Plot", cex.main=0.9); grid()


#####

## SILHOUETTE ####
# Suited for PAM and Average Linkage
fviz_nbclust(x, FUNcluster = hcut, k.max = 20, hc_method = "complete",
             hc_metric="euclidean", method = "silhouette")

pasw <- NA
pclusk <- list()
psil <- list()
# Look at K between 2 and 30:
for (k in 2:30){
  # PAM clustering:
  pclusk[[k]] <- pam(diss,k)
  # Computation of silhouettes (partition vector as argument):
  psil[[k]] <- silhouette(pclusk[[k]],dist=diss)
  # ASW needs to be extracted:
  pasw[k] <- summary(psil[[k]])$avg.width}

plot(1:30,pasw,type="l",xlab="Number of clusters",ylab="ASW Plot")
plot(psil[[5]])
#####

## VISUALIZATION ####
## PCA is connected to varcov and euclidean dist so good for kmeans
## Made for continuous variables - struggle with high dimensional (tiny % of variation)
princomp(x)
plot(x=pca$scores[,1], y=pca$scores[,2], 
     col = ($cluster+1), pch=($cluster+1),
     cex=1.1, xlab="PC1 (44.3%)", ylab = "PC2 (27.8%)", main = "PC1 / PC2 - plot")

clusym[]; clucols(mod$cluster)

smacof::mds(diss, ndim = 2, type = "ratio")
plot(mds$conf, type = "p", asp = 1, main="MDS Plot -"
     xlab=paste("Stress : ", as.character(mds$stress*100),"%"))
# Intended to represent diss. in a euclidean way, as good as possible
# Artifacts of MDS may show results in a misleading way (depending on stress)
# Good level of stress: under 10% 
# Over 20% there is quite a bit of info not represented

# I decided to use the Jaccard distance also for AFLP genes:
vveronica <- dist(t(x),method="binary")
varclust <- hclust(vveronica,method="average")
# As a clustering this is pretty messy,
# but still it can be used to impose an order of genes.

heatmap(as.matrix(x),Rowv=as.dendrogram(average),
        Colv=as.dendrogram(varclust),
        col=grey(seq(1,0,-0.01)))   # OR

# heatmap, rows ordered by clusters,
# columns by earlier variable clustering
heatmap(veronicam[order(veronicabernm$flexout[[6]]@cluster),],
        Rowv=NA,Colv=as.dendrogram(varclust),
        RowSideColors=palette()[veronicabernm$flexout[[6]]@cluster]
        [order(veronicabernm$flexout[[6]]@cluster)],
        col=c(0,1),scale="none")

clusplot(pam)
#####

## MODEL BASED CLUSTERING ####

mclust::Mclust(x,G=1:15, modelNames = c("VVV"))
summary(mclust$BIC)
plot(mclust)


molive$classification
# Clustering vector
molive$parameters
# Estmated parameters
molive$z
# Matrix of posterior probabilities p_ik that point i was generated
# by mixture component k
prod(mclust$parameters$variance$shape)==1

factoextra::fviz_mclust_bic(mclust, model.names = NULL, shape = 1, lwd=3,
                            color = "model", palette = NULL, legend = NULL,cex=2,
                            main = "Mixture Model Selection", xlab = "Number of Components",
                            ylab = "BIC")

### To have a look at the best covmat models, just rerun the function
selected_mix=Mclust(swdbcc, G=1:10, verbose = T, modelNames = c('VVV','VEV','EVV'))

#####

## MIXTURES OF SKEW AND HEAVY-TAILED ####
# Different mixtures fit data in different ways and it's hard to say what's best
# Mixtures can well approximate each other, it's hard to choose among them
# Skew shapes can be interpreted as one cluster or as one symmetric core plus others

library(EMMIXskew)
for (i in 1:12){
  print(i)
  tryattempts <- 3
  trycounter <- 1
  tst <- try(skewmix[[i]] <- EmSkew(x,g=i,distr="mst",ncov=3))
  while((is.null(tst) | class(tst)=="try-error") & trycounter<tryattempts+1){
    print("Error, try again")
    tst <- try(skewmix[[i]] <- EmSkew(x,g=i,distr="mst",ncov=3))
    trycounter <- trycounter+1
  }
  trycounter <- 1
  while((is.null(tst) | class(tst)=="try-error") & trycounter<tryattempts+1){
    print("Error, try again")
    tst <- try(skewmix[[i]] <- EmSkew(x,g=i,distr="mst",ncov=4))
    trycounter <- trycounter+1
  }
  trycounter <- 1
  while((is.null(tst) | class(tst)=="try-error") & trycounter<tryattempts+1){
    print("Error, try again")
    tst <- try(skewmix[[i]] <- EmSkew(x,g=i,distr="mst",ncov=2))
    trycounter <- trycounter+1
  }
  bicvals[i] <- skewmix[[i]]$bic
  #ariarea[i] <- adjustedRandIndex(skewmix[[i]]$clust,oliveoil$macro.area)
  #ariregion[i] <- adjustedRandIndex(skewmix[[i]]$clust,oliveoil$region)
}



#####

## LATENT CLASS ANALYSIS ####
## For categorical data (also binary)
## Previously we used AHC on SM-J Dissimilarities + MDS

set.seed(1234)
fpc::flexmixedruns(x, continuous=0, discrete=ncol(x), n.cluster=1:10, simruns = 100,
                   verbose = T, allout = F)

which.min(out$bicvals) # Which 1:10 model is the best?
out$optimalk; out$optsummary
plot(1:10,out$bicvals,typ="l", xlab="Number of clusters", ylab="BIC")

str(out$flexout, max.level = 2)
# if allout=TRUE, flexout[[]] list of flexmix output objects for all numbers of components
out$flexout[[k]]@cluster # Clustering
out$flexout[[k]]@prior # Mixing parameters
str(out$flexout[[k]]@components) # Model Object Structure
out$flexout[[k]]@components$Comp.1[[1]]@parameters$pp # zeta par for 1st Comp
out$flexout[[k]]@components[[5]][[1]]@parameters$pp




f <- cbind(x1,x2,x3,x4)~1
set.seed(1234)
poLCA::poLCA(f, x, nclass=3, maxiter = 5000, nrep = 70)
# nclass K is fixed in advance
a$predclass # Clustering
a$bic
#####

## ROBUST STATISTICS ####

mad(x)
robustbase::huberM(x,k=1.5)
huber$s    # MAD by default
huber$SE   # standard error
robustbase::covMcd(x, alpha = 0.75)
# Do not use RAW components

# library robustbase has a plot.mcd function and one could
# use plot(mcdd) for outlier diagnostic plots, but this has some problems
# - need to add tol=1e-20 because otherwise gives an error.
plot(1:nrow(x),sqrt(mcdd$mah),type="n",xlab="Observation",
     ylab="Squared robust Mahalanobis distance")
text(1:170,sqrt(mcdd$mah),rownames(dortmund),cex=0.7)
abline(sqrt(qchisq(0.99,7)),0,col=2)

plot(sqrt(mcdd75$mah),sqrt(mcdd$mah),xlim=c(0,30),ylim=c(0,30),
     xlab="Squared robust Mahalanobis distance (alpha=0.75)",
     ylab="Squared robust Mahalanobis distance (alpha=0.5)")
abline(sqrt(qchisq(0.99,7)),0,col=2)
abline(v=sqrt(qchisq(0.99,7)),col=2)

robustbase::lmrob(y~x1+x2+x3, method="MM",data=regdata3)   # MM-estimator
par(mfrow=c(2,3))
plot(mm)
plot(1:nrow(x),mm$rweights)
summary(mm)

#####


