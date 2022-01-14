### EXERCISES 1-4 - MS&BDA 
### ZURLO GIOVANNI 03/10/2021

## EXERCISE 1 #### 

library(pdfCluster)
data("oliveoil")
str(oliveoil)
olivex=oliveoil[,3:10]
solive=scale(olivex)

# UNSCALED K=3
set.seed(1234)
olive3=kmeans(olivex, 3, nstart = 100, iter.max = 100)
pro=table(olive3$cluster,oliveoil$macro.area)

# SCALED K=3
set.seed(1234)
solive3=kmeans(solive, 3, nstart = 100, iter.max = 100)
table(solive3$cluster,oliveoil$macro.area)

# UNSCALED K=9
set.seed(1234)
olive9=kmeans(olivex, 9, nstart = 100, iter.max = 100)
table(olive9$cluster,oliveoil$region)

# SCALED K=9
set.seed(1234)
solive9=kmeans(solive, 9, nstart = 100, iter.max = 100)
table(solive9$cluster,oliveoil$region)

#####

## EXERCISE 3 ####

boston<-read.table(file.choose(),header = T,sep="")
str(boston)
boston2=boston[,-4]
sboston<-scale(boston2)

# Pair plot on principal components
pca.boston=princomp(sboston)
# Most variability is explained by the first 3 components
biplot(pca.boston, cex=.5)
# Not very useful...
pca.scores=pca.boston$scores[,1:4]
# Jointly explain 80% of the standardize data sample
library(fpc)
pairs(pca.scores, cex=1.2, col=(boston$chas+8), pch=(boston$chas+1),
      main='Pair Plot on PC Scores', sub='First 4 PC explain 79% of total sample variability')
unique(boston$rad)
pairs(pca.scores, cex=1, col=(boston$rad), pch=clusym[boston$rad])
# symbol n represents districts with index of accessibility to radial high = 24

hierarc=hclust(method = 'ward.D2', dist(sboston))
plot(hierarc)
abline(h=25, col=2, lwd=2)

library(factoextra)
par(mfrow =c(2,2))
set.seed(1234)
fviz_nbclust(sboston, kmeans, method = "wss", k.max = 15)
fviz_nbclust(sboston, kmeans, method = "silhouette", k.max = 15)

set.seed(1234)
ksboston=kmeans(sboston,centers = 5,nstart = 100)
sboston_clu5=ksboston$cluster
table(ksboston$cluster,boston$rad)

centers=matrix(NA,5,13)
for (i in 1:5) {
centers[i,]=as.vector(colMeans(boston2[sboston_clu5==i,]))  
}
centers=data.frame(centers, row.names = 1:5)
colnames(centers)=colnames(boston2)
# First and third clusters represent districts with index of access rad = 24 
table(ksboston$cluster,boston$rad)
table(ksboston$cluster,boston$chas)
colMeans(boston2[boston2$rad==24,])
summary(boston2[boston2$rad==24,])

pairs(pca.scores, cex=1, col=(sboston_clu5), pch=clusym[sboston_clu5],
main='Plot Matrix on PCA Scores (pch = clusters)')

# CLUSGAP
prova=gapnc(sboston,K.max = 10,B=100, iter.max=50,nstart=90)
prova
prova$gapout
plot(prova$gapout,main = 'ClusGap Plot')

# CLUS PLOT
library(cluster)
clusplot()
#####

## EXERCISE 4####

library(pracma)
kmpp <- function(X, k) {
  n <- nrow(X)
  C <- numeric(k)
  C[1] <- sample(1:n, 1)
  for (i in 2:k) {
    dm <- distmat(X, X[C, ])
    pr <- apply(dm, 1, min); pr[C] <- 0
    C[i] <- sample(1:n, 1, prob = pr)
  }
  kmeans(X, X[C, ])
}

set.seed(1234)
solive9=kmeans(solive, 9, nstart = 1, iter.max = 100)

set.seed(1234)
kmpp.solive9=kmpp(solive,9)
(c(solive9$tot.withinss,kmpp.solive9$tot.withinss))

set.seed(1234)
sboston5=kmeans(sboston, 5, nstart = 1, iter.max = 100)

set.seed(1234)
library(pracma)
kmpp.sboston5=kmpp(sboston,5)
(c(sboston5$tot.withinss,kmpp.sboston5$tot.withinss))

#####

