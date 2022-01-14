### EXERCISES 7-4 - MS&BDA 
### ZURLO GIOVANNI 29/11/2021

## FLEXMIXEDRUNS TEST
library(flexmix)
library(fpc)
library(poLCA)
library(nomclust)
library(smacof)
library(mclust)

data(carcinoma)
?carcinoma
head(carcinoma)

out=flexmixedruns(carcinoma,discrete = ncol(carcinoma),continuous = 0,simruns = 50, verbose = 
                    T,allout = T, n.cluster =1:10)

out$bicvals # Which 1:10 model is the best?
which.min(out$bicvals)

plot(1:10,out$bicvals,typ="l",
     xlab="Number of clusters", ylab="BIC")

str(out$flexout[[6]], max.level = 2)
out$flexout[[6]]@prior # Mixing parameters
str(out$flexout[[6]]@components) # Response prob/zeta parameters
out$flexout[[6]]@components[[5]][[1]]@parameters$pp





## EXERCISE 1 ####
library(poLCA)
data(election)
?election
# Two sets of six questions with four responses, asking opinions about various traits (moral, caring, dishonest, intelligent)
# describing presidential candidates Al Gore and George W. Bush
# The responses are (1) Extremely well; (2) Quite well; (3) Not too well; (4) Not well 
# + 5 covariates in the end

election12 <- election[,1:12]
electioncomplete <- election12[complete.cases(election12),] # 1311 obs.

electionwithna <- election12 
for (i in 1:12){
  levels(electionwithna[,i]) <- c(levels(election12[,i]),"NA") 
  electionwithna[is.na(election12[,i]),i] <- "NA"
} # 1785 obs

nomclust::sm.complete=sm(electioncomplete)
nomclust::sm.withna=sm(electionwithna)
mds.complete=mds(sm.complete,type = 'ratio')
mds.complete$stress # 29.8%
mds.withna=mds(sm.withna,type = 'ratio') # 31%

# ARI Comparison restricted to the 1311 obs
# clustering_to_compare = clustering_with_na[complete.case(election12)]

## a) ####
f <- cbind(MORALG,CARESG,KNOWG,LEADG,DISHONG,INTELG,
           MORALB,CARESB,KNOWB,LEADB,DISHONB,INTELB)~1
set.seed(1234)
a <- poLCA(f,election12,nclass=3, maxiter = 5000, nrep = 70) #BIC: 34218
# Uses complete cases, equivalent to running it on electioncomplete

## b) ####
set.seed(1234)
b <- poLCA(f,electionwithna,nclass=3,maxiter = 5000, nrep = 50) #BIC: 52863

b.c=b$predclass[complete.cases(election12)]

par(mfrow=c(1,2))
plot(mds.complete$conf[,1],mds.complete$conf[,2],col=a$predclass,cex=0.75,pch=20,xlab="MDS Stress = 29.8%",ylab="",main="poLCA electioncomplete")
plot(mds.withna$conf[,1],mds.withna$conf[,2],col=b$predclass,cex=0.75,pch=20,xlab="MDS Stress = 31%",ylab="",main="poLCA electionwithna")

## c) ####
set.seed(1234)
c=flexmixedruns(electioncomplete, discrete = 12 ,continuous = 0,simruns = 70, verbose = 
                    T,allout = T, n.cluster =3, control = list(iter.max = 5000)) # BIC: 34219

c$flexout[[3]]
d$flexout[[3]]@cluster

## d) ####
set.seed(1234)
d=flexmixedruns(electionwithna, discrete = 12 ,continuous = 0,simruns = 70, verbose = 
                  T,allout = T, n.cluster =3, control = list(iter.max = 5000)) # BIC: 52863

d.c=d$flexout[[3]]@cluster[complete.cases(election12)]

par(mfrow=c(1,2))
plot(mds.complete$conf[,1],mds.complete$conf[,2],
     col=c$flexout[[3]]@cluster,cex=0.75,pch=20,
     xlab="MDS Stress = 29.8%",ylab="",main="flexmixedruns electioncomplete")
plot(mds.withna$conf[,1],mds.withna$conf[,2],
     col=d$flexout[[3]]@cluster, cex=0.75,pch=20,
     xlab="MDS Stress = 31%",ylab="",main="flexmixedruns electionwithna")
# cluster ordering

## e) ####
set.seed(1234)
e=pam(sm.complete,diss = T, k=3, nstart = 100, trace.lev = 1)

## f) ####
set.seed(1234)
f.mod=pam(sm.withna,diss = T, k=3, nstart = 100, trace.lev = 1)

f.c=f.mod$clustering[complete.cases(election12)]

par(mfrow=c(1,2))
plot(mds.complete$conf[,1],mds.complete$conf[,2],
     col=e$clustering,cex=0.75,pch=20,
     xlab="MDS Stress = 29.8%",ylab="",main="PAM electioncomplete")
plot(mds.withna$conf[,1],mds.withna$conf[,2],
     col=f$clustering, cex=0.75,pch=20,
     xlab="MDS Stress = 31%",ylab="",main="PAM electionwithna")

## g) ####

gower=daisy(election12,metric = 'gower')
# 'dissimilarity' num [1:1592220]
attributes(gower)$NA.message
head(which(is.na(gower)))
gower[which(is.na(gower))]=0
gower=as.matrix(gower)
rownames(gower)=as.integer(rownames(gower))
#index=as.integer(rownames(gower[which(is.na(gower)),]))
index=rownames(gower[rowSums(is.na(gower)) == 0, colSums(is.na(gower)) == 0, drop = FALSE])
str(index) #OK
# BUT The dissimilarity matrix CANNOT be calculated if
# the 'data' argument contains NA values.

gower.mod = as.matrix(gower)
gower.mod = gower.mod[rowSums(is.na(gower.mod)) == 0, colSums(is.na(gower.mod)) == 0, drop = FALSE]
gower.mod = as.dist(gower.mod)
# 'dist' Named num [1:1073845]
# dissimilarities for 1466 are maintained (smaller loss of info)

mds.gow=mds(gower.mod,type = 'ratio')
mds.gow$stress

plot(mds.gow$conf[,1],mds.gow$conf[,2],
     col=g$clustering,cex=0.95,pch=20,
     xlab="MDS Stress = 29.8%",ylab="",main="PAM (3) - Gower DissMat Sample")


set.seed(1234)
g=pam(gower.mod,diss = T,k=3, nstart = 150, trace.lev = 1)
g.c=g$clustering[complete.cases(election12)]

## h) ####
set.seed(1234)
h=flexmixedruns(electionwithna, discrete = 12 ,continuous = 0,simruns = 50, verbose = 
                  T,allout = T, n.cluster=1:10, control = list(iter.max = 5000))

h$optsummary # K = 8 is optimal

plot(1:10,h$bicvals,typ="l",xlab="Number of clusters",ylab="BIC")
h.c=h$flexout[[8]]@cluster[complete.cases(election12)]

plot(mds.withna$conf[,1],mds.withna$conf[,2], 
     col=h$flexout[[8]]@cluster, cex=1,pch=clusym[h$flexout[[8]]@cluster],
     xlab="MDS Stress = 31%",ylab="",main="flexmixedruns (K = 8) electionwithna")


## COMPARISON ####

results=list()
results[[1]]=a$predclass
results[[2]]=b.c
results[[3]]=c$flexout[[3]]@cluster
results[[4]]=d.c
results[[5]]=e$clustering
results[[6]]=f.c
results[[7]]=rep(NA,1311)
results[[8]]=h.c

labels=c('a','b','c','d','e','f','g','h')
results=setNames(object=results,labels)

ARI=matrix(NA,nrow = 8,ncol = 8, dimnames = list(labels,labels))

for (i in 1:8) {
  for (j in 1:8) {
    ARI[i,j]=adjustedRandIndex(results[[i]],results[[j]])
  }
}

library(kableExtra)
ARI |>
kable("html", caption = 'Clusterings Pairwise Comparisons (ARI)', row.names = T) |>
  kable_styling(full_width = F, position = "center")






## EXERCISE 2 ####
#vveronica <- dist(t(veronicam),method="binary")
library(nomclust)

election.t=as.data.frame(t(electionwithna))
#f=factor(c("NA", "4 Not well at all","3 Not too well","2 Quite well","1 Extremely well" ))
levels(f)<-list("NA"=0,"4 Not well at all"=1, "3 Not too well"=2, "2 Quite well"=3,"1 Extremely well"=4)


election.t[,2]=as.factor(election.t[,2])
levels(election.t[,2])<-list("NA"=0,"4 Not well at all"=1, "3 Not too well"=2, "2 Quite well"=3,"1 Extremely well"=4)

for (i in 1:ncol(election.t)) {
  election.t[,i]=as.factor(election.t[,i])
  levels(election.t[,i])<-list("NA"=1,"4 Not well at all"=2, "3 Not too well"=3, "2 Quite well"=4,"1 Extremely well"=5)
}

j.vars=dist(election.t,method="binary")
j.vars=sm(election.t)
varclust <- hclust(j.vars,method="average")
# heatmap, rows ordered by clusters,
# columns by earlier variable clustering

heatmap.df=electionwithna
for (i in 1:ncol(heatmap.df)) {
  levels(heatmap.df[,i])<-list("NA"=1,"4 Not well at all"=2, "3 Not too well"=3, "2 Quite well"=4,"1 Extremely well"=5)
  heatmap.df[,i]=as.integer(heatmap.df[,i])
}

heatmap.df=as.matrix(heatmap.df)

library(RColorBrewer)
  # col=colorRampPalette(brewer.pal(8, "Blues"))(25)
heatmap(x=heatmap.df[order(f.mod$clustering),],
        Rowv=NA,
        Colv=as.dendrogram(varclust),
        RowSideColors=palette()[f.mod$clustering][order(f.mod$clustering)],
        col=c('lavenderblush3','darkorange4','darkorange2','olivedrab3','olivedrab4'),
        scale="none") # scale column
legend(x='right',cex=1,
       legend =c("NA","Not well at all","Not too well","Quite well","Extremely well"),
       fill =c('lavenderblush3','darkorange4','darkorange2','olivedrab3','olivedrab4') )


## EXERCISE 3####
# Five variables are binary, three variables have three categories
# and two variables have five categories
p=10
m=c(rep(2,5),rep(3,3),rep(5,2))
sum(m-1)
prod(m-1) #128 df
# would need at least 13*p sample


mj=rep(2,3)
mj

bin=c(0,1)
obs=expand.grid(bin,bin,bin)
p_y.1=rep(0.1,3)
p_y.2=rep(0.5,3)
p_y.3=rep(0.8,3)

######
g_x.y1=numeric()
g_x.y1[1]=1/3*(0.9)^3
g_x.y1[2]=1/3*0.1*(0.9)^2
g_x.y1[3]=1/3*0.1*(0.9)^2
g_x.y1[4]=1/3*((0.1)^2)*0.9
g_x.y1[5]=1/3*0.1*(0.9)^2
g_x.y1[6]=1/3*((0.1)^2)*0.9
g_x.y1[7]=1/3*((0.1)^2)*0.9
g_x.y1[8]=1/3*(0.1)^3

g_x.y2=numeric()
g_x.y2[1]=1/3*(0.5)^3
g_x.y2[2]=1/3*(0.5)^3
g_x.y2[3]=1/3*(0.5)^3
g_x.y2[4]=1/3*(0.5)^3
g_x.y2[5]=1/3*(0.5)^3
g_x.y2[6]=1/3*(0.5)^3
g_x.y2[7]=1/3*(0.5)^3
g_x.y2[8]=1/3*(0.5)^3

g_x.y3=numeric()
g_x.y3[1]=1/3*(0.2)^3
g_x.y3[2]=1/3*0.8*(0.2)^2
g_x.y3[3]=1/3*0.8*(0.2)^2
g_x.y3[4]=1/3*((0.8)^2)*0.2
g_x.y3[5]=1/3*0.8*(0.2)^2
g_x.y3[6]=1/3*((0.8)^2)*0.2
g_x.y3[7]=1/3*((0.8)^2)*0.2
g_x.y3[8]=1/3*(0.8)^3

f_x=(g_x.y1+g_x.y2+g_x.y3)

probs=round(cbind(obs,g_x.y1,g_x.y2,g_x.y3,f_x),4)
colnames(probs)=c("X1","X2","X3","g(x|y=1)","g(x|y=2)","g(x|y=3)","f(x)")

library(kableExtra)
probs |>
  kable("html", caption = 'Components Contribution & f(x)', row.names = F) |>
  kable_styling(full_width = F, position = "center")
######

g_x.y1=numeric()
g_x.y1[1]=(0.9)^3
g_x.y1[2]=0.1*(0.9)^2
g_x.y1[3]=0.1*(0.9)^2
g_x.y1[4]=((0.1)^2)*0.9
g_x.y1[5]=0.1*(0.9)^2
g_x.y1[6]=((0.1)^2)*0.9
g_x.y1[7]=((0.1)^2)*0.9
g_x.y1[8]=(0.1)^3

g_x.y2=numeric()
g_x.y2[1]=(0.5)^3
g_x.y2[2]=(0.5)^3
g_x.y2[3]=(0.5)^3
g_x.y2[4]=(0.5)^3
g_x.y2[5]=(0.5)^3
g_x.y2[6]=(0.5)^3
g_x.y2[7]=(0.5)^3
g_x.y2[8]=(0.5)^3

g_x.y3=numeric()
g_x.y3[1]=(0.2)^3
g_x.y3[2]=0.8*(0.2)^2
g_x.y3[3]=0.8*(0.2)^2
g_x.y3[4]=((0.8)^2)*0.2
g_x.y3[5]=0.8*(0.2)^2
g_x.y3[6]=((0.8)^2)*0.2
g_x.y3[7]=((0.8)^2)*0.2
g_x.y3[8]=(0.8)^3

f_x=(g_x.y1/3 + g_x.y2/3 + g_x.y3/3)

probs=round(cbind(obs,g_x.y1,g_x.y2,g_x.y3,f_x),4)
colnames(probs)=c("X1","X2","X3","g(x|y=1)","g(x|y=2)","g(x|y=3)","f(x)")

library(kableExtra)
probs |>
  kable("html", caption = 'Components Contribution & f(x)', row.names = F) |>
  kable_styling(full_width = F, position = "center")
