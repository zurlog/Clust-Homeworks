# Estimate parameters of 1d finite normal mixture model
# with 2 components for Waiting time observations from
# the Old Faithful data set:

library(mclust)
waiting <- faithful$waiting
n <- length (waiting)
waiting.Mclust <- Mclust (waiting ,model="V",G=2)

# Plot densities:
x <- seq (from=min(waiting),to=max(waiting),length=1000)
den1 <- dnorm (x,mean=waiting.Mclust$parameters$mean[1],
               sd=sqrt(waiting.Mclust$parameters$variance$sigmasq[1]))
den2 <- dnorm (x,mean=waiting.Mclust$parameters$mean[2],
               sd=sqrt(waiting.Mclust$parameters$variance$sigmasq[2]))
tau1 <- waiting.Mclust$parameters$pro [1]
tau2 <- waiting.Mclust$parameters$pro [2]
dens <- tau1*den1 + tau2*den2
plot (x,dens ,type="l",xlab="y",ylab="Probability Density",
      ylim=c(-max(dens)/10,1.05*max(dens)),
      main="Density for 1-dim 2-component normal mixture model",lwd =2)
lines (x,tau1*den1 ,col="red")
lines (x,tau2*den2 , col="blue")
legend (x=min(x),y=max(dens),legend=c("Mixture","Component 1","Component 2"),
        col=c("black","red","blue"),lty=c(1,1,1),lwd=c(2,1,1))

# The density is clearly bimodal, lower in the middle between the two humps
# so there is some separation between the two mixture components
# separation is not total, thereis some overlap between the points from the two

# Simulate and plot a sample from the distribution:
sim.results <- simV (waiting.Mclust$parameters ,n,seed=0)
ysim <- sim.results [,2]
groupsim <- sim.results[,"group"]
ysim1 <- ysim[groupsim ==1]
ysim2 <- ysim[groupsim ==2]
points (ysim1 ,rep(-max(dens)/20,length(ysim1)),col="red",pch =19)
points (ysim2 ,rep(-max(dens)/20,length(ysim2)),col="blue",pch =19)

