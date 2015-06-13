
#  Required libraries
library(rjags)
library(ggmcmc)
library(ggplot2)
library(ggthemes)
library(pander)
library(Cairo)
library(plyr)
library(MASS)
library(scales)
library(plyr)
require(gdata)
require(runjags)
require(gdata)

data.1= read.table(file="..//data//multinomial_cat.dat",header=TRUE,na.strings = "",sep="\t")
data.1$Galtype2<-trim(data.1$Galtype2)


galtype<-match(data.1$Galtype2,c("E","E/S0","S","S0","Im"))
Ntype<-length(unique(data.1$Galtype2))


typeSne<-match(trim(data.1$SNtype),c("Ia","Ib","Ib/c","Ic","II","IIn"))
bar<-as.numeric(data.1$bar)-1
jags.data <- list(Y= typeSne,
                  N = nrow(data.1),
                  mag_g = data.1$mag_g,
                  galtype = galtype,
                  bar=bar,
                  Ntype=Ntype
)

model<-"model{
## priors
#tau.R<-pow(sdBeta,-1)
#sdBeta ~ dgamma(0.001,0.001)
# Random intercept 
#for (k in 1:Ntype){
#ranef[k]~ddexp(0,tau.R)
#}



for(k in 2:6){
for(j in 1:2){
beta[j,k]~dnorm(0,1e-5)
}
}

## prior for coefficients
for(j in 1:2){
beta[j,1]~dnorm(0,1e6)
}


## Likelihood
for(i in 1:N){
Y[i]~dcat(p[i,1:6])
    for(j in 1:6){
#z[i,j]<-beta[j,1]+ranef[galtype[i]]
z[i,j]<-beta[1,j]+beta[2,j]*bar[i]
expz[i,j]<-exp(z[i,j])
p[i,j]<-expz[i,j]/sum(expz[i,1:Ntype])
                     }

}

}"


params <- c("beta","p")

inits<-function(){list(beta=structure(.Data=c(rep(NA,2),runif(2*(6-1),-1,1)),.Dim=c(2,6)))}


inits1=inits()
inits2=inits()
inits3=inits()



library(parallel)
cl <- makeCluster(3)
jags.mlogit <- run.jags(method="rjparallel", 
                       data = jags.data, 
#                      inits = list(inits1,inits2,inits3),
                       model=model,
                       n.chains = 3,
                       adapt=2500,
                       monitor=c(params),
                       burnin=20000,
                       sample=40000,
                       summarise=FALSE,
                       plots=FALSE
)

jagssamples <- as.mcmc.list(jags.mlogit)