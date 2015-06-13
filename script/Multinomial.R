
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
                  Ntype=Ntype,
                  b0 = rep(0, 2),
                  B0 = diag(0.00001, 2)
)

model<-"model{
## priors
#tau.R<-pow(sdBeta,-1)
#sdBeta ~ dgamma(0.001,0.001)
# Random intercept 
#for (k in 1:Ntype){
#ranef[k]~ddexp(0,tau.R)
#}

## prior for coefficients
for(k in 1:2){
beta[1,k]<-0
}
for(j in 2:6){
beta[j,1:2]~dmnorm(b0[], B0[,])
}

## Likelihood
for(i in 1:N){
    for(j in 1:6){
#z[i,j]<-beta[j,1]+ranef[galtype[i]]
z[i,j]<-beta[j,1]
expz[i,j]<-exp(z[i,j])
p[i,j]<-expz[i,j]/sum(expz[i,1:Ntype])
                     }
Y[i]~dcat(p[i,1:6])
}

}"


params <- c("beta","ranef","p")

inits1=list(beta.0=rnorm(1,0,1),beta.1=rnorm(1,0,1),beta.2=rnorm(1,0,1))
inits2=list(beta.0=rnorm(1,0,1),beta.1=rnorm(1,0,1),beta.2=rnorm(1,0,1))
inits3=list(beta.0=rnorm(1,0,1),beta.1=rnorm(1,0,1),beta.2=rnorm(1,0,1))

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