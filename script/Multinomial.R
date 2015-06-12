
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


typeSne<-match(trim(data.1$SNtype),c("Ia","Ib","Ib/c","Ic","II","IIn"))-1
bar<-as.numeric(data.1$bar)-1
jags.data <- list(Y= typeSne,
                  N = nrow(data.1),
                  mag_g = data.1$mag_g,
                  galtype = galtype,
                  bar=bar,
                  Ntype=Ntype
)

model<-"{
## priors

tau.R<-pow(sdBeta,-1)
sdBeta ~ dgamma(0.001,0.001)
  
# Random intercept 
for (j in 1:Ntype){
ranef[j]~ddexp(0,tau.R)
}

for (k in 2:Ntype){ 

## prior for coefficients
for( j in 1:2){
beta[j,k]~dnorm(0,1e6)
}
}
## Likelihood

for(i in 1:N)
Y~dcat(p[i,1:Ntype])
for(k in 1:Ntype){
z[i,k]<-beta[1,k]+ranef[galtype[i]]
expz[i,k]<-exp(z[i,k])
p[i,k]<-expz[i,k]/sum(expz[i,1:Ntype])
}

for(j in 1:2)
{
beta[j,1]<-0
}
  
}"