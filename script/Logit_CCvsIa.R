# Bayesian  Logistic Regression using JAGS

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
require(runjags)
#Read the already clean dataset

data.0= read.table(file="..//data//clean_cat.dat",header=TRUE)

galtype<-match(SN_cat4$Galtype2,c("S","S0","E/S0","E","Im"))
Ntype<-length(unique(SN_cat4$Galtype2))


typeSne<-as.numeric(SN_cat4$SNtype)-1
jags.data <- list(Y= typeSne,
                 N = nrow(SN_cat4),
                 mag_g = data.0$mag_g,
                 galtype = galtype,
                 Ntype=Ntype
                 )


model<-"model{
#1. Priors

tau.R<-pow(sdBeta,-1)
sdBeta ~ dgamma(0.01,0.01)

# Random intercept 
for (j in 1:Ntype){
ranef[j]~ddexp(0,tau.R)
}

beta.0~dnorm(0,0.001)
beta.1~dnorm(0,0.001)

#2. Likelihood
for (i in 1:N){
Y[i] ~ dbern(p[i])
logit(p[i]) <-  eta[i]
eta[i] <- beta.0+beta.1*mag_g[i]+ranef[galtype[i]]
#3. Prediction
prediction[i]~dbern(p[i])
}
}"

params <- c("beta.0","beta.1","prediction","ranef","p")

inits1=list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1))
inits2=list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1))
inits3=list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1))

library(parallel)
cl <- makeCluster(3)
jags.logit <- run.jags(method="rjparallel", method.options=list(cl=cl),
                     data = jags.data, 
                     inits = list(inits1,inits2,inits3),
                     model=model,
                     n.chains = 3,
                     adapt=2500,
                     monitor=c(params),
                     burnin=20000,
                     sample=40000,
                     summarise=FALSE,
                     thin=5,
                     plots=FALSE
)

jagssamples <- as.mcmc.list(jags.logit)
summary<-extend.jags(jags.logit,drop.monitor=c("prediction"), summarise=TRUE)


require(ggmcmc)
L.factors <- data.frame(
  Parameter=paste("ranef[", seq(1:5), "]", sep=""),
  Label=c("S","S0","E/S0","E","Im"))
head(L.factors)
ranef_post<-ggs(jagssamples,par_labels=L.factors,family=c("ranef"))


pdf("..//figures/ranef.pdf",width=7,height=7)
ggs_caterpillar(ranef_post)+theme_stata()+ylab("")
dev.off()


pdf("..//figures/density2.pdf",width=7,height=10)
ggs_density(ranef_post)+theme_stata()+ylab("")
dev.off()


# Diagnostics Confusion Matrix (very unlikely to be good)
predtype<-summary(as.mcmc.list(jags.logit, vars="prediction"))
predtype<-predtype$quantiles
require(mlearning)
require(caret)


xtab <- table(predtype[,3], typeSne)
confusionMatrix(xtab)

# probabilities
prob<-summary(as.mcmc.list(jags.logit, vars="p"))
prob<-prob$quantiles