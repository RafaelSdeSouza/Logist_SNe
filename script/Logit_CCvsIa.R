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
#Read clean dataset

data.0= read.table(file="..//data//clean_cat.dat",header=TRUE)

galtype<-match(SN_cat4$Galtype2,c("S","S0","E/S0","E","Im"))
Ntype<-length(unique(SN_cat4$Galtype2))

X<-model.matrix(~ Galtype2-1, data = SN_cat4)
K<-ncol(X)

typeSne<-as.numeric(SN_cat4$SNtype)-1
jags.data <- list(Y= typeSne,
                 N = nrow(SN_cat4),
                 X = X,
                 b0 = rep(0, K),
                 B0 = diag(0.00001, K))


model<-"model{
#1. Priors
tau.R<-pow(sdBeta,-1)
sdBeta ~ dgamma(0.01,0.01)

for(j in 1:5){
beta[j] ~ ddexp(0,tau.R)
}
beta.0~dnorm(0,1e-6)
#2. Likelihood
for (i in 1:N){
Y[i] ~ dbern(p[i])
logit(p[i]) <-  eta[i]
eta[i] <- beta.0+inprod(beta[], X[i,])
#3. Prediction
prediction[i]~dbern(p[i])
}
}"
inits<-function () {
  list(beta = rnorm(K, 0, 0.1))}
params <- c("beta","prediction")

inits1=inits()
inits2=inits()
inits3=inits()

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
                     sample=30000,
                     summarise=FALSE,
                     thin=5,
                     plots=FALSE
)

jagssamples <- as.mcmc.list(jags.logit)
summary<-extend.jags(jags.logit,drop.monitor=c("prediction"), summarise=TRUE)


require(ggmcmc)
L.factors <- data.frame(
  Parameter=paste("beta[", seq(1:6), "]", sep=""),
  Label=c("beta.0","S","S0","E/S0","E","Im"))
head(L.factors)
beta_post<-ggs(jagssamples,par_labels=L.factors,family=c("beta"))


pdf("..//figures/betas.pdf",width=7,height=7)
ggs_caterpillar(beta_post)+theme_stata()+ylab("")
dev.off()


pdf("..//figures/density.pdf",width=7,height=10)
ggs_density(beta_post)+theme_stata()+ylab("")
dev.off()


# Diagnostics Confusion Matrix (very unlikely to be good)
predtype<-summary(as.mcmc.list(jags.logit, vars="prediction"))
predtype<-predtype$quantiles
require(mlearning)
require(caret)


xtab <- table(predtype[,3], typeSne)
confusionMatrix(xtab)

confusionMatrix(xtab, prevalence = 0.25) 