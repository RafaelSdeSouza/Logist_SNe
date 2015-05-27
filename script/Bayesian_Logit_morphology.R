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

#Read dataset

data.0= read.fwf(file="..//data//snsdss.dat",width = c(7,1, 11, 9, 1,9,1,8,2,1,
                                                       11,1,1,26,2,14,10,9,10,1,
                                                       1,1,2,6,5,5,6,6,6,7,7))


# Select few variables for test. For now, galaxy morphology and SN type

SN_cat<-data.frame(SNtype=data.0[,11],Galtype=data.0[,19])

SN_cat2<-na.omit(SN_cat)
require(gdata)
SN_cat2$SNtype<-trim(SN_cat2$SNtype)

SN_cat3<-SN_cat2[which(SN_cat2$SNtype=="Ia"|SN_cat2$SNtype=="II"),]
SN_cat3$SNtype<-droplevels(SN_cat3$SNtype)

# Start the logit model

# Define data for JAGS
#X<-model.matrix(~ M_gal + log_sSFR + g.r+
#                r.i+i.z, data = SNedata3)

X<-model.matrix(~ Galtype, data = SN_cat3)
  
K<-ncol(X)

typeSne<-as.numeric(SN_cat3$SNtype)-1
jags.data <- list(Y= typeSne,
                 N = nrow(SN_cat3),
                 X = X,
                 b0 = rep(0, K),
                 B0 = diag(0.00001, K))


model<-"model{
#1. Priors
beta ~ dmnorm(b0[], B0[,])
#2. Likelihood
for (i in 1:N){
Y[i] ~ dbern(p[i])
logit(p[i]) <-  eta[i]
eta[i] <- inprod(beta[], X[i,])
#3. Prediction
prediction[i]~dbern(p[i])
}
}"
inits<-function () {
  list(beta = rnorm(K, 0, 0.1))}
params <- c("beta","prediction")

jags.logit<-jags.model(
  data = jags.data, 
  inits = inits(), 
  textConnection(model),
  n.chains = 3,
  n.adapt=1000
  )
update(jags.logit, 20000)
posterior.logit <- coda.samples(jags.logit, params, n.iter = 50000)

require(ggmcmc)
beta_post<-ggs(posterior.logit ,family=c("beta"))
ggs_density(beta_post)

jagssamples <- jags.samples(jags.logit, params, n.iter = 50000)
predtype<-summary(as.mcmc.list(jagssamples$prediction))
predtype<-predtype$quantiles

require(mlearning)
require(caret)
SNe_conf<-confusion(predtype[,3], typeSne)

xtab <- table(predtype[,3], typeSne)
confusionMatrix(xtab)

confusionMatrix(xtab, prevalence = 0.25) 