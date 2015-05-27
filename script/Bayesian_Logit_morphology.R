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

# Host galaxy properties

#Host galaxy log(Mass): logMassFSPS
#Host galaxy age: ageFSPS
#Host galaxy log(SSFR) (FSPS):logSSFRFSPS
#Host galaxy dereddened u magnitude:dereduhost
#Host galaxy dereddened g magnitude:deredghost
#Host galaxy dereddened r magnitude:deredrhost
#Host galaxy dereddened i magnitude:deredihost 
#Host galaxy dereddened z magnitude:deredzhost
#SN-host separation (arcseconds):separationhost

# SNe properties

# Classification: Classification

#Select subset

data.1=data.0[,c("logMassFSPS","ageFSPS","logSSFRFSPS","dereduhost","deredghost",
                      "deredrhost","deredihost","deredzhost","separationhost",
                      "Classification")]
#Let's use colours u-g, g-r, r-i, i-z,

SNedata<-data.frame(M_gal=data.1$logMassFSPS,age=data.1$ageFSPS,
                    log_sSFR=data.1$logSSFRFSPS,u.g=data.1$dereduhost-data.1$deredghost,
                    g.r=data.1$deredghost-data.1$deredrhost,
                    r.i=data.1$deredrhost-data.1$deredihost,
                    i.z=data.1$deredihost-data.1$deredzhost,
                    Sep_host=data.1$separationhost,Type=data.1$Classification)
#Let's use only Type Ia and II for now (no AGNs)

#SNedata2<-SNedata[which(SNedata$Type=="pSNIa" | 
#                SNedata$Type=="pSNII" | SNedata$Type=="SNIa"
#              | SNedata$Type=="SNIa?" | SNedata$Type=="SNII" |
#                SNedata$Type=="zSNIa" |
#                SNedata$Type=="zSNII"),]

SNedata2<-SNedata[which(SNedata$Type=="zSNIa"|SNedata$Type=="zSNII"),]
# Now collapse all to SNeIa or II

library(plyr)
SNedata2$type_bin<-droplevels(revalue(SNedata2$Type, c("pSNIa"="SNIa", "SNIa?"="SNIa", "zSNIa"="SNIa",
                         "pSNII"="SNII","zSNII"="SNII")))
levels(SNedata2$type_bin)

# First data is almost ready, now let's remove missing data for simplicity. 
#But we will include a treatment for this in the final model. 
complete.cases(SNedata2)
SNedata3<-na.omit(SNedata2)
# Start the logit model

# Define data for JAGS
#X<-model.matrix(~ M_gal + log_sSFR + g.r+
#                r.i+i.z, data = SNedata3)

X<-model.matrix(~ g.r+r.i+i.z, data = SNedata3)
  
K<-ncol(X)

typeSne<-as.numeric(SNedata3$type_bin)-1
jags.data <- list(Y= typeSne,
                 N = nrow(SNedata3),
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