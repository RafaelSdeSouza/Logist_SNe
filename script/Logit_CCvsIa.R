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
require(gdata)
require(runjags)
require(gdata)
#Read the already clean dataset

data.1= read.table(file="..//data//logit_cat.dat",header=TRUE,na.strings = "",sep="\t")
data.1$Galtype2<-trim(data.1$Galtype2)

# Transform in absolute magnitude

data.1$Mag<-AbMag(data.1$mag_g,data.1$redshift)

galtype<-match(data.1$Galtype,c("E","E/S0","S","S0","Im"))
Ntype<-length(unique(data.1$Galtype))


typeSne<-match(trim(data.1$SNtype),c("Ia","CC"))-1
bar<-as.numeric(data.1$bar)-1
jags.data <- list(Y= typeSne,
                 N = nrow(data.1),
                 Mag = data.1$Mag,
                 galtype = galtype,
#                 bar=bar,
                 Ntype=Ntype
                 )


model<-"model{
#1. Priors

tau.R<-pow(sdBeta,-1)
sdBeta ~ dgamma(0.001,0.001)
#alpha1~dnorm(0,0.00001)
#tau.z~dgamma(0.001,0.001)

# Random intercept 
for (j in 1:Ntype){
ranef[j]~ddexp(0,tau.R)
#ranef[j]~dnorm(alpha1,tau.z)
}

beta.0~dnorm(0,0.001)
beta.1~dnorm(0,0.001)
#beta.2~dnorm(0,0.001)

#2. Likelihood
for (i in 1:N){
Y[i] ~ dbern(pi[i])
logit(pi[i]) <-  eta[i]
eta[i] <- beta.0+beta.1*Mag[i]+ranef[galtype[i]]
#eta[i] <- beta.0+beta.1*Mag[i]+ranef[galtype[i]]
#3. Prediction
prediction[i]~dbern(pi[i])
}
}"

params <- c("beta.0","beta.1","prediction","ranef","pi")

inits1=list(beta.0=rnorm(1,0,1),beta.1=rnorm(1,0,1))
inits2=list(beta.0=rnorm(1,0,1),beta.1=rnorm(1,0,1))
inits3=list(beta.0=rnorm(1,0,1),beta.1=rnorm(1,0,1))

library(parallel)
cl <- makeCluster(3)
jags.logit <- run.jags(method="rjparallel", 
                     data = jags.data, 
                     inits = list(inits1,inits2,inits3),
                     model=model,
                     n.chains = 3,
                     adapt=2500,
                     monitor=c(params),
                     burnin=20000,
                     sample=40000,
                     summarise=FALSE,
                     plots=FALSE
)

jagssamples <- as.mcmc.list(jags.logit)
summary<-extend.jags(jags.logit,drop.monitor=c("prediction","pi"), summarise=TRUE)
print(summary)


require(ggmcmc)
L.factors <- data.frame(
  Parameter=paste("ranef[", seq(1:5), "]", sep=""),
  Label=c("E","E/S0","S","S0","Im"))
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

require(caret)


xtab <- table(predtype[,3], typeSne)
confusionMatrix(xtab)

# probabilities
prob<-summary(as.mcmc.list(jags.logit, vars="pi"))
prob<-prob$quantiles

library(pROC)
ROCF<-data.frame(True=typeSne,Predicted=prob[,3])
F1 <-roc(ROCF$True,ROCF$Predicted)

pdf("..//figures/ROC_GLM.pdf")
plot.roc(F1, col="blue",auc.polygon.col="blue", print.auc=TRUE)
dev.off()

gd<-cbind(data.1,prob,deparse.level = 2)

#gd<-read.fwf("../data/ggdata.dat",width=c(11,13,13,14,27))
colnames(gd)<-c("SNtype","Galtype2","logSSFRF","lw2","lw1","prob","up1","up2")
gd$Galtype2<-trim(as.factor(gd$Galtype2))
#gd$bar<-trim(gd$bar)
#gd$bar<-as.numeric(gd$bar)
#gd$bar<-as.factor(gd$bar)
library(plyr)
#gd$bar<-revalue(gd$bar, c("1"="No", "2"="Yes"))
#gd$bar<-as.factor(gd$bar)
pdf("..//figures/probs_GLM2.pdf",height = 10,width = 12)
ggplot(data=gd,aes(x=logSSFRF,y=prob,colour=Galtype2))+
  geom_point(size=3)+
  theme_stata()+xlab("u-band")+ylab("Probability of CC event")
#+scale_x_reverse()+
  
  scale_color_gdocs(name="Galaxy type")+
#  scale_shape_manual(values=c(3,19),name="Bar")+
  theme(strip.background = element_rect(fill="gray95"),
  legend.position="bottom",plot.title = element_text(hjust=0.5),
                                                       axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=20),
                                                       strip.text.x=element_text(size=20),
                                                       axis.title.x=element_text(vjust=-0.25),
 text = element_text(size=20))
dev.off()



