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
jags.data <- list(y= typeSne,
                  N = nrow(data.1),
                 galtype = galtype,
                  bar=bar,
                  J=6,
                  Ntype = Ntype
                  )


model<- "model{
  for(i in 1:N){
    y[i] ~ dcat(p[i, 1:J])
    
    for (j in 1:J){
      log(q[i,j]) <-  beta[1,j] + beta[2,j]*bar[i]+ranef[galtype[i],j]
      
      p[i,j] <- q[i,j]/sum(q[i,1:J])
    }   # close J loop
    
  }  # close N loop

# Priors

#tau.R<-pow(sdBeta,-1)
#sdBeta ~ dgamma(0.001,0.001)

# Random intercept model
for (k in 1:Ntype){
#ranef[j]~dnorm(0,1e-5)
ranef[k,1]~dbern(0)
for(j in 2:J){
ranef[k,j]~ddexp(0,1/1.64)
}
}

  for(k in 1:2){
   beta[k,1]~dbern(0) ## MUST set the first set of covariates (for the first outcome category) to 0
    for(j in 2:J){
      beta[k,j] ~ dnorm(0, 0.1)
    }  # close J loop
  }  # close K loop
}"  # close model loop 




params<-c("beta","p","ranef")

inits<-function(){list(beta=structure(.Data=c(rep(NA,2),runif(2*(6-1),-1,1)),.Dim=c(2,6)))}


inits1=inits()
inits2=inits()
inits3=inits()
chain_inits<-list(inits1,inits2,inits3)




library(parallel)
cl <- makeCluster(3)
jags.mlogit <- run.jags(method="rjparallel", 
                       data = jags.data, 
                      inits = chain_inits,
                       model=model,
                       n.chains = 3,
                       adapt=1000,
                       monitor=c(params),
                       burnin=2000,
                       sample=5000,
                       summarise=FALSE,
                       plots=FALSE
)
summary<-extend.jags(jags.mlogit,drop.monitor=c("p"), summarise=TRUE)
print(summary)

jagssamples <- as.mcmc.list(jags.mlogit)
L.factors <- data.frame(
  Parameter=c(paste("beta[1,", seq(1:6), "]", sep=""),paste("beta[2,", seq(1:6), "]", sep="")),
  Label=c("Ia","Ib","Ib/c","Ic","II",
        "IIn","Bar-Ia","Bar-Ib","Bar-Ib/c","Bar-Ic","Bar-II",
        "Bar-IIn"))
#head(L.factors)

beta_post<-ggs(jagssamples,par_labels=L.factors,family=c("beta"))
pdf("..//figures/multi_beta.pdf",width=8,height=10)
ggs_caterpillar(beta_post)+theme_stata()+ylab("")
dev.off()

L.factors2<-data.frame(
  Parameter=c(paste("ranef[1,",seq(1:6), "]", sep=""),paste("ranef[2,",seq(1:6), "]", sep=""),paste("ranef[3,",seq(1:6), "]", sep=""),
              paste("ranef[4,",seq(1:6), "]", sep=""),paste("ranef[5,",seq(1:6), "]", sep="")),
  Label=c("E-Ia","E-Ib","E-Ibc","E-Ic","E-II","E-IIn","E/S0-Ia","E/S0-Ib","E/S0-Ibc","E/S0-Ic","E/S0-II","E/S0-IIn",
          "S-Ia","S-Ib","S-Ibc","S-Ic","S-II","S-IIn","S0-Ia","S0-Ib","S0-Ibc","S0-Ic","S0-II","S0-IIn",
          "Im-Ia","Im-Ib","Im-Ibc","Im-Ic","Im-II","Im-IIn"),
  supernova=as.factor(rep(c("Ia","Ib","Ibc","Ic","II","IIn"),5)),galaxies=as.factor(rep(c("E","E/S0","S","S0","Im"),each=6)))

cof_post<-ggs(jagssamples,par_labels=L.factors2,family="ranef")
pdf("..//figures/multi_ranef.pdf",width=12,height=20)
ggs_caterpillar(cof_post)+theme_stata()+ylab("")
dev.off()





# probabilities
prob<-summary(as.mcmc.list(jags.mlogit, vars="p"))
prob<-prob$quantiles


