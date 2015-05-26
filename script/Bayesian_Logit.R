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

data.0= read.csv(file="..//data//sdsssn_master.dat2.txt",header=TRUE,dec=".",sep="")


# Select few variables for test (avoiding double observations for now). Using only FSPS
# Let's ignore the errors for now

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

data.1=SNe_data[,c("logMassFSPS","ageFSPS","logSSFRFSPS","dereduhost","deredghost",
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

SNedata2<-SNedata[which(SNedata$Type=="pSNIa" | 
                SNedata$Type=="pSNII" | SNedata$Type=="SNIa"
              | SNedata$Type=="SNIa?" | SNedata$Type=="SNII" |
                SNedata$Type=="zSNIa" |
                SNedata$Type=="zSNII"),]
# Now collapse all to SNeIa or II

library(plyr)
SNedata2$type_bin<-droplevels(revalue(SNedata2$Type, c("pSNIa"="SNIa", "SNIa?"="SNIa", "zSNIa"="SNIa",
                         "pSNII"="SNII","zSNII"="SNII")))
levels(SNedata2$type_bin)

# First data is almost ready, now let's remove missing data for simplicity. 
#But we will include a treatment for this in the final model. 


                      