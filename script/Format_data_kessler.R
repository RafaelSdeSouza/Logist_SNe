# Script to prepare data set for logit analysis
require(MASS)
require(gdata)
require(plyr)
#Read dataset

data.0= read.table(file="..//data/match",header=TRUE,na.strings = "",sep="")
write.matrix(data.0,"..//data/sako.dat",sep = "\t")

# Select few variables for test. For now, galaxy morphology and SN type

SN_cat<-data.frame(SNtype=data.0$SNtype,Galtype=data.0$Galtype2,logSSFRF=data.0$logSSFRFSPS)

SN_cat2<-na.omit(SN_cat)
require(gdata)
SN_cat2$SNtype<-trim(SN_cat2$SNtype)

SN_cat3<-SN_cat2[which(SN_cat2$SNtype=="Ia"|SN_cat2$SNtype=="II"|SN_cat2$SNtype=="IIn"|
                        SN_cat2$SNtype=="Ib"|SN_cat2$SNtype=="Ib/c"|SN_cat2$SNtype=="Ic"),]
SN_cat3$SNtype<-droplevels(SN_cat3$SNtype)
SN_cat3$SNtype<-revalue(SN_cat2$SNtype,c("Ia"="Ia","Ib"="CC","Ib/c"="CC","Ic"="CC","II"="CC"))

SN_cat3<-na.omit(SN_cat3)
# Start the logit model

# Define data for JAGS
#X<-model.matrix(~ M_gal + log_sSFR + g.r+
#                r.i+i.z, data = SNedata3)


SN_cat3$Galtype<-trim(SN_cat3$Galtype)

SN_cat3$Galtype2<-revalue(SN_cat3$Galtype, c("E pec" ="E", "E/S0 pec"="E/S0","S pec"="S","S0 pec"="S0",
                                             "S0/a"="S0","S0/a pec"="S0","Sa"="S","Sa pec"="S","Sab"="S","Sab pec"="S","Sb"="S","Sb pec"="S",
                                             "Sbc"="S","Sbc pec"="S","Sc"="S","Sc pec"="S","Scd"="S",
                                             "Scd pec"="S","Sd"="S",
                                             "Sd pec"="S","Sdm"="S","Sdm pec"="S","Sm"="S"))
SN_cat4<-SN_cat3[which(SN_cat3$Galtype2=="E"|SN_cat3$Galtype2=="E/S0"|SN_cat3$Galtype2=="Im"|
                         SN_cat3$Galtype2=="S"|SN_cat3$Galtype2=="S0"),]
SN_cat4$Galtype2<-droplevels(SN_cat4$Galtype2)

write.matrix(SN_cat4[,-2],"..//data/multinomial_cat.dat",sep = "\t")


write.matrix(SN_cat3,"..//data/logit_combined.dat",sep = "\t")