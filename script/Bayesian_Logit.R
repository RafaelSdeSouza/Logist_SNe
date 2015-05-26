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

SNe_data= read.csv(file="..//data//sdsssn_master.dat2.txt",header=TRUE,dec=".",sep="")
dim(GCS)
