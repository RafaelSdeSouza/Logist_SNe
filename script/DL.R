# Constants
pi<-3.1415926


# Physical constants (cgs)
c <- 3*1e5#speed of light km/s

pc<-3.08*1e18# 1 pc in cm 


#Cosmological parameters
H0 <-71 #Hubble constant
h.0<-0.71#dimensionless Hubble constant
H0.t<-7.37205e-13
Omega.m <-0.27#Omega matter
Omega.l <-0.73#Omega lambda
Omega.b<-0.044

H.z.inv<-function(z){1/(H0*(Omega.m*(1+z)**3+Omega.l)**0.5)}#1/H(z)
dL<-function(zmax){1e6*3e5*(1+zmax)*integrate(H.z.inv,0,zmax)$value}#Luminosity Distance in pc

AbMag<-function(mag,z){mag-5*(log(dL(z),10)-1)}