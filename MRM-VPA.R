library(ecodist)
library(vegan)

#1. read environmental variables data and calculate distance matrix
env=read.csv("env.csv",header = T,row.names = 1)
pHdis<-vegdist(scale(env$pH),method="euclidean")
SWCdis<-vegdist(scale(env$SWC),method="euclidean")
SOCdis<-vegdist(scale(env$SOC.g.kg.1),method="euclidean")
ANdis<-vegdist(scale(env$AN.mg.kg.1),method="euclidean")
MBCdis<-vegdist(scale(env$MBC.mg.kg.1),method="euclidean")
MATdis<-vegdist(scale(env$MAT..),method="euclidean")

#2. read host community data (resample_motu) and calculate distance matrix
host=read.csv("resampled_motu.csv",header = T, row.names = 1)
host=as.data.frame(t(host))
hostdis=vegdist(host,method="bray")

#3. read viral community data (resample_votu) and calculate distance matrix
votu=read.csv("resampled_votu.csv",header = T,row.names = 1)
votu=as.data.frame(t(votu))
votudis=vegdist(votu,method="bray")

#4. read geographic distance data and calculate distance matrix
geo=read.csv("geodis.csv",header = T,row.names = 1)
geodis=vegdist(geo,method="euclidean")
min.dist=15
lngeodis=log(geodis+min.dist) #log transform


#5.MRM-VPA analysis
MRM(as.dist(votudis) ~ pHdis + SWCdis +SOCdis + ANdis +MBCdis + MATdis) #env
MRM(as.dist(votudis) ~ as.dist(hostdis)) #host
MRM(as.dist(votudis) ~ as.dist(lngeodis)) #goe

MRM(as.dist(votudis) ~ as.dist(hostdis)+as.dist(lngeodis)+pHdis + SWCdis +SOCdis + ANdis +MBCdis + MATdis) #total

MRM(as.dist(votudis) ~ as.dist(hostdis)+as.dist(lngeodis)) #host+geo
MRM(as.dist(votudis) ~ as.dist(hostdis)+pHdis + SWCdis +SOCdis + ANdis +MBCdis + MATdis) #host+env
MRM(as.dist(votudis) ~ as.dist(lngeodis)+pHdis + SWCdis +SOCdis + ANdis +MBCdis + MATdis) #geo+env
