##0. read resampled votu dataset
otu <- read.csv("resampled_votu.csv",row.names = 1,header=T)
otu <- as.data.frame(t(otu))

##1. calculate alpha diversity
library(vegan)
shannon <- diversity(otu,"shannon")
simpson <- diversity(otu,"simpson")
SR <- specnumber(otu,MARGIN = 1)
Pielou <- shannon/log(SR) 
diversity <- data.frame(shannon,simpson,Pielou) 
diversity <- scale(diversity())

#note "0" represents "control" treatment and "1" represent "warming" treatment."0.05,0.2,0.45 and 0.8" represent different soil depth.
diversity$treatment <- rep(c(rep("0",4),rep("1",4)),4)
diversity$depth <- rep(c(0.05,0.20,0.45,0.80),8)
diversity$block <- c(rep(1,8),rep(2,8),rep(3,8),rep(4,8))

##2. treatment effects by linear mixed model
library(lme4)
library(lmerTest)

m1 <- lmer(shannon ~ treatment*depth+(1 |block),data = diversity)
summary(m1)

m2 <- lmer(simpson ~ treatment*depth+(1 |block),data = diversity)
summary(m2)

m3 <- lmer(Pielou ~ treatment*depth+(1 |block),data = diversity)
summary(m3)