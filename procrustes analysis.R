#0. read resampled votu and motu dataset
motu <- read.csv("resampled_motu.csv",row.names = 1)
votu <- read.csv("resampled_votu.csv",row.names = 1)

motu <- as.data.frame(t(motu))
votu <- as.data.frame(t(votu))

#1. calculate the distance matrix
library(vegan)
dist_m <- vegdist(motu,method="bray")
dist_v <- vegdist(votu,method="bray")

#2. PCoA analysis
library(ape)
pcoa_m <- pcoa(dist_m) 
pcoa_v <- pcoa(dist_v) 

#3. Procrustes analysis
vare.proc <- procrustes( pcoa_m$vectors, pcoa_v$vectors)
vare.proc
summary(vare.proc)
plot(vare.proc)
plot(vare.proc, kind=2)
residuals(vare.proc)

protest<- protest(X = pcoa_m$vectors, Y = pcoa_v$vectors,  permutations = 999)
protest#p=0.001,cor=0.847


#4. related plot
mloc<-vare.proc$X
vloc<-vare.proc$Yrot

mloc=as.data.frame(mloc)
vloc=as.data.frame(vloc)

dat <- data.frame(x_viral=vloc$V1,y_viral=vloc$V2,x_micro=mloc$Axis.1, y_micro=mloc$Axis.2)
dat$treatment <- rep(c(rep("CK",4),rep("W",4)),4)
dat$depth <- rep(c(5,20,45,80),8)
dat$depth <- as.factor(dat$depth)

library(ggplot2)
ggplot(dat) +
  geom_segment(aes(x=x_viral,y=y_viral,xend=x_micro,yend=y_micro,colour=treatment),linewidth=1)+ 
  geom_point(size=4.5, colour = "gray24", aes(x=x_micro, y=y_micro, fill=treatment,shape=depth)) +
  scale_shape_manual(values=c(21,24,23,22))+
  scale_colour_manual(values=c("#2b8cbe","#e37d67"))+
  scale_fill_manual(values=c("#2b8cbe","#e37d67"))+
  xlab("Dimension 1")+
  ylab("Dimension 2")+
  theme_linedraw()+
  theme( panel.grid.minor = element_blank()) 