#0. read resampled votu dataset
otu <- read.csv("resampled_votu.csv",row.names = 1,header=T)
otu <- as.data.frame(t(otu))

#1. calculate the distance matrix
library(vegan)
library(stats)
data_dist <- vegdist(otu,method="bray")

#2. PCoA analysis
data_pcoa <- cmdscale(data_dist,k=(nrow(otu)-1),eig=TRUE) 
sample_site <- data.frame(data_pcoa$points)[1:2] 
sample_site$treatment <- rep(c("control","control","control","control","warming","warming","warming","warming"),4)
sample_site$depth <- rep(c("0-10","10-30","30-60","60-100"),8)
names(sample_site)[1:4] <- c('PCoA1', 'PCoA2','treatment',"depth")

#3. PERMANOVA
result <- adonis2(formula = otu ~ (sample_site$treatment * sample_site$depth),data =sample_site, permutations = 999, method = "bray") 
