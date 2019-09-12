library(ape)
library(igraph)
library(ggplot2)
library(reshape2)
library(dplyr)

#read file in
setwd("~/DogProject_Jaz/Fst")
df = read.delim("~/DogProject_Jaz/Fst/pairwiseWeighted.fst")
orderGrEql30 = read.table("~/DogProject_Jaz/Fst/breeds_grEql30.txt")

#make data framewith comparisons the other way from current
df2 = df %>% select(breed1 = breed2, breed2 = breed1,WeightedFst)

#merge together
mergeDF = rbind.data.frame(df,df2)

#Set breeds as factor so Wolves and dogs group together
mergeDF$breed1 = factor(mergeDF$breed1, levels=orderGrEql30$V1)
mergeDF$breed2 = factor(mergeDF$breed2, levels=orderGrEql30$V1)

#plot
ggplot(data = mergeDF, aes(x=breed1, y=breed2,fill=WeightedFst), colour = "white") +  geom_tile() + scale_fill_gradient(low = "blue", high = "red") + theme_bw() + theme(axis.text.x = element_text(size  = 24,angle=40, vjust=1, hjust=1), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24),legend.title=element_text(size=24), legend.text=element_text(size=18)) 


#Plot as cladogram
g = graph.data.frame(fst, directed=FALSE)
fstMatrix=get.adjacency(g, attr="WeightedFst", sparse=FALSE)
#write.csv(fstMatrix,file="FstMatrix.csv")
fstMatrixDistance=nj(dist(fstMatrix))
plot(fstMatrixDistance,type="cladogram",cex=1,lwd=10,col='red')