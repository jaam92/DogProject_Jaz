#Set working directory and load libraries
setwd("~/Documents/DogProject_Jaz/LocalRscripts/IBDSegs/")
library(dplyr)
library(data.table)
library(cowplot)
library(ggplot2)

#Read Files in
dfMerge = fread("~/Documents/DogProject_Jaz/LocalRscripts/IBDSegs/IBDSeq/MergedFitakCornell_allChroms_Haplotypes_IBDSeq.ibd")
dfStronen = fread("~/Documents/DogProject_Jaz/LocalRscripts/IBDSegs/Stronen2013Wolves_allChroms_Haplotypes_IBDSeq.ibd")
popmapMerge = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/BreedAndCladeInfo_mergedFitakCornell.txt")
popmapStronen = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/Stronen2013_SamplesUsed.txt")
orderPops = read.table("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/OrderPops.txt")
#orderCluster = read.table("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/OrderCluster.txt")

#Assign breed,clade,cluster etc. and bind dog and wolf together
shortPopmapWolves = popmapStronen %>% select(IID_Col2,Population,Cluster)
#shortPopmapWolves$Population = shortPopmapWolves$Cluster #make the cluster the population for wolves
names(shortPopmapWolves)[1] = "dogID"
names(shortPopmapWolves)[2] = "breed"
names(shortPopmapWolves)[3] = "clade"
mergedPopmap = rbind.data.frame(popmapMerge, shortPopmapWolves)
mergedPopmap$breed = gsub("large_munsterlander","munsterlander_large", mergedPopmap$breed)

#Merge data frames together
dogWolfDF = rbind.data.frame(dfMerge, dfStronen)

#Figure out breeds of dogs and wolves
dogWolfDF$Breed1 = mergedPopmap$breed[match(dogWolfDF$V1, mergedPopmap$dogID)]
dogWolfDF$Breed2 = mergedPopmap$breed[match(dogWolfDF$V3, mergedPopmap$dogID)]

#Length and counts of IBDSegs per indiv and number of indivs per Breed
dogWolfDF$segLen = dogWolfDF$V7 - dogWolfDF$V6
CountIndvperBreed = dogWolfDF %>% filter(Breed1 == Breed2) %>% group_by(Breed1) %>% summarise(sampSize = n_distinct(V1)) %>% mutate(normConstant = (choose(2*as.numeric(sampSize), 2)) - as.numeric(sampSize)) #constant comes from Nakatsuka et al. 2017 Nature Genetics)
MeanIBDperBreed = dogWolfDF %>% filter(Breed1 == Breed2) %>% group_by(Breed1) %>% tally() %>% mutate(meanIBDperIndiv = as.numeric(n)/as.numeric(CountIndvperBreed$sampSize[match(Breed1, CountIndvperBreed$Breed1)]))
#Generate a data frame with scores per breed
#Only using IBD segments within breed and greater than 4Mb 
#chose 4Mb rather than 3Mb because that's the cutoff used to reliably call segments from array data with IBDSeq

IBDScoreDF = dogWolfDF %>% filter(segLen >= 4e6 & Breed1 == Breed2) %>% select(segLen, V1, Breed1) %>% group_by(Breed1) %>% summarise(PopIBDScore = sum(as.numeric(segLen))) %>% mutate(MeanIBDperIndivCount = MeanIBDperBreed$meanIBDperIndiv[match(Breed1, MeanIBDperBreed$Breed1)], sampSize = CountIndvperBreed$sampSize[match(Breed1, CountIndvperBreed$Breed1)], NormPopScore = as.numeric(PopIBDScore)/as.numeric(CountIndvperBreed$normConstant[match(Breed1, CountIndvperBreed$Breed1)])) %>% filter(sampSize > 1) %>% as.data.frame()

IBDScoreDF$RelToEURO = IBDScoreDF$NormPopScore/IBDScoreDF[grep("EURO", IBDScoreDF$Breed1),]$NormPopScore
IBDScoreDF$RelToLab = IBDScoreDF$NormPopScore/IBDScoreDF[grep("labrador_retriever", IBDScoreDF$Breed1),]$NormPopScore

#Save final table
#write.table(IBDScoreDF, "IBDScoresPerPopulation.txt",sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#Set Populations and Clusters as factor so Wolves and dogs group together
IBDScoreDF$Breed1 = factor(IBDScoreDF$Breed1, levels=orderPops$V1)

#Plot IBDScores 
mid = mean(IBDScoreDF$RelToLab)
PlotRelToLab = ggplot(IBDScoreDF, aes(y=RelToLab, x=Breed1, fill=RelToLab)) + geom_bar(stat="identity") + coord_flip() + theme_bw() + labs(y = "Within Population IBD Score Relative to Labrador Retriever (Normalized by Sample Size)", x="Breed") + geom_hline(yintercept = 1, linetype="dashed", colour = "black") + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 10), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=20)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=18)) + scale_y_continuous(expand = c(0,0)) + scale_fill_gradient2(name= "Fold-Enrichment",midpoint=mid, low="blue", mid="white",high="red", space ="Lab")

PlotRelToEuroWolf = ggplot(IBDScoreDF, aes(y=RelToEURO, x=Breed1)) + geom_bar(stat="identity") + coord_flip() + theme_bw() + labs(y = "Within Population IBD Score Relative to European Gray Wolf (Normalized by Sample Size)", x="Breed") + geom_hline(yintercept = 1, linetype="dashed", colour = "blue") + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 10), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=20)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=18)) + scale_y_continuous(expand = c(0,0))

