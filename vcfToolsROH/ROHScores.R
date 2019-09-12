#Set working directory and load libraries
setwd("~/DogProject_Jaz/ROH")
library(ggplot2)
library(reshape2)
library(dplyr)
library(randomcoloR)
library(data.table)
library(ggpubr)
library(cowplot)

#Read Files in
#dfMerge = read.delim("~/DogProject_Jaz/ROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt")
dfMerge = read.delim("~/DogProject_Jaz/ROH/Stronen2013SitesOnlyTrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt")
dfStronen = read.delim("~/DogProject_Jaz/ROH/Stronen2013_WolvesTrueROH_propCoveredwithin1SDMean_allChroms.txt")
popmapMerge = read.delim("~/DogProject_Jaz/BreedCladeInfo/BreedAndCladeInfo_mergedFitakCornell.txt")
popmapStronen = read.delim("~/DogProject_Jaz/BreedCladeInfo/Stronen2013_SamplesUsed.txt")
orderPops = read.table("~/DogProject_Jaz/BreedCladeInfo/OrderPops.txt")
orderCluster = read.table("~/DogProject_Jaz/BreedCladeInfo/OrderCluster.txt")

#Concatenate dataframes and aggregate ROH
WolfDog = rbind.data.frame(dfMerge,dfStronen)
aggregateROH = WolfDog %>% select(AUTO_LEN, INDV) %>% group_by(INDV) %>% summarise(totalLen = sum(AUTO_LEN)) %>% mutate(totalLenGb = totalLen/10^9) %>% as.data.frame()
CountROHperIndiv = WolfDog %>% select(AUTO_LEN, INDV) %>% group_by(INDV) %>% tally()
aggregateROH$countROH = CountROHperIndiv$n[match(aggregateROH$INDV, CountROHperIndiv$INDV)]

#Assign breed,clade,cluster etc. and bind dog and wolf together
shortPopmapWolves = popmapStronen %>% select(IID_Col2,Population,Cluster)
#shortPopmapWolves$Population = shortPopmapWolves$Cluster #make the cluster the population for wolves
names(shortPopmapWolves)[1] = "dogID"
names(shortPopmapWolves)[2] = "breed"
names(shortPopmapWolves)[3] = "clade"
mergedPopmap = rbind.data.frame(popmapMerge, shortPopmapWolves)
mergedPopmap$breed = gsub("large_munsterlander","munsterlander_large", mergedPopmap$breed)

##make sure names are unique (they are)
#length(unique(mergedPopmap$dogID))
#length(mergedPopmap$dogID)

#Add Population and Cluster information
aggregateROH$Population = mergedPopmap$breed[match(aggregateROH$INDV,mergedPopmap$dogID)]
aggregateROH$Cluster = mergedPopmap$clade[match(aggregateROH$INDV,mergedPopmap$dogID)]
CountsPerBreed = aggregateROH %>% group_by(Population) %>% tally()
CountsPerCluster = aggregateROH %>% group_by(Cluster) %>% tally()
Breeds_grThaneql30 = CountsPerBreed %>% filter(n >= 30) %>% as.data.frame()

#Make ROH score for each clade and population with at least 30 individuals
CountsPerBreed$NormConstant = (choose(2*CountsPerBreed$n, 2)) - CountsPerBreed$n #From Nakatsuka et al. 2017 Nature Genetics

ROHScoreDF = aggregateROH %>% group_by(Population) %>% summarise(PopScore = sum(as.numeric(totalLen)), MeanCountROH = mean(as.numeric(countROH))) %>% mutate(sampSize = CountsPerBreed$n[match(Population, CountsPerBreed$Population)], NormPopScore = as.numeric(PopScore)/as.numeric(CountsPerBreed$NormConstant[match(Population, CountsPerBreed$Population)])) %>% filter(sampSize > 1) %>% as.data.frame()


