#Set working directory and load libraries
setwd("~/Documents/DogProject_Jaz/ROH")
library(ggplot2)
library(dplyr)
library(data.table)
library(cowplot)

#Read Files in
dfMerge = read.delim("~/Documents/DogProject_Jaz/ROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt")
#dfMerge = read.delim("~/Documents/DogProject_Jaz/ROH/Stronen2013SitesOnlyTrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt")
dfStronen = read.delim("~/Documents/DogProject_Jaz/ROH/Stronen2013_WolvesTrueROH_propCoveredwithin1SDMean_allChroms.txt")
popmapMerge = read.delim("~/Documents/DogProject_Jaz/BreedCladeInfo/BreedAndCladeInfo_mergedFitakCornell.txt")
popmapStronen = read.delim("~/Documents/DogProject_Jaz/BreedCladeInfo/Stronen2013_SamplesUsed.txt")
orderPops = read.table("~/Documents/DogProject_Jaz/BreedCladeInfo/OrderPops.txt")
#orderCluster = read.table("~/Documents/DogProject_Jaz/BreedCladeInfo/OrderCluster.txt")

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


#Make ROH score for each clade and each population 
CountsPerBreed$NormConstant = (choose(2*CountsPerBreed$n, 2)) - CountsPerBreed$n #From Nakatsuka et al. 2017 Nature Genetics

ROHScoreDF = aggregateROH %>% group_by(Population) %>% summarise(PopROHScore = sum(as.numeric(totalLen)), MeanROHperIndivCount = mean(as.numeric(countROH))) %>% mutate(sampSize = CountsPerBreed$n[match(Population, CountsPerBreed$Population)], NormPopScore = as.numeric(PopROHScore)/as.numeric(CountsPerBreed$NormConstant[match(Population, CountsPerBreed$Population)])) %>% filter(sampSize > 1) %>% as.data.frame()
ROHScoreDF$RelToEURO = ROHScoreDF$NormPopScore/ROHScoreDF[grep("EURO", ROHScoreDF$Population),]$NormPopScore
ROHScoreDF$RelToLab = ROHScoreDF$NormPopScore/ROHScoreDF[grep("labrador_retriever", ROHScoreDF$Population),]$NormPopScore

#Save final table
#write.table(ROHScoreDF, "ROHScoresPerPopulation.txt",sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#Set Populations and Clusters as factor so Wolves and dogs group together
ROHScoreDF$Population = factor(ROHScoreDF$Population, levels=orderPops$V1)

#Plot Data
mid = mean(ROHScoreDF$RelToLab)
PlotROHRel2Lab = ggplot(ROHScoreDF, aes(y=RelToLab, x=Population, fill=RelToLab)) + geom_bar(stat="identity") + coord_flip() + theme_bw() + labs(y = "Within Population ROH Score Relative to Labrador Retriever (Normalized by Sample Size)", x="Breed") + geom_hline(yintercept = mid, linetype="dashed", colour = "black") + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 10), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=20)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=18)) + scale_y_continuous(expand = c(0,0)) + scale_fill_gradient(name= "Fold-Enrichment", low="blue" ,high="red")

PlotROHRel2EuroWolf = ggplot(ROHScoreDF, aes(y=RelToEURO, x=Population)) + geom_bar(stat="identity") + coord_flip() + theme_bw() + labs(y = "Within Population ROH Score Relative to European Gray Wolf (Normalized by Sample Size)", x="Breed") + geom_hline(yintercept = 1, linetype="dashed", colour = "blue") + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 10), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=20)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=18)) + scale_y_continuous(expand = c(0,0))

