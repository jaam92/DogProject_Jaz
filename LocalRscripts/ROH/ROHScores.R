#Set working directory and load libraries
setwd("~/Documents/DogProject_Jaz/LocalRscripts/ROH")
library(tidyverse)

#Read Files in
dfMerge = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/ROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt")
#dfMerge = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/ROH/Stronen2013SitesOnlyTrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt")
dfStronen = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/ROH/Stronen2013_WolvesTrueROH_propCoveredwithin1SDMean_allChroms.txt") 
popmapMerge = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/BreedAndCladeInfo_mergedFitakCornell.txt")
popmapStronen = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/Stronen2013_SamplesUsed.txt")
orderPops = read.table("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/OrderPops.txt")
#orderCluster = read.table("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/OrderCluster.txt")

#Concatenate dataframes 
WolfDog = rbind.data.frame(dfMerge,dfStronen)

shortPopmapWolves = popmapStronen %>% 
  select(IID_Col2,Population,Cluster) %>%
  rename("dogID" = IID_Col2, "breed" = Population, "clade" = Cluster)

#Merge everything
mergedPopmap = rbind.data.frame(popmapMerge, shortPopmapWolves) %>%
  mutate(breed = gsub("large_munsterlander","munsterlander_large", breed)) 

##make sure names are unique (they are)
#length(unique(mergedPopmap$dogID))
#length(mergedPopmap$dogID)

#Count ROH
CountROHperIndiv = WolfDog %>% 
  select(AUTO_LEN, INDV) %>% 
  group_by(INDV) %>% 
  tally()

#Aggregate ROH and add Population and Cluster information
aggregateROH = WolfDog %>% 
  select(AUTO_LEN, INDV) %>% 
  group_by(INDV) %>% 
  summarise(totalLen = sum(AUTO_LEN)) %>% 
  mutate(totalLenGb = totalLen/10^9,
         countROH = CountROHperIndiv$n[match(INDV, CountROHperIndiv$INDV)],
         Population = mergedPopmap$breed[match(INDV,mergedPopmap$dogID)],
         Cluster = mergedPopmap$clade[match(INDV,mergedPopmap$dogID)])

#Count of ROH per breed
CountsPerBreed = aggregateROH %>% 
  group_by(Population) %>% 
  tally() %>%
  ungroup() %>%
  mutate(NormConstant = (choose(2*n, 2)) - n) #Make ROH score for each clade and each population from Nakatsuka et al. 2017 Nature Genetics


#Compute ROH Score
ROHScoreDF = aggregateROH %>% 
  group_by(Population) %>% 
  summarise(PopROHScore = sum(as.numeric(totalLen)), 
            MeanROHperIndivCount = mean(as.numeric(countROH))) %>% 
  mutate(sampSize = CountsPerBreed$n[match(Population, CountsPerBreed$Population)], 
         NormPopScore = as.numeric(PopROHScore)/as.numeric(CountsPerBreed$NormConstant[match(Population, CountsPerBreed$Population)])) %>% 
  filter(sampSize > 1)

#Make scores reltive to European wolves and labs
ROHScoreDF$RelToEURO = ROHScoreDF$NormPopScore/ROHScoreDF[grep("EURO", ROHScoreDF$Population),]$NormPopScore
ROHScoreDF$RelToLab = ROHScoreDF$NormPopScore/ROHScoreDF[grep("labrador_retriever", ROHScoreDF$Population),]$NormPopScore

#Save final table
#write.table(ROHScoreDF, "ROHScoresPerPopulation.txt",sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#Plot Data
ROHScoreDF$Population = factor(ROHScoreDF$Population, levels=orderPops$V1)#Set Populations and Clusters as factor so Wolves and dogs group together
mid = mean(ROHScoreDF$RelToLab)
PlotROHRel2Lab = ggplot(ROHScoreDF, aes(y=RelToLab, x=Population, fill=RelToLab)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() + 
  labs(y = "Within Population ROH Score Relative to Labrador Retriever (Normalized by Sample Size)", x="Breed") + 
  geom_hline(yintercept = mid, linetype="dashed", colour = "black") + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 10), 
        plot.title=element_text(size=26, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_gradient(name= "Fold-Enrichment", low="blue" ,high="red")


