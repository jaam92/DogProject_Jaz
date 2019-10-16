#####Load Libraries
library(dplyr)
library(tidyr)

######Read Files in
OMIA = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/OMIA/processedCausalVarsOMIA.txt")
#LifeSpanData = read.delim("~/Documents/Documents/DogProject_Jaz/LocalRscripts/AKC/Adams2010_BreedLifeSpan_addBreeds.txt")
AKC = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/AKC/AKC_breedPopularity_1926thru2005.txt", check.names = F)
IBDScores = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/IBDSegs/IBDScoresPerPopulation.txt")
ROHScores = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/ROH/ROHScoresPerPopulation.txt")
popmapMerge = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/BreedAndCladeInfo_mergedFitakCornell.txt")
orderPops = read.table("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/OrderPops.txt")
#orderCluster = read.table("~/Documents/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/OrderCluster.txt")

###Create Data Frames to Process
###Set Populations and Clusters as factor so Wolves and dogs group together
BreedsWithCausalVars = OMIA %>%
  count(Breed) #count causal vars associated with each breed

#IBD dataframe
FinalIBDScores = IBDScores %>%
  mutate(PopIBDScore = PopIBDScore/10^6, NormPopScore = NormPopScore/10^6,
         CausalVars = BreedsWithCausalVars$n[match(Breed1, BreedsWithCausalVars$Breed)],
         OverallPopularityRank = AKC$popularity[match(Breed1, AKC$breed)]) %>%
  rename(Population=Breed1) %>%
  mutate(CausalVars = replace_na(CausalVars, 0),
         Population = factor(Population, levels=orderPops$V1))

#ROH dataframe
FinalROHScores = ROHScores %>%
  mutate(PopROHScore = PopROHScore/10^6,
         NormPopScore = NormPopScore/10^6,
         CausalVars = BreedsWithCausalVars$n[match(Population, BreedsWithCausalVars$Breed)],
         OverallPopularityRank = AKC$popularity[match(Population, AKC$breed)]) %>%
  mutate(CausalVars = replace_na(CausalVars, 0),
         Population = factor(Population, levels=orderPops$V1))

#Combine and only keep with both an ROH and IBD Score
comboDF = merge(FinalROHScores, FinalIBDScores, by ="Population") %>%
  select("Population", "PopROHScore", "NormPopScore.x", "PopIBDScore","NormPopScore.y") %>%
  mutate(Clade = popmapMerge$clade[match(Population, popmapMerge$breed)],
         CausalVars = BreedsWithCausalVars$n[match(Population, BreedsWithCausalVars$Breed)]) %>%
  rename(NormPopScore_ROH = NormPopScore.x, NormPopScore_IBD = NormPopScore.y) %>%
  mutate(CausalVars = replace_na(CausalVars, 0),
         Population = factor(Population, levels=orderPops$V1))

#Popularity data frame
PopularityDF = popmapMerge %>%
  select("breed", "clade") %>%
  group_by(breed) %>%
  sample_n(1) %>%
  mutate(OverallPopularityRank = AKC$popularity[match(breed, AKC$breed)],
         CausalVars = BreedsWithCausalVars$n[match(breed, BreedsWithCausalVars$Breed)]) %>%
  filter(!is.na(OverallPopularityRank)) %>%
  mutate(CausalVars = replace_na(CausalVars, 0)) %>%
  rename(Population=breed)

#####Correlations
corrROHvsIBD = lm(PopROHScore~PopIBDScore, data = comboDF)
corrROHScorevsIBDScore = lm(NormPopScore_ROH~NormPopScore_IBD, data = comboDF)

corrROHScorecausVars = lm(CausalVars~NormPopScore, data = FinalROHScores)
corrIBDScorecausVars = lm(CausalVars~NormPopScore, data = FinalIBDScores)

corrROHCountcausVars = lm(CausalVars~MeanROHperIndivCount, data = FinalROHScores)
corrIBDCountcausVars = lm(CausalVars~MeanIBDperIndivCount, data = FinalIBDScores)

corrPopularitycausVars = lm(CausalVars~OverallPopularityRank, data = PopularityDF)
corrPopularityROHScore = lm(NormPopScore~OverallPopularityRank, data = FinalROHScores)
corrPopularityIBDScore = lm(NormPopScore~OverallPopularityRank, data = FinalIBDScores)

#Does correcting for popularity improve associations with ROH and IBD
multiVarcorrIBDScorecausVars = lm(CausalVars~NormPopScore + OverallPopularityRank, data = FinalIBDScores)
multiVarcorrIBDScorecausVars_plusInt = lm(CausalVars~NormPopScore + OverallPopularityRank + OverallPopularityRank*NormPopScore, data = FinalIBDScores)

multiVarcorrROHScorecausVars = lm(CausalVars~NormPopScore + OverallPopularityRank, data = FinalROHScores)