#####Load Libraries
library(tidyverse)

######Read Files in
OMIA = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/OMIA/processedCausalVarsOMIA.txt", stringsAsFactors = F) 
OMIA_nonFitness = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/OMIA/processedCausalVarsOMIA_nonFitnessRelated.txt", stringsAsFactors = F) %>%
  filter(Breed != "numerous_breeds") 
DisPrev = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/OMIA/AKC_DiseasePrev_PointEstWiles2017.txt", stringsAsFactors = F)
#LifeSpanData = read.delim("~/Documents/Documents/DogProject_Jaz/LocalRscripts/AKC/Adams2010_BreedLifeSpan_addBreeds.txt", stringsAsFactors = F)
AKC = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/AKC/AKC_breedPopularity_1926thru2005.txt", check.names = F, stringsAsFactors = F)
IBDScores = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/IBDSegs/IBDScoresPerPopulation.txt", stringsAsFactors = F) 
ROHScores = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/ROH/ROHScoresPerPopulation.txt", stringsAsFactors = F) 
popmapMerge = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/BreedAndCladeInfo_mergedFitakCornell.txt", stringsAsFactors = F) 
orderPops = read.table("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/OrderPops.txt", stringsAsFactors = F) %>%
  mutate(V1 = gsub("(?<=^|_)([a-z])", "\\U\\1", V1, perl=TRUE),
         V1 = gsub("_Dog", "_dog", V1, perl = TRUE))
#orderCluster = read.table("~/Documents/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/OrderCluster.txt", stringsAsFactors = F)

###Create Data Frames to Process
###Set Populations and Clusters as factor so Wolves and dogs group together
BreedsWithCausalVars = OMIA %>%
  count(Breed) #count causal vars associated with each breed

BreedsWithCausalVars_nonFitness = OMIA_nonFitness %>%
  count(Breed) #count causal vars associated with non-fitness related traits

#IBD dataframe
FinalIBDScores = IBDScores %>%
  mutate(PopIBDScore = PopIBDScore/10^6, NormPopScore = NormPopScore/10^6,
         CausalVars = BreedsWithCausalVars$n[match(Breed1, BreedsWithCausalVars$Breed)],
         CausalVars_nonFitness = BreedsWithCausalVars_nonFitness$n[match(Breed1, BreedsWithCausalVars_nonFitness$Breed)],
         OverallPopularityRank = AKC$popularity[match(Breed1, AKC$breed)]) %>%
  rename(Population=Breed1) %>%
  mutate(CausalVars = replace_na(CausalVars, 0),
         CausalVars_nonFitness = replace_na(CausalVars_nonFitness, 0),
         Population = factor(Population, levels=orderPops$V1))

#ROH dataframe
FinalROHScores = ROHScores %>%
  mutate(PopROHScore = PopROHScore/10^6,
         NormPopScore = NormPopScore/10^6,
         CausalVars = BreedsWithCausalVars$n[match(Population, BreedsWithCausalVars$Breed)],
         CausalVars_nonFitness = BreedsWithCausalVars_nonFitness$n[match(Population, BreedsWithCausalVars_nonFitness$Breed)],
         OverallPopularityRank = AKC$popularity[match(Population, AKC$breed)]) %>%
  mutate(CausalVars = replace_na(CausalVars, 0),
         CausalVars_nonFitness = replace_na(CausalVars_nonFitness, 0),
         Population = factor(Population, levels=orderPops$V1))

#Combine and only keep with both an ROH and IBD Score
comboDF = merge(FinalROHScores, FinalIBDScores, by ="Population") %>%
  select("Population", "PopROHScore", "NormPopScore.x", "PopIBDScore","NormPopScore.y") %>%
  mutate(Clade = popmapMerge$clade[match(Population, popmapMerge$breed)],
         CausalVars = BreedsWithCausalVars$n[match(Population, BreedsWithCausalVars$Breed)],
         CausalVars_nonFitness = BreedsWithCausalVars_nonFitness$n[match(Population, BreedsWithCausalVars_nonFitness$Breed)]) %>%
  rename(NormPopScore_ROH = NormPopScore.x, NormPopScore_IBD = NormPopScore.y) %>%
  mutate(CausalVars = replace_na(CausalVars, 0),
         CausalVars_nonFitness = replace_na(CausalVars_nonFitness, 0),
         Population = gsub("(?<=^|_)([a-z])", "\\U\\1", Population, perl=TRUE),
         Population = gsub("_Dog", "_dog", Population, perl = TRUE),
         Population = factor(Population, levels=orderPops$V1)) 


#Popularity data frame
PopularityDF = popmapMerge %>%
  select("breed", "clade") %>%
  group_by(breed) %>%
  sample_n(1) %>%
  mutate(OverallPopularityRank = AKC$popularity[match(breed, AKC$breed)],
         CausalVars = BreedsWithCausalVars$n[match(breed, BreedsWithCausalVars$Breed)],
         CausalVars_nonFitness = BreedsWithCausalVars_nonFitness$n[match(breed, BreedsWithCausalVars_nonFitness$Breed)]) %>%
  filter(!is.na(OverallPopularityRank)) %>%
  ungroup() %>%
  mutate(CausalVars = replace_na(CausalVars, 0),
         CausalVars_nonFitness = replace_na(CausalVars_nonFitness, 0),
         breed = gsub("(?<=^|_)([a-z])", "\\U\\1", breed, perl=TRUE),
         breed = gsub("_Dog", "_dog", breed, perl = TRUE)) %>%
  rename(Population=breed)
  

comboDF_noWolves = comboDF %>%  
  filter(!grepl("grayWolf",Clade)) %>% #remove labelled wolves
  na.omit() %>% #remove Stronen Wolves (no clade)
  mutate(OverallPopularityRank = PopularityDF$OverallPopularityRank[match(Population, PopularityDF$Population)])

#Add ROH and IBD Scores to the Disease Prevalance data
TDisPrev = t(DisPrev[,2:ncol(DisPrev)]) %>%
  as.data.frame() #Transpose data
colnames(TDisPrev) = DisPrev[,1] #Set the column headings from the first column in the original table
DiseasePrev = TDisPrev %>%
  rownames_to_column(var = "Population") %>%
  mutate(NormPopROHScore = comboDF_noWolves$NormPopScore_ROH[match(Population, comboDF_noWolves$Population)],
         NormPopIBDScore = comboDF_noWolves$NormPopScore_IBD[match(Population, comboDF_noWolves$Population)]) %>%
  na.omit()
rm(TDisPrev)

#####Correlations
corrROHvsIBD = lm(PopROHScore~PopIBDScore, data = comboDF)
corrROHScorevsIBDScore = lm(NormPopScore_ROH~NormPopScore_IBD, data = comboDF)

corrROHScorecausVars = lm(CausalVars~NormPopScore_ROH, data = comboDF_noWolves)
corrIBDScorecausVars = lm(CausalVars~NormPopScore_IBD, data = comboDF_noWolves)
corrROHScorecausVars_nonFitness = lm(CausalVars_nonFitness~NormPopScore_ROH, data = comboDF_noWolves)
corrIBDScorecausVars_nonFitness = lm(CausalVars_nonFitness~NormPopScore_IBD, data = comboDF_noWolves)

corrPopularitycausVars = lm(CausalVars~OverallPopularityRank, data = PopularityDF)
corrPopularitycausVars_nonFitness = lm(CausalVars_nonFitness~OverallPopularityRank, data = PopularityDF)
corrPopularityROHScore = lm(NormPopScore_ROH~OverallPopularityRank, data = comboDF_noWolves)
corrPopularityIBDScore = lm(NormPopScore_IBD~OverallPopularityRank, data = comboDF_noWolves)

