#####Load Libraries
library(dplyr)
library(tidyr)

######Read Files in
DisPrev = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/OMIA/AKC_DiseasePrev_PointEstWiles2017.txt")
IBDScores = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/IBDSegs/IBDScoresPerPopulation.txt")
ROHScores = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/ROH/ROHScoresPerPopulation.txt")
popmapMerge = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/BreedAndCladeInfo_mergedFitakCornell.txt")
orderPops = read.table("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/OrderPops.txt")


######Plot Linear Regression Function###
ggplotRegression = function (fit) {
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(size = 2) + 
    stat_smooth( method = 'lm', col = "blue") +  
    theme_bw() + 
    labs(title = bquote(R^2== ~.(signif(summary(fit)$adj.r.squared, 5))~"&"~"p"==~.(signif(summary(fit)$coef[2,4], 5))))
  
}

###Create Data Frames to Process
###Set Populations and Clusters as factor so Wolves and dogs group together
BreedsWithCausalVars = OMIA %>%
  count(Breed) #count causal vars associated with each breed

#IBD dataframe
FinalIBDScores = IBDScores %>%
  mutate(PopIBDScore = PopIBDScore/10^6, NormPopScore = NormPopScore/10^6) %>%
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

#plotting 
plottingCols = colnames(DiseasePrev)[2:10]
pdf("test.pdf", height = 8, width = 8)
for (i in plottingCols) {
  print(i)
  #linearMod = lm(i~NormPopROHScore, data = DiseasePrev)
  #print(ggplotRegression(linearMod) + 
  #        xlab("ROH Score in Mb (Normalized)"))
}
dev.off()
