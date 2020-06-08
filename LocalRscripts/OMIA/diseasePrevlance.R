#####Load Libraries
library(dplyr)
library(tidyr)

######Read Files in
WilesData = read.delim("~/DogProject_Jaz/LocalRscripts/OMIA/AKC_DiseasePrev_PointEstWiles2017.txt", check.names = F)
IBDScores = read.delim("~/DogProject_Jaz/LocalRscripts/IBDSegs/IBDScoresPerPopulation.txt")
ROHScores = read.delim("~/DogProject_Jaz/LocalRscripts/ROH/ROHScoresPerPopulation.txt")
AKC = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/AKC/AKC_breedPopularity_1926thru2005.txt", check.names = F)
popmapMerge = read.delim("~/DogProject_Jaz/LocalRscripts/BreedCladeInfo/BreedAndCladeInfo_mergedFitakCornell.txt")
orderPops = read.table("~/DogProject_Jaz/LocalRscripts/BreedCladeInfo/OrderPops.txt")


######Plot Linear Regression Function###
ggplotRegression = function (fit) {
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(size = 2) + 
    stat_smooth( method = 'lm', col = "blue") +  
    theme_bw() + 
    labs(title = bquote(R^2== ~.(signif(summary(fit)$adj.r.squared, 5))~"&"~"p"==~.(signif(summary(fit)$coef[2,4], 5))))
  
}

###Create Data Frames to Process
#IBD dataframe
FinalIBDScores = IBDScores %>%
  mutate(PopIBDScore = PopIBDScore/10^6, NormPopScore = NormPopScore/10^6) %>%
  rename(Population=Breed1) %>%
  mutate(Population = factor(Population, levels=orderPops$V1))

#ROH dataframe
FinalROHScores = ROHScores %>%
  mutate(PopROHScore = PopROHScore/10^6,
         NormPopScore = NormPopScore/10^6,
         OverallPopularityRank = AKC$popularity[match(Population, AKC$breed)]) %>%
  mutate(Population = factor(Population, levels=orderPops$V1))

#Combine and only keep with both an ROH and IBD Score
comboDF = merge(FinalROHScores, FinalIBDScores, by ="Population") %>%
  select("Population", "PopROHScore", "NormPopScore.x", "PopIBDScore","NormPopScore.y") %>%
  mutate(Clade = popmapMerge$clade[match(Population, popmapMerge$breed)],) %>%
  rename(NormPopScore_ROH = NormPopScore.x, NormPopScore_IBD = NormPopScore.y) %>%
  mutate(Population = factor(Population, levels=orderPops$V1))

#Popularity data frame
PopularityDF = popmapMerge %>%
  select("breed", "clade") %>%
  group_by(breed) %>%
  sample_n(1) %>%
  mutate(OverallPopularityRank = AKC$popularity[match(breed, AKC$breed)]) %>%
  filter(!is.na(OverallPopularityRank)) %>%
  rename(Population=breed)

comboDF_noWolves = comboDF %>%  
  filter(!grepl("grayWolf",Clade)) %>% #remove labelled wolves
  na.omit() %>% #remove Stronen Wolves (no clade)
  mutate(OverallPopularityRank = PopularityDF$OverallPopularityRank[match(Population, PopularityDF$Population)])

#Add ROH and IBD Scores to the Disease Prevalance data
TWilesData = t(WilesData[,2:ncol(WilesData)]) %>%
  as.data.frame() #Transpose data
colnames(TWilesData) = WilesData[,1] #Set the column headings from the first column in the original table
DiseasePrev = TWilesData %>%
  rownames_to_column(var = "Population") %>%
  mutate(NormPopROHScore = comboDF_noWolves$NormPopScore_ROH[match(Population, comboDF_noWolves$Population)],
         NormPopIBDScore = comboDF_noWolves$NormPopScore_IBD[match(Population, comboDF_noWolves$Population)]) %>%
  na.omit()
rm(TWilesData)



x = DiseasePrev %>% 
  dplyr::select(-c("Population","NormPopROHScore", "NormPopIBDScore")) %>%  # exclude outcome, leave only predictors 
  map(~lm(DiseasePrev$NormPopROHScore ~ .x, data = DiseasePrev)) %>% 
  map(summary) %>% 
  map(c("coefficients")) %>% 
  map_dbl(8)  # 8th element is the p-value 

y = DiseasePrev %>% 
  dplyr::select(-c("NormPopROHScore", "NormPopIBDScore")) %>%  # exclude outcome, leave only predictors 
  map(~lm(DiseasePrev$NormPopROHScore ~ .x, data = DiseasePrev)) %>% 
  map(summary) %>% 
  map_dbl("r.squared") 


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
