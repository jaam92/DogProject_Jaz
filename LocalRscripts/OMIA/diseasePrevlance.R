#####Load Libraries
library(tidyverse)

######Read Files in
WilesData = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/OMIA/AKC_DiseasePrev_PointEstWiles2017.txt", check.names = F, stringsAsFactors = F)
IBDScores = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/IBDSegs/IBDScoresPerPopulation.txt", stringsAsFactors = F)
ROHScores = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/ROH/ROHScoresPerPopulation.txt", stringsAsFactors = F)
AKC = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/AKC/AKC_breedPopularity_1926thru2005.txt", check.names = F, stringsAsFactors = F)
popmapMerge = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/BreedAndCladeInfo_mergedFitakCornell.txt")
orderPops = read.table("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/OrderPops.txt")


######Function to run linear regressions the summaryStat should be ROH or IBD depending on predictor you want to use #### 
LinearRegResults = function(summaryStat){
  if (summaryStat == "ROH") {
    Rsquared = DiseasePrev %>% 
      select(-c("Population","NormPopROHScore", "NormPopIBDScore")) %>%  # exclude predictors and population, leave only outcomes
      map(~lm(.x ~ DiseasePrev$NormPopROHScore, data = DiseasePrev)) %>% 
      map(summary) %>% 
      map_dbl("adj.r.squared") %>%
      as.data.frame() %>%
      rownames_to_column("trait") %>%
      rename("AdjRsquared" = ".") 
  }else{
    Rsquared = DiseasePrev %>% 
      select(-c("Population","NormPopROHScore", "NormPopIBDScore")) %>% 
      map(~lm(.x ~ DiseasePrev$NormPopIBDScore, data = DiseasePrev)) %>% 
      map(summary) %>% 
      map_dbl("adj.r.squared") %>%
      as.data.frame() %>%
      rownames_to_column("trait") %>%
      rename("AdjRsquared" = ".")
  } #decide on ROH or IBD
  
  #Merge the adjusted R-squared with p-value
  LinearReg = DiseasePrev %>% 
    select(-c("Population","NormPopROHScore", "NormPopIBDScore")) %>%  
    map(~lm(.x ~ DiseasePrev$NormPopROHScore, data = DiseasePrev)) %>% 
    map(summary) %>% 
    map(c("coefficients")) %>% 
    map_dbl(8) %>% #the 8th value is the p-value in the coefficients
    as.data.frame() %>%
    rownames_to_column("trait") %>%
    rename("pvalue" = ".") %>%
    left_join(Rsquared)
  
  return(LinearReg)
}
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

LinearRegROH = LinearRegResults("ROH")
LinearRegBoth = LinearRegResults("IBD") %>%
  mutate(Predictor = "IBD") %>%
  rbind.data.frame(LinearRegROH)

#Plot Linear Regression Rsquared and pvalue
ggplot(LinearRegBoth %>% filter(pvalue < 0.2), aes(x=reorder(trait, -AdjRsquared), y=AdjRsquared, colour=Predictor, size=pvalue<0.05)) +
  geom_point() +
  theme_bw() + 
  scale_colour_manual( values = c("ROH"="#E69F00", "IBD"="#56B4E9")) + 
  labs(x= "Trait", y = bquote('Adjusted'~R^2)) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 0.5, vjust=0.7),
        axis.text.y = element_text(size = 16), 
        plot.title=element_text(size=20, face = "bold", hjust=0.5), 
        axis.title=element_text(size=18),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18))

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
