#Load Libraries
library(tidyverse)
library(data.table)
library(ggplot2)
library(DescTools)

#Read Files in 
setwd("~/Documents/DogProject_Jaz/AKC")
AKC = read.table("~/Documents/DogProject_Jaz/AKC/AKC_1926thru2005_Ne.txt", check.names = F, header = T, fill = T)

#Make everything proportional
propAKC = apply(AKC[,-c(1)], 2, function(x) (as.numeric(x)/sum(x, na.rm = T))) %>% as.data.frame() #rank data as proportion of reg dogs

#Generate Rank based on popularity
#most popular with largest number and least popular with smallest number
rankedDF = apply(propAKC, 2, FUN=rank) %>% as.data.frame() #rank data

#Set rank to NA if AKC number is NA
AllColumns = colnames(rankedDF) #grab column names
for (i in AllColumns){
  newCol = ifelse(is.na(AKC[,i]), NA, rankedDF[,i])
  rankedDF[,i] = newCol
} 

#Add column with harmonic mean and add breed id
rankedDF$HarmonicMean = apply(rankedDF, 1, function(x)Hmean(x, na.rm = T))
FinalAKC = cbind.data.frame(AKC$Breed, rankedDF)
names(FinalAKC)[1] = "Breed"

#Find Overall Popularity and add it to Final DF
OverallPopularity = FinalAKC %>% select("Breed", "HarmonicMean") %>% filter(HarmonicMean != "NaN") %>% mutate(OverallPopularity = dense_rank(as.numeric(HarmonicMean)))
FinalAKC$OverallPopularity = OverallPopularity$OverallPopularity[match(FinalAKC$Breed, OverallPopularity$Breed)]

#write out
#write.table(FinalAKC, "AKC_breedPopularity_1926thru2005.txt", sep = "\t", row.names = F, col.names = T, quote = F)



