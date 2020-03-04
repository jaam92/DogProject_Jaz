#load libraries
library(tidyverse)
library(data.table)
library(ade4)
library(ggpubr)

#Load breed info, ROHs, and kinship matrix
popmapDryad = read.delim("~/DogProject_Jaz/LocalRscripts/BreedCladeInfo/breeds_dryad.txt")

rohs = read.delim(file = "~/DogProject_Jaz/LocalRscripts/ROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt") %>%
  mutate(breed = popmapDryad$breed[match(INDV, popmapDryad$dogID)],
         group = ifelse(is.na(breed), "remove", "dog")) 

pcRelateMat = readRDS("~/DogProject_Jaz/LocalRscripts/CaseControlROH/pcRelateMatrix_allIndivs.rds")

#Subset kinship matrix
#keep only those rows and columns that correspond to dogs of interest
kinshipMat = pcRelateMat %>%
  as.matrix()
kinshipMat = kinshipMat[, colnames(kinshipMat) %in% popmapDryad$dogID]
kinshipMat = kinshipMat[ rownames(kinshipMat) %in% popmapDryad$dogID,]
rm(pcRelateMat)
kinshipMatOrder = rownames(kinshipMat) 

#Find all possible pair of inidividuals
allIndivs = rohs %>%
  filter(group != "remove" & INDV%in%kinshipMatOrder) %>%
  select(INDV) %>%
  distinct() %>%
  mutate()

#make data frame with individuals of interest
ROHperIndiv = rohs %>%
  filter(INDV %in% allIndivs$INDV) %>%
  group_by(INDV) %>%
  summarise(ROHBurden = as.numeric(sum(AUTO_LEN))) %>% 
  ungroup() %>%
  mutate(ROHBurdenNorm = as.numeric((ROHBurden-min(ROHBurden))/(max(ROHBurden)-min(ROHBurden)))) %>%
  select(INDV, ROHBurdenNorm) 

#all pairwise combinations of individuals (order does not matter)
dat = t(combn(allIndivs$INDV, 2)) %>% 
  as.data.table() %>% 
  #distinct(V1, V2) %>%
  as.data.frame() %>%
  mutate(sim = 1-abs(ROHperIndiv$ROHBurdenNorm[match(V1, ROHperIndiv$INDV)] - ROHperIndiv$ROHBurdenNorm[match(V2, ROHperIndiv$INDV)]), #1 - absolute value of the difference in ROH burden between two individuals (to get similarity)
         breed1 = popmapDryad$breed[match(V1, popmapDryad$dogID)],
         breed2 = popmapDryad$breed[match(V2, popmapDryad$dogID)],
         indicator = ifelse(breed1 == breed2, "same", "different")) 

#Make a similarity matrix with ROH burden
#First, make a data frame with each individual to themselves 
self = allIndivs %>%
  mutate(V2 = INDV,
         sim = as.numeric(1)) %>%
  rename("V1" = "INDV")

#Make correlation matrix 
#Add self matrix to comparison matrix then reorder by order in kinship file
orderedData = dat %>%
  select(V1, V2, sim) %>%
  rbind.data.frame(self) %>%
  arrange(match(V1, kinshipMatOrder)) %>% 
  arrange(match(V2, kinshipMatOrder)) 

#Reshape the kinship matrix and normalize it by bounding between 0 and 1
longData = melt(kinshipMat) %>%
  mutate(V1 = as.character(Var1),
         V2 = as.character(Var2),
         kinshipNorm = as.numeric((value-min(value))/(max(value)-min(value)))) %>%
  select(V1, V2, kinshipNorm)

#plot the matrices
test = orderedData %>%
  inner_join(longData)

#ROH burden
rohHeatMap = ggplot(test, aes(x = V1, y = V2)) + 
  geom_raster(aes(fill = sim)) + 
  scale_fill_gradient(low="grey90", high="red", name = "Similarity") +
  labs(x="dogID1", y="dogID2", title = "ROH Burden") +
  theme_bw() + 
  theme(plot.title=element_text(size=18, face = "bold", hjust=0.5),
        axis.title=element_text(size=16),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

#Kinship
kinshipHeatMap = ggplot(test, aes(x = V1, y = V2)) + 
  geom_raster(aes(fill = kinshipNorm)) + 
  scale_fill_gradient(low="grey90", high="red", name = "Similarity") +
  labs(x="dogID1", y="dogID2", title = "Kinship") +
  theme_bw() + 
  theme(plot.title=element_text(size=18, face = "bold", hjust=0.5),
        axis.title=element_text(size=16),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
