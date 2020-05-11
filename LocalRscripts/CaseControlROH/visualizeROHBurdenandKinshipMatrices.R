#load libraries
library(tidyverse)
library(data.table)
library(ade4)
library(ggpubr)

#Functions
#fxn for plotting
ggplotRegression = function (fit) {
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(size = 2) + 
    stat_smooth( method = 'lm', col = "blue") +  
    theme_bw() + 
    labs(title = bquote(R^2== ~.(signif(summary(fit)$adj.r.squared, 5))~"&"~"p"==~.(signif(summary(fit)$coef[2,4], 5))))
  
}

#Load breed info, ROHs, and kinship matrix
popmapDryad = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/breeds_dryad.txt")

rohs = read.delim(file = "~/Documents/DogProject_Jaz/LocalRscripts/ROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFile_Cornell_allChroms_vcfToolsROH_rmROHlessThan50snps_HaywardDataOnly.txt") %>%
  mutate(breed = popmapDryad$breed[match(INDV, popmapDryad$dogID)],
         group = ifelse(is.na(breed), "remove", "dog")) 

pcRelateMat = readRDS("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/pcRelateMatrix_allIndivs.rds")

grmROHs = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/ROHGRM.ped", stringsAsFactors = F) %>% 
  mutate(rohGRMNorm = (Overlaps-min(Overlaps))/(max(Overlaps)-min(Overlaps))) %>%
  #select(-c(Overlaps)) %>%
  rename(V1 = ID1, V2 = ID2)

#Compute bounded ROH burden for each individual
ROHperIndiv = rohs %>%
  group_by(INDV) %>%
  summarise(ROHBurden = as.numeric(sum(AUTO_LEN))) %>% 
  ungroup() %>%
  mutate(ROHBurdenNorm = as.numeric((ROHBurden-min(ROHBurden))/(max(ROHBurden)-min(ROHBurden)))) %>%
  select(INDV, ROHBurdenNorm) 

#Subset kinship matrix
#keep only those rows and columns that correspond to dogs of interest
kinshipMat = pcRelateMat %>%
  as.matrix()
kinshipMat = kinshipMat[, colnames(kinshipMat) %in% popmapDryad$dogID]
kinshipMat = kinshipMat[ rownames(kinshipMat) %in% popmapDryad$dogID,]
rm(pcRelateMat)
kinshipMatOrder = rownames(kinshipMat) 

#Reshape the kinship matrix and normalize it by bounding between 0 and 1
longData = melt(kinshipMat) %>%
  mutate(V1 = as.character(Var1),
         V2 = as.character(Var2),
         kinshipNorm = as.numeric((value-min(value))/(max(value)-min(value))), 
         rohNorm = 1-abs(ROHperIndiv$ROHBurdenNorm[match(V1, ROHperIndiv$INDV)] - ROHperIndiv$ROHBurdenNorm[match(V2, ROHperIndiv$INDV)])) %>%
  select(V1, V2, kinshipNorm, rohNorm)

#correlate them ROH and kinship
#linear.model = lm(rohNorm ~ kinshipNorm, data = longData)
#summary(linear.model)

#ggplotRegression(linear.model)

#Add kinship based on shared ROH
mergedDF = longData %>%
  select(-c(kinshipNorm)) %>%
  left_join(grmROHs) %>%
  mutate(breed1 = popmapDryad$breed[match(V1,popmapDryad$dogID)],
         breed2 = popmapDryad$breed[match(V2,popmapDryad$dogID)]) %>%
  arrange(breed1)

mergedDF[is.na(mergedDF)] <- as.numeric(1)

#Kinship
kinshipHeatMap = ggplot(mergedDF, aes(x = V1, y = V2)) + 
  geom_raster(aes(fill = grmROHs)) + 
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

#ROH burden
rohHeatMap = ggplot(mergedDF, aes(x = V1, y = V2)) + 
  geom_raster(aes(fill = rohNorm)) + 
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

#plot them together
ggarrange(rohHeatMap, kinshipHeatMap, 
          nrow = 1, 
          common.legend = T, 
          legend = "bottom")


#These lines of code will get you all the pairwise compairsons in the hayward data and adds the file extension to the end so that the file names can be fed into bedtools
#x = longData %>% 
#  select(V1,V2) %>% 
#  filter(V1 != V2) %>% 
#  mutate(V1 = paste0(V1, ".bed"), 
#         V2 = paste0(V2, ".bed"), 
#         V1 = gsub("-N/A", "", V1), 
#         V2 = gsub("-N/A", "", V2)) 
#write.table(x, file = "~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/fnames_haywardComps.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#convert to distance matrices and mantel test
#finalROHMatrix = longData %>%
#  select(V1, V2, rohNorm) %>%
#  acast(V1~V2, value.var="rohNorm")

#roh.dist = dist(finalROHMatrix)
#kinship.dist = dist(kinshipMat)
#mantel.rtest(kinship.dist, roh.dist, nrepet = 1000)
