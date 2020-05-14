#load libraries
library(tidyverse)
library(reshape2)
library(data.table)
library(ade4)
library(ggpubr)
library(mgsub)

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

#add N/A back to the IDs of Boxers
grmROHs = read.table("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/ROHGRM.ped", stringsAsFactors = F) %>%
  mutate(V1 = str_replace_all(grmROHs$V1, c("A4-B2-FEMALE" = "A4-B2-FEMALE-N/A","B1-B13-FEMALE" = "B1-B13-FEMALE-N/A", "C3-B4-FEMALE" = "C3-B4-FEMALE-N/A", "C5-HUC14-MALE" = "C5-HUC14-MALE-N/A", "D4-HUC3-FEMALE" = "D4-HUC3-FEMALE-N/A", "E5-HUC6-FEMALE" = "E5-HUC6-FEMALE-N/A", "F1-HUC13-MALE" = "F1-HUC13-MALE-N/A","F3-HUC16-MALE"= "F3-HUC16-MALE-N/A", "G1-HUC7-FEMALE" = "G1-HUC7-FEMALE-N/A", "H5-HUC12-MALE" = "H5-HUC12-MALE-N/A")),
         V2 = str_replace_all(grmROHs$V2, c("A4-B2-FEMALE" = "A4-B2-FEMALE-N/A","B1-B13-FEMALE" = "B1-B13-FEMALE-N/A", "C3-B4-FEMALE" = "C3-B4-FEMALE-N/A", "C5-HUC14-MALE" = "C5-HUC14-MALE-N/A", "D4-HUC3-FEMALE" = "D4-HUC3-FEMALE-N/A", "E5-HUC6-FEMALE" = "E5-HUC6-FEMALE-N/A", "F1-HUC13-MALE" = "F1-HUC13-MALE-N/A","F3-HUC16-MALE"= "F3-HUC16-MALE-N/A", "G1-HUC7-FEMALE" = "G1-HUC7-FEMALE-N/A", "H5-HUC12-MALE" = "H5-HUC12-MALE-N/A")))

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
  #left join the data with count of total bp of shared ROH
  #replace NA with approximate length of genome so comparing self to self will be 1
  #bound sharing of genome between 0 and 1
  #last add breed info
mergedDF = longData %>%
  select(-c(kinshipNorm)) %>%
  left_join(grmROHs) %>%
  mutate(V3 = ifelse(is.na(V3), as.numeric(3e9), V3), 
         rohGRMNorm = (V3-min(V3))/(max(V3)-min(V3)),
         breed1 = popmapDryad$breed[match(V1,popmapDryad$dogID)],
         breed2 = popmapDryad$breed[match(V2,popmapDryad$dogID)]) %>%
  select(-c(V3)) %>%
  arrange(breed1)

#Kinship
kinshipHeatMap = ggplot(mergedDF, aes(x = V1, y = V2)) + 
  geom_raster(aes(fill = rohGRMNorm)) + 
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
finalgrmROHMatrix = mergedDF %>%
  select(V1, V2, rohGRMNorm) %>%
  acast(V1~V2, value.var="rohGRMNorm")

#saveRDS(finalgrmROHMatrix, "~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/rohGRM_allIndivs.rds")

finalROHBurdenMatrix = mergedDF %>%
  select(V1, V2, rohNorm) %>%
  acast(V1~V2, value.var="rohNorm")

roh.dist = dist(finalROHBurdenMatrix)
kinship.dist = dist(finalgrmROHMatrix)
mantel.rtest(kinship.dist, roh.dist, nrepet = 1000)
