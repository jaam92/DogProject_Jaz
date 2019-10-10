#Plot manuscript figures
library(tidyverse)
library(ggpubr)

#Read Files in 
IBDScores = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/IBDSegs/IBDScoresPerPopulation.txt") %>%
  rename("Population" = Breed1)
ROHScores = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/ROH/ROHScoresPerPopulation.txt")
orderPops = read.table("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/OrderPops.txt")
#orderCluster = read.table("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/OrderCluster.txt")

#Rename column and merge
comboDF = merge(ROHScores, IBDScores, by ="Population", all = T)

#Set Populations and Clusters as factor so Wolves and dogs group together
IBDScores$Population = factor(IBDScores$Population, levels=orderPops$V1)
ROHScores$Population = factor(ROHScores$Population, levels=orderPops$V1)
comboDF$Population = factor(comboDF$Population, levels=orderPops$V1)
ROHScores$logRelToLab = log(ROHScores$RelToLab) #log ROH Scores

#Plot ROH Data
midROH = mean(ROHScores$logRelToLab)
PlotROHRel2Lab = ggplot(ROHScores, aes(y=logRelToLab, x=Population, fill=logRelToLab)) + 
  geom_bar(stat="identity") + 
  coord_flip() + theme_bw() + 
  labs(y = "Within Population ROH Score Relative to Labrador Retriever (Normalized by Sample Size)", x="Breed") + 
  geom_hline(yintercept = midROH, linetype="dashed", colour = "black")  + 
  geom_hline(yintercept = mid, linetype="dashed", colour = "black") + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 10), 
        plot.title=element_text(size=26, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_gradient2(name= "Log Fold-Enrichment",midpoint=midROH, low="blue", mid="white",high="red", space ="Lab") 

#Plot IBD Data
midIBD = mean(IBDScores$RelToLab)
PlotIBDRelToLab = ggplot(IBDScores, aes(y=RelToLab, x=Population, fill=RelToLab)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() + 
  labs(y = "Within Population IBD Score Relative to Labrador Retriever (Normalized by Sample Size)", x="Breed") + 
  geom_hline(yintercept = midIBD, linetype="dashed", colour = "black") + 
  geom_hline(yintercept = mid, linetype="dashed", colour = "black") + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 10), 
        plot.title=element_text(size=26, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_gradient2(name= "Fold-Enrichment", midpoint=midIBD, low="blue", mid="white",high="red", space ="Lab")

#Multiplot
ggarrange(PlotROHRel2Lab, 
         PlotIBDRelToLab,
         ncol = 2,
         labels=c("A","B"))