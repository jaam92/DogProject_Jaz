#Plot manuscript figures
library(dplyr)
library(data.table)
library(cowplot)
library(ggplot2)

#Read Files in 
IBDScores = read.delim("~/Documents/DogProject_Jaz/IBDSegs/IBDScoresPerPopulation.txt")
ROHScores = read.delim("~/Documents/DogProject_Jaz/ROH/ROHScoresPerPopulation.txt")
orderPops = read.table("~/Documents/DogProject_Jaz/BreedCladeInfo/OrderPops.txt")
#orderCluster = read.table("~/Documents/DogProject_Jaz/BreedCladeInfo/OrderCluster.txt")

#Rename column and merge
names(IBDScores)[1] = "Population"
comboDF = merge(ROHScores, IBDScores, by ="Population", all = T)

#Set Populations and Clusters as factor so Wolves and dogs group together
IBDScores$Population = factor(IBDScores$Population, levels=orderPops$V1)
ROHScores$Population = factor(ROHScores$Population, levels=orderPops$V1)
comboDF$Population = factor(comboDF$Population, levels=orderPops$V1)
ROHScores$logRelToLab = log(ROHScores$RelToLab)

#Plot Data
midROH = mean(ROHScores$logRelToLab)
PlotROHRel2Lab = ggplot(ROHScores, aes(y=logRelToLab, x=Population, fill=logRelToLab)) + geom_bar(stat="identity") + coord_flip() + theme_bw() + labs(y = "Within Population ROH Score Relative to Labrador Retriever (Normalized by Sample Size)", x="Breed") + geom_hline(yintercept = midROH, linetype="dashed", colour = "black") + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 10), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=20)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=18)) + scale_y_continuous(expand = c(0,0)) + scale_fill_gradient2(name= "Log Fold-Enrichment",midpoint=midROH, low="blue", mid="white",high="red", space ="Lab") 


midIBD = mean(IBDScores$RelToLab)
PlotIBDRelToLab = ggplot(IBDScores, aes(y=RelToLab, x=Population, fill=RelToLab)) + geom_bar(stat="identity") + coord_flip() + theme_bw() + labs(y = "Within Population IBD Score Relative to Labrador Retriever (Normalized by Sample Size)", x="Breed") + geom_hline(yintercept = midIBD, linetype="dashed", colour = "black") + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 10), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=20)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=18)) + scale_y_continuous(expand = c(0,0)) + scale_fill_gradient2(name= "Fold-Enrichment",midpoint=midIBD, low="blue", mid="white",high="red", space ="Lab")

plot_grid(PlotROHRel2Lab, PlotIBDRelToLab)

#Plot REl to Wolf
PlotROHRel2EuroWolf = ggplot(ROHScores, aes(y=RelToEURO, x=Population)) + geom_bar(stat="identity") + coord_flip() + theme_bw() + labs(y = "Within Population ROH Score Relative to European Gray Wolf (Normalized by Sample Size)", x="Breed") + geom_hline(yintercept = 1, linetype="dashed", colour = "blue") + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 10), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=20)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=18)) + scale_y_continuous(expand = c(0,0))

PlotIBDRelToEuroWolf = ggplot(IBDScores, aes(y=RelToEURO, x=Population)) + geom_bar(stat="identity") + coord_flip() + theme_bw() + labs(y = "Within Population IBD Score Relative to European Gray Wolf (Normalized by Sample Size)", x="Breed") + geom_hline(yintercept = 1, linetype="dashed", colour = "blue") + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 10), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=20)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=18)) + scale_y_continuous(expand = c(0,0))