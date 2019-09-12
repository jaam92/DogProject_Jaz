#Set working directory and load libraries
setwd("~/Documents/DogProject_Jaz/LocalRscripts/ROH")
library(ggplot2)
library(reshape2)
library(dplyr)
library(randomcoloR)
library(data.table)
library(ggpubr)
library(cowplot)

#Read Files in
dfMerge = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/ROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt")
dfStronen = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/ROH/Stronen2013_WolvesTrueROH_propCoveredwithin1SDMean_allChroms.txt")
popmapMerge = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/BreedAndCladeInfo_mergedFitakCornell.txt")
popmapStronen = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/Stronen2013_SamplesUsed.txt")
orderPops = read.table("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/OrderPops.txt")
orderCluster = read.table("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/OrderCluster.txt")

#Concatenate dataframes
WolfDog = rbind.data.frame(dfMerge,dfStronen)

#Split into length intervals, find length for each indiv, add column in Gb, add range column
aggIntROH = WolfDog %>% select(AUTO_LEN, INDV) %>% filter(AUTO_LEN < 1000000) %>% group_by(INDV) %>% summarise(totalLen = sum(AUTO_LEN)) %>% mutate(totalLenGb = totalLen/10^9, Range = "[0.1-1)Mb") %>% as.data.frame()

aggLongROH = WolfDog %>% select(AUTO_LEN, INDV) %>% filter(AUTO_LEN >= 1000000 & AUTO_LEN < 10000000) %>% group_by(INDV) %>% summarise(totalLen = sum(AUTO_LEN)) %>% mutate(totalLenGb = totalLen/10^9, Range = "[1-10)Mb")  %>% as.data.frame()

aggVLongROH = WolfDog %>% select(AUTO_LEN, INDV) %>% filter(AUTO_LEN >= 10000000) %>% group_by(INDV) %>% summarise(totalLen = sum(AUTO_LEN)) %>% mutate(totalLenGb = totalLen/10^9, Range = "[10-63)Mb") %>% as.data.frame()

#Merged interval data frames
FinalDF = rbind(aggIntROH, aggLongROH, aggVLongROH)
FinalDF$Range = factor(FinalDF$Range, levels=c("[10-63)Mb","[1-10)Mb","[0.1-1)Mb"))

#Assign breed,clade,cluster etc. and bind dog and wolf together
shortPopmapWolves = popmapStronen %>% select(IID_Col2,Population,Cluster)
#shortPopmapWolves$Population = shortPopmapWolves$Cluster #make the cluster the population for wolves
names(shortPopmapWolves)[1] = "dogID"
names(shortPopmapWolves)[2] = "breed"
names(shortPopmapWolves)[3] = "clade"
mergedPopmap = rbind.data.frame(popmapMerge, shortPopmapWolves)
CountsPerBreed = popmapMerge %>% group_by(breed) %>% tally()
CountsPerCluster = popmapMerge %>% group_by(clade) %>% tally()

##make sure names are unique (they are)
#length(unique(mergedPopmap$dogID))
#length(mergedPopmap$dogID)

#Add Population and Cluster information
FinalDF$Population = mergedPopmap$breed[match(FinalDF$INDV,mergedPopmap$dogID)]
FinalDF$Cluster = mergedPopmap$clade[match(FinalDF$INDV,mergedPopmap$dogID)]

#Remove inidivudals without Cluster and outliers
FinalDF_rmNA = FinalDF %>% filter(!is.na(Cluster) & Cluster!= "Outlier")

#Set Populations and Clusters as factor so Wolves and dogs group together
FinalDF_rmNA$Population = factor(FinalDF_rmNA$Population, levels=orderPops$V1)
FinalDF_rmNA$Cluster = factor(FinalDF_rmNA$Cluster, levels=orderCluster$V1)

#Find mean per cluster and population
MeanPerPopulationDF = FinalDF_rmNA %>% select(Population,Range,totalLen) %>% group_by(Population,Range) %>% summarise_all(mean) %>%as.data.frame()
MeanPerClusterDF = FinalDF_rmNA %>% select(Cluster,Range,totalLen) %>% group_by(Cluster,Range) %>% summarise_all(mean) %>%as.data.frame()

#Plot by Cluster
MeanROHperPopulation = ggplot(MeanPerPopulationDF, aes(x=Population, y=totalLen/10^6, fill=Range)) + geom_bar(stat="identity") + theme_bw() + coord_flip() + scale_fill_manual(values = c("[0.1-1)Mb"= "yellow", "[10-63)Mb" = "red", "[1-10)Mb" = "orange"), breaks = c("[0.1-1)Mb","[1-10)Mb", "[10-63)Mb"), name = "Range") + ylab("Mean ROH Length per Bin (Mb)") + xlab("Breed") + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 10), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=20)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=18))

MeanROHperCluster = ggplot(MeanPerClusterDF, aes(x=Cluster, y=totalLen/10^6, fill=Range)) + geom_bar(stat="identity") + theme_bw() + coord_flip() + scale_fill_manual(values = c("[0.1-1)Mb"= "yellow", "[10-63)Mb" = "red", "[1-10)Mb" = "orange"), breaks = c("[0.1-1)Mb","[1-10)Mb", "[10-63)Mb"), name = "Range") + ylab("Mean ROH Length per Bin (Mb)") + xlab("Clade") + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 14), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=20))

#####Plot FROH
FROHdf = WolfDog %>% select(AUTO_LEN, INDV) %>% group_by(INDV) %>% summarise(totalLen = sum(AUTO_LEN)) %>% mutate(FROH = totalLen/2500000000) %>% as.data.frame() #Dog genome is 2.5 gigabases
FROHdf$Population = mergedPopmap$breed[match(FROHdf$INDV,mergedPopmap$dogID)]
FROHdf$Cluster = mergedPopmap$clade[match(FROHdf$INDV,mergedPopmap$dogID)]
FROHdf_rmNA = FROHdf %>% filter(!is.na(Cluster) & Cluster!= "Outlier")

#Set Populations and Clusters as factor so Wolves and dogs group together
FROHdf_rmNA$Population = factor(FROHdf_rmNA$Population, levels=orderPops$V1)
FROHdf_rmNA$Cluster = factor(FROHdf_rmNA$Cluster, levels=orderCluster$V1)

#expand color palette of choice to hve number of colors equal to number of clades
colourCount_pop = length(unique(FROHdf_rmNA$Population)) 
palette = distinctColorPalette(colourCount_pop)

#Plot
FROH_perCluster = ggplot(FROHdf_rmNA, aes(x=Cluster, y=FROH, colour=Cluster)) + geom_boxplot(size=1) + geom_point(size=0.5) + scale_colour_manual(values = palette, na.value="grey") + theme_bw() + labs(x = "Cluster", y=expression(F[ROH])) + theme(axis.text.x = element_text(size  = 24,angle=40, vjust=1, hjust=1), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=24), legend.text=element_text(size=18), legend.position = "none") + guides(fill = guide_legend(nrow = 4))

FROH_perPopulation = ggplot(FROHdf_rmNA, aes(x=Population, y=FROH, colour=Population)) + geom_boxplot(size=1) + geom_point(size=0.5) + scale_colour_manual(values = palette, na.value="grey") + theme_bw() + labs(x = "Population", y=expression(F[ROH])) + theme(axis.text.x = element_text(size  = 10, angle = 40, vjust=1, hjust=1), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=24), legend.text=element_text(size=18), legend.position = "none") + guides(fill = guide_legend(nrow = 4))


#Multiplot the ROH and FROH with Clade
#ggarrange(FROH_perPopulation + xlab(NULL), ggarrange(MeanROHperPopulation, MeanROHperCluster + xlab(NULL), ncol = 2, labels = c("B", "C"),common.legend = TRUE, legend = "right"), nrow = 2,labels = "A") 
plot_grid(FROH_perCluster + theme(axis.text.x = element_text(size=18)), MeanROHperCluster)


