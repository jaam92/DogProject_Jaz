####Based on the results from the PCA going to remove ID-14 wolf from all downstream analyses ####


#Download packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("gdsfmt")
#biocLite("SNPRelate")

#Load libraries
#library(GWASTools)
library(gdsfmt)
library(SNPRelate)
library(ggpubr)
library(randomcoloR)
library(tidyverse)

#Load Sample information
setwd("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH")
popmapMerge = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/BreedAndCladeInfo_mergedFitakCornell.txt")
orderPops = read.table("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/OrderPops.txt")

#PLINK BED files for unpruned data aka random sample of 100000 snps
#bed.fn = ("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/MergedFile_CornellCanineFitak_allIndivs.bed")
#bim.fn = ("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/MergedFile_CornellCanineFitak_allIndivs.bim")
#fam.fn = ("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/MergedFile_CornellCanineFitak_allIndivs.fam")

#convert
#snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/MergedFile_CornellCanineFitak_allIndivs.gds")
#snpgdsSummary("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/MergedFile_CornellCanineFitak_allIndivs.gds")

#Open file
genofile = snpgdsOpen("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/MergedFile_CornellCanineFitak_allIndivs.gds")
sampIds = read.gdsn(index.gdsn(genofile, "sample.id")) #grab sample ids 
famIds = gsub(".*-","",sampIds) #make family ids

#LD prune
snpset = snpgdsLDpruning(genofile, sample.id=sampIds, method="corr", slide.max.n=50, ld.threshold=0.3, maf = 0.05, autosome.only = F)
snpset.id = unlist(snpset)

#Run PCA with SNPs at frequecy higher than 5% and missingness less than 10%
pca = snpgdsPCA(genofile, snp.id=snpset.id, autosome.only = F, num.thread = 2)

#Look at percentage of variance explained by each PC
pc.percent = pca$varprop*100
pc = head(round(pc.percent,2))
pc

#Refomat 
pca$sample.id = gsub('(.*)-\\w+', '\\1', pca$sample.id)
pca$sample.id = gsub('-chinese_shar', '', pca$sample.id)
pca$sample.id = gsub('-curly', '', pca$sample.id)
pca$sample.id = gsub('-flat', '', pca$sample.id)

#Reformat the population map files
popmapMerge$breed = gsub("large_munsterlander","munsterlander_large", popmapMerge$breed)
popmapDog = popmapMerge[!(grepl("Wolf",popmapMerge$clade)),]
popmapDog$Type = ifelse(popmapDog$clade != "Village", "Breed dog", "Village dog")
popmapWolf = popmapMerge[grep("Wolf",popmapMerge$clade),]
popmapWolf$Type = ifelse(popmapWolf$clade !="grayWolf_Europe", "North American wolf", "European wolf")
popmapMaster = rbind.data.frame(popmapDog,popmapWolf)
sample.id = as.character(popmapMerge$dogID)
population = as.character(popmapMerge$breed)

#Make data frame with first four pcs
df_PCA = data.frame(sample.id = pca$sample.id, 
                    population = factor(as.character(popmapMaster$Type))[match(pca$sample.id, sample.id)], 
                    clade = factor(as.character(popmapMaster$clade))[match(pca$sample.id, sample.id)],
                    EV1 = pca$eigenvect[,1], 
                    EV2 = pca$eigenvect[,2], 
                    EV3 = pca$eigenvect[,3], 
                    EV4 = pca$eigenvect[,4],stringsAsFactors = FALSE)

#plot the nice version
allSampsPC1vPC2 = ggplot(df_PCA, aes(y=EV2, x=EV1, colour=population)) +
  geom_point(size=2) +
  scale_colour_manual(values = c("Breed dog"= "red", "Village dog"="black", "North American wolf" = "gray80", "European wolf" = "blue"), name="Population") + 
  labs(y=bquote('PC2' ~'('~.(pc[2])~'%'~')'), x=bquote('PC1'~'('~.(pc[1])~'%'~')')) +
  theme_bw() + 
  theme(axis.text.x = element_text(size  = 24), 
        axis.text.y = element_text(size  = 24), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=18),
        legend.text=element_text(size=18), 
        legend.position = "bottom")

allSampsPC2vPC3 = ggplot(df_PCA, aes(y=EV3, x=EV2, colour=population)) +
  geom_point(size=2) +
  scale_colour_manual(values = c("Breed dog"= "red", "Village dog"="black", "North American wolf" = "gray80", "European wolf" = "blue"), name="Population") + 
  labs(y=bquote('PC3' ~'('~.(pc[3])~'%'~')'), x=bquote('PC2'~'('~.(pc[2])~'%'~')')) +
  theme_bw() + 
  theme(axis.text.x = element_text(size  = 24), 
        axis.text.y = element_text(size  = 24), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=18),
        legend.text=element_text(size=18), 
        legend.position = "bottom")

print(allSampsPC1vPC2)
print(allSampsPC2vPC3)

#Evaluate data frame with breed and clade information
#expand color palette of choice to have number of colors equal to number of clades
colourCount = length(unique(df_PCA$clade))
palette = distinctColorPalette(colourCount)

#plot the nice version
allSampsPC1vPC2byClade = ggplot(df_PCA, aes(y=EV2, x=EV1, colour=clade)) +
  geom_point(size=2) + 
  scale_colour_manual(values = palette) +
  labs(y=bquote('PC2' ~'('~.(pc[2])~'%'~')'), x=bquote('PC1'~'('~.(pc[1])~'%'~')')) +
  theme_bw() + 
  theme(axis.text.x = element_text(size  = 24), 
        axis.text.y = element_text(size  = 24), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=18),
        legend.text=element_text(size=18), 
        legend.position = "bottom")

print(allSampsPC1vPC2byClade)

ggarrange(allSampsPC1vPC2 + theme(legend.text=element_text(size=12)), 
          allSampsPC1vPC2byClade + theme(legend.text=element_text(size=12)))

###PCA on North American populations###
gdsSampIDs = read.gdsn(index.gdsn(genofile, "sample.id"))
sampleid_American = gdsSampIDs[grep("_NorthAmerica", gdsSampIDs)]
samp.id.America = unlist(sampleid_American)

#LD prune
snpset_America = snpgdsLDpruning(genofile, sample.id = samp.id.America, method="corr", slide.max.n=50, ld.threshold=0.3, maf = 0.05, autosome.only = F)
snpset.id.America = unlist(snpset_America)

#Run PCA
pca_America = snpgdsPCA(genofile, snp.id = snpset.id.America, sample.id = samp.id.America, autosome.only = F)

#Look at percentage of variance explained by each PC
pc.percent_America = pca_America$varprop*100
pc_America = head(round(pc.percent_America,2))
pc_America
pca_America$sample.id = gsub('(.*)-\\w+', '\\1', pca_America$sample.id)

#Make data frame with first two pcs
df_PCA_America = data.frame(sample.id = pca_America$sample.id, 
                            population = factor(as.character(popmapMaster$Type))[match(pca_America$sample.id, sample.id)], 
                            clade = factor(as.character(popmapMaster$clade))[match(pca_America$sample.id, sample.id)], 
                            EV1 = pca_America$eigenvect[,1], 
                            EV2 = pca_America$eigenvect[,2], 
                            EV3 = pca_America$eigenvect[,3], 
                            EV4 = pca_America$eigenvect[,4], 
                            stringsAsFactors = FALSE)

#plot the nice version
PC1vPC2_America = ggplot(df_PCA_America, aes(y=EV2, x=EV1, colour=population)) +
  geom_point(size=2) +
  labs(y=bquote('PC2' ~'('~.(pc_America[2])~'%'~')'), x=bquote('PC1'~'('~.(pc_America[1])~'%'~')')) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(size  = 24), 
        axis.text.y = element_text(size  = 24), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=24), 
        legend.text=element_text(size=24))

PC2vPC3_America = ggplot(df_PCA_America, aes(y=EV3, x=EV2, colour=population))+
  geom_point(size=2) +
  labs(y=bquote('PC3' ~'('~.(pc_America[3])~'%'~')'), x=bquote('PC2'~'('~.(pc_America[2])~'%'~')')) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(size  = 24), 
        axis.text.y = element_text(size  = 24), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=24), 
        legend.text=element_text(size=24))

print(PC1vPC2_America)
print(PC2vPC3_America)

