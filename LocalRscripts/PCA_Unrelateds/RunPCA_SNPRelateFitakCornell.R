####Based on the results from the PCA going to remove ID-14 wolf from all downstream analyses ####


#Download packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("SNPRelate")

#Load libraries
library(mgsub)
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
popmapDog$Type = ifelse(popmapDog$clade != "Village", "breed dog", "village dog")
popmapWolf = popmapMerge[grep("Wolf",popmapMerge$clade),]
popmapWolf$Type = ifelse(popmapWolf$clade !="grayWolf_Europe", "North American wolf", "European wolf")
popmapMaster = rbind.data.frame(popmapDog,popmapWolf)
popmapMaster$Type = ifelse(is.na(popmapMaster$Type), "breed dog", popmapMaster$Type) #fill all the N/As in for breed dogs
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

newOrder_Legend = c("breed dog", "village dog", "European wolf", "North American wolf")
df_PCA$population = factor(df_PCA$population, levels = newOrder_Legend)

#plot results with dogs and wolves
allSampsPC1vPC2 = ggplot(df_PCA, aes(y=EV2, x=EV1, colour=population)) +
  geom_point(size=2) +
  scale_colour_manual(values = c("breed dog"= "#56B4E9", "village dog"="black", "North American wolf" = "#E69F00", "European wolf" = "#009E73"), name="Population") + 
  labs(y=bquote('PC2' ~'('~.(pc[2])~'%'~')'), x=bquote('PC1'~'('~.(pc[1])~'%'~')')) +
  theme_bw() + 
  theme(axis.text.x = element_text(size  = 24), 
        axis.text.y = element_text(size  = 24), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=18),
        legend.text=element_text(size=18), 
        legend.position = "bottom")

allSampsPC2vPC3 = ggplot(df_PCA, aes(y=EV3, x=EV2, colour=population)) +
  geom_point(size=2)  +
  scale_colour_manual(values = c("breed dog"= "#56B4E9", "village dog"="black", "North American wolf" = "#E69F00", "European wolf" = "#009E73"), name="Population") + 
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

#Run just the dog data
gdsSampIDs = read.gdsn(index.gdsn(genofile, "sample.id"))
sampleid_dog = gdsSampIDs[!(grepl("grayWolf",gdsSampIDs))]
samp.id.dog = unlist(sampleid_dog)

#LD prune
snpset_dog = snpgdsLDpruning(genofile, sample.id = samp.id.dog, method="corr", slide.max.n=50, ld.threshold=0.3, maf = 0.05, autosome.only = F)
snpset.id.dog = unlist(snpset_dog)

#Run PCA
pca_dog = snpgdsPCA(genofile, snp.id = snpset.id.dog, sample.id = samp.id.dog, autosome.only = F)

#Look at percentage of variance explained by each PC
pc.percent_dog = pca_dog$varprop*100
pc_dog = head(round(pc.percent_dog,2))
pc_dog
pca_dog$sample.id = gsub('(.*)-\\w+', '\\1', pca_dog$sample.id)

#Make data frame with first two pcs
df_PCA_dog = data.frame(sample.id = pca_dog$sample.id, 
                        clade = factor(as.character(popmapMaster$clade))[match(pca_dog$sample.id, sample.id)],
                        breed = popmapMaster$breed[match(pca_dog$sample.id, sample.id)],
                        EV1 = pca_dog$eigenvect[,1], 
                        EV2 = pca_dog$eigenvect[,2], 
                        EV3 = pca_dog$eigenvect[,3], 
                        EV4 = pca_dog$eigenvect[,4], 
                        stringsAsFactors = FALSE) %>%
  mutate(clade = gsub("_"," ", clade),
         clade = replace_na(clade, "Unassigned"))

#Evaluate data frame with breed and clade information
#expand color palette of choice to have number of colors equal to number of clades
palette = c("#6DDFD1","#74C5E0","#DDE651","#9BA577","#AB779F","#6EE29A","#78E053","#ABE3BD","#5763E7","#F12364","#E8DEC1","#942BF7","#E53FDE","#70A1E5","#D1967D","#D5A6E6","#B75ACA","#F087C5","#EA4AB5","#D34785","#634C96","#C7E4E3","#CDC6E4","#9F8AE6","#E9B0C1","#E887EC","#9A4EE0","#EA6B37","#DDAE56","#CDE18D","#D4615F","#83909A")

#plot the nice version
allSampsPC1vPC2byClade = ggplot(df_PCA_dog, aes(y=EV2, x=EV1, colour=clade)) +
  geom_point(size=2) + 
  scale_colour_manual(values = palette, na.value = "black", name = "Clade") +
  labs(y=bquote('PC2' ~'('~.(pc_dog[2])~'%'~')'), x=bquote('PC1'~'('~.(pc_dog[1])~'%'~')')) +
  theme_bw() + 
  theme(axis.text.x = element_text(size  = 24), 
        axis.text.y = element_text(size  = 24), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=18),
        legend.text=element_text(size=18), 
        legend.position = "bottom")

print(allSampsPC1vPC2byClade)


#color by roh burden
rohs = read.delim(file = "~/Documents/DogProject_Jaz/LocalRscripts/ROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt") %>%
  group_by(INDV) %>%
  summarise(ROHBurden = as.numeric(sum(AUTO_LEN))) %>% 
  ungroup() %>%
  mutate(ROHBurdenNorm = as.numeric((ROHBurden-min(ROHBurden))/(max(ROHBurden)-min(ROHBurden)))) %>%
  select(INDV, ROHBurdenNorm) %>%
  mutate(breed = popmapMaster$breed[match(INDV, popmapMaster$dogID)])

df_PCA_dog$ROHBurden = rohs$ROHBurdenNorm[match(df_PCA_dog$sample.id, rohs$INDV)]

rohBurden = ggplot(df_PCA_dog, aes(y=EV2, x=EV1, colour=ROHBurden)) +
  geom_point(size=2) + 
  scale_colour_gradientn(colours = terrain.colors(10), name = "ROH Burden") +
  labs(y=bquote('PC2' ~'('~.(pc_dog[2])~'%'~')'), x=bquote('PC1'~'('~.(pc_dog[1])~'%'~')')) +
  theme_bw() + 
  theme(axis.text.x = element_text(size  = 24), 
        axis.text.y = element_text(size  = 24), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=18),
        legend.text=element_text(size=18), 
        legend.position = "right")

print(rohBurden)

###PCA on wolf populations###
sampleid_wolf = gdsSampIDs[grep("grayWolf", gdsSampIDs)]
samp.id.wolf = unlist(sampleid_wolf)

#LD prune
snpset_wolf = snpgdsLDpruning(genofile, sample.id = samp.id.wolf, method="corr", slide.max.n=50, ld.threshold=0.3, maf = 0.05, autosome.only = F)
snpset.id.wolf = unlist(snpset_wolf)

#Run PCA
pca_wolf = snpgdsPCA(genofile, snp.id = snpset.id.wolf, sample.id = samp.id.wolf, autosome.only = F)

#Look at percentage of variance explained by each PC
pc.percent_wolf = pca_wolf$varprop*100
pc_wolf = head(round(pc.percent_wolf,2))
pc_wolf
pca_wolf$sample.id = gsub('(.*)-\\w+', '\\1', pca_wolf$sample.id)

#Make data frame with first two pcs
df_PCA_Wolf = data.frame(sample.id = pca_wolf$sample.id, 
                            population = factor(as.character(popmapMaster$Type))[match(pca_wolf$sample.id, sample.id)], 
                            clade = factor(as.character(popmapMaster$clade))[match(pca_wolf$sample.id, sample.id)], 
                            location = popmapMaster$breed[match(pca_wolf$sample.id, sample.id)],
                            EV1 = pca_wolf$eigenvect[,1], 
                            EV2 = pca_wolf$eigenvect[,2], 
                            EV3 = pca_wolf$eigenvect[,3], 
                            EV4 = pca_wolf$eigenvect[,4], 
                            stringsAsFactors = FALSE) 

#replace location labels
df_PCA_Wolf$location = mgsub(df_PCA_Wolf$location, pattern=c("EURO", "WO_BC", "WO_ID", "WO_INTAK", "WO_MN", "WO_MAT", "X", "MB", "WO_SEAK", "WO_WO", "WO_LUPA", "GR", "AR"), replacement=c("Europe", "British Columbia", "Idaho", "Interior Alaska", "Minnesota", "Montana", "Mexico","McBride", "Southeast Alaska", "Wyoming", "Europe", "Ghost Ranch", "Aragon"))
"#D3DDD5"

#plot the nice version
wolfPalette =c("#7BA6D6","#D6A9D0","firebrick2","#78D7D6","purple","#896ED7","dodgerblue2","royalblue4","forestgreen","#DD60B2","#C9A981","darkmagenta")

PC1vPC2_wolf = ggplot(df_PCA_Wolf, aes(y=EV2, x=EV1, colour=population)) +
  geom_point(size=2) +
  labs(y=bquote('PC2' ~'('~.(pc_wolf[2])~'%'~')'), x=bquote('PC1'~'('~.(pc_wolf[1])~'%'~')')) +
  scale_colour_manual(values = c("North American wolf" = "gray80", "European wolf" = "blue"), name="Population") +
  theme_bw() +
  theme(axis.text.x = element_text(size  = 24), 
        axis.text.y = element_text(size  = 24), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=24), 
        legend.text=element_text(size=24))

PC1vPC2_wolf = ggplot(df_PCA_Wolf, aes(y=EV3, x=EV2, colour=location))+
  geom_point(size=2) +
  labs(y=bquote('PC3' ~'('~.(pc_wolf[3])~'%'~')'), x=bquote('PC2'~'('~.(pc_wolf[2])~'%'~')')) +
  scale_colour_manual(values = wolfPalette, na.value = "black", name = "Location") +
  theme_bw() +
  theme(axis.text.x = element_text(size  = 24), 
        axis.text.y = element_text(size  = 24), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=24), 
        legend.text=element_text(size=24),
        legend.position = "bottom")

print(PC1vPC2_wolf)
print(PC2vPC3_wolf)
