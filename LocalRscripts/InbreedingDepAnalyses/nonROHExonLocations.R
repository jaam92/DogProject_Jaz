#load libraries 
library("tidyverse")

#set working directory and load files
setwd("~/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/")

ensemblGenes = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/ForAbi_EnsemblGenes_CanFam3.1_SingleTranscript.bed", col.names = c("chrom", "exonStart", "exonStop","GeneName"), stringsAsFactors = F)

#find the space in between exons within a gene and length of exons
ExonSpace = ensemblGenes %>%
  group_by(GeneName) %>%
  mutate(spaceBtwn = exonStart - lag(exonStop, default = first(exonStart)),
         exonLength = exonStop - exonStart) %>%
  as.data.frame()

TSSLength = ExonSpace %>% 
  group_by(chrom, GeneName) %>% 
  summarise(TSSLength = sum(exonLength)) 

#Figure out which exons do not have ROH
regionNonOverlaps_V = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/VettingResults/ExonRegion_NonOverlapsROH_cornellData.bed",  col.names = c("chrom", "exonStart", "exonStop","GeneName"), stringsAsFactors = F) %>%
        mutate(nonROH = "1") %>%
  right_join(., ensemblGenes) %>%
  mutate(data = "VCFTools",
         nonROH = ifelse(is.na(nonROH), as.numeric(0), nonROH)) #replace NA with 0 

regionNonOverlaps = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/VettingResults/ExonRegion_NonOverlapsROH_cornellData_plink.bed",  col.names = c("chrom", "exonStart", "exonStop","GeneName"), stringsAsFactors = F) %>%
  mutate(nonROH = "1") %>%
  right_join(., ensemblGenes) %>%
  mutate(data = "PLINK",
         nonROH = ifelse(is.na(nonROH), as.numeric(0), nonROH)) %>%
  rbind.data.frame(regionNonOverlaps_V)

#remove genes that are covered by roh
allGenesInfo = left_join(regionNonOverlaps, ExonSpace) %>%
  mutate(propGene = ifelse(nonROH == "1", round(exonLength/TSSLength$TSSLength[match(GeneName, TSSLength$GeneName)], digits = 1), as.numeric(0)),
         distTSS = lag(spaceBtwn) + spaceBtwn,
         propBasedOnTSS = ifelse(is.na(distTSS), "0", round(distTSS/TSSLength$TSSLength[match(GeneName, TSSLength$GeneName)], digits = 1))) 

NoROH = allGenesInfo %>%
  group_by(chrom,GeneName, data) %>%
  summarise(coverage = sum(propGene)) %>%
  filter(coverage > 0) 

plotDF = allGenesInfo %>%
  left_join(.,NoROH) %>%
  na.omit(coverage) %>%
  mutate(coverage = ifelse(coverage > 1, as.numeric(1), coverage))
  

#plot results

splitOnData = function(dataFrame, rohCaller){
  inputDF = dataFrame %>%
    filter(data == rohCaller)
  ggplot(inputDF, aes(x=exonStart, y=nonROH)) +
  geom_point() +
  labs(y = "nonROH status", x = "exon location") +
  facet_wrap(~GeneName, scales = "free", nrow = 3) +  
  theme_bw()+ 
  theme(axis.text.x = element_text( hjust= 0.5, vjust=0.75, size=14, angle = 40), 
        axis.text.y = element_text(size =20), 
        plot.title=element_text(size = 24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=24),
        strip.text = element_text(size=14)) 
}

plinkNonROH = splitOnData(plotDF, "PLINK")
vcfToolsNonROH = splitOnData(plotDF, "VCFTools")


ggplot(plotDF, aes(x=data,y=coverage)) +
  geom_violin() +
  coord_flip() +
  theme_bw() + 
  theme(axis.text.x = element_text( hjust= 0.5, vjust=0.75, size=14, angle = 40), 
        axis.text.y = element_text(size =20), 
        plot.title=element_text(size = 24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=24),
        strip.text = element_text(size=14)) 
