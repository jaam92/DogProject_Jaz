#load libraries 
library("tidyverse")

#set working directory and load files
setwd("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/")

ensemblGenes = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/ForAbi_EnsemblGenes_CanFam3.1_SingleTranscript.bed", col.names = c("chrom", "exonStart", "exonStop","GeneName"), stringsAsFactors = F)


regionNonOverlaps = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/ExonRegion_NonOverlapsROH_cornellData.bed",  col.names = c("chrom", "exonStart", "exonStop","GeneName"), stringsAsFactors = F) %>%
        mutate(nonROH = "1")

#gives the space in between exons within a gene
ExonSpace = ensemblGenes %>%
  group_by(GeneName) %>%
  mutate(diff = exonStart - lag(exonStop, default = first(exonStart))) %>%
  as.data.frame()

#merge data frames to figure out which exons do not have ROH
nonROHExonLocation = ensemblGenes %>%
  filter(GeneName%in%regionNonOverlaps$GeneName) %>%
  left_join(., regionNonOverlaps) 

nonROHExonLocation[is.na(nonROHExonLocation)] <- 0 #replace NA with 0


#Find length from transcription start site
TSSLength = merge(nonROHExonLocation, ExonSpace) %>% 
        group_by(chrom, GeneName) %>% 
        summarise(TSSLength = sum(diff)) 

TSSInfo = merge(nonROHExonLocation, ExonSpace) %>% 
        group_by(chrom, GeneName) %>%
        slice(which.min(exonStart)) %>%
        select(chrom, exonStart, GeneName) %>%
        mutate(TSSLength = TSSLength$TSSLength[match(GeneName, TSSLength$GeneName)]) %>%
        ungroup()

#plotting location of nonROH exons
plotDF = left_join(nonROHExonLocation, ExonSpace) %>%
  mutate(TSSInfo = TSSInfo$TSSLength[match(GeneName, TSSInfo$GeneName)],
         distTSS = lag(diff) + diff,
         propGene = ifelse(is.na(distTSS), "0", round(distTSS/TSSLength$TSSLength[match(GeneName, TSSLength$GeneName)], digits = 2)))

ggplot(plotDF, aes(x=exonStart, y=nonROH)) +
  geom_point() +
  labs(y = "nonROH status", x = "exon location") +
  facet_wrap(~GeneName, scales = "free", nrow = 3) +  
  theme_bw()+ 
  theme(axis.text.x = element_text( hjust= 0.5, vjust=0.75, size=14, angle = 40), 
        axis.text.y = element_text(size =20), 
        plot.title=element_text(size = 24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=24),
        strip.text = element_text(size=14)) 


#New plots
nonROHGenes = unique(nonROHExonLocation$GeneName)

genesWithinOneMb <- function(geneOfInterest){
  coord = fullGenomeAnnot %>%
    filter(GeneName == "CENPS")
  start = min(coord$exonStart) - 500000
  end = min(coord$exonStop) + 500000
  chromo = unique(coord$chrom)
  genomicContext = fullGenomeAnnot %>%
    filter(chrom%in%chromo & exonStart >= start &  exonStop <= end) %>%
    arrange(exonStart)
}
