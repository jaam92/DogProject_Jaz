#load libraries 
library("tidyverse")

#set working directory and load files
setwd("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/")
ensemblGenes = read.table("~/Documents/DogProject_Jaz/InbreedingDepression/ForAbi_EnsemblGenes_CanFam3.1_SingleTranscript.bed", col.names = c("chrom", "exonStart", "exonStop","GeneName"), stringsAsFactors = F)
genes = read.delim()
regionNonOverlaps = read.table("~/Documents/DogProject_Jaz/InbreedingDepression/vcftools/ExonRegion_NonOverlapsROH_vcfTools.bed",  col.names = c("chrom", "exonStart", "exonStop","GeneName"), stringsAsFactors = F) %>%
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


#plot the location of exons in roh or not
TSSLength = merge(nonROHExonLocation, ExonSpace) %>% 
        group_by(chrom, GeneName) %>% 
        summarise(TSSLength = sum(diff)) 
TSSInfo = merge(nonROHExonLocation, ExonSpace) %>% 
        group_by(chrom, GeneName) %>%
        slice(which.min(exonStart)) %>%
        select(chrom, exonStart, GeneName) %>%
        mutate(TSSLength = TSSLength$TSSLength[match(GeneName, TSSLength$GeneName)]) %>%
        ungroup()
#rm(TSSLength)

plotDF = merge(ExonSpace, TSSInfo)

plotDF$distTSS = lag(plotDF$diff) + plotDF$diff
plotDF$propGene = ifelse(is.na(plotDF$distTSS), "0", round(plotDF$distTSS/TSSLength$TSSLength[match(plotDF$GeneName, TSSLength$GeneName)], digits = 2))

ggplot(plotDF, aes(x = propGene, fill=GeneName)) +
        geom_density() + 
        theme_bw() +
        theme(legend.position = "none")
