#load libraries 
library("dplyr")
library("ggplot2")

#set working directory and load files
setwd("~/Desktop/abidata")
ensemblGenes = read.table("ForAbi_EnsemblGenes_CanFam3.1_SingleTranscript.bed", col.names = c("chrom", "exonStart", "exonStop","GeneName"), stringsAsFactors = F)

regionNonOverlaps = read.table("ExonRegion_NonOverlapsROH.bed",  col.names = c("chrom", "exonStart", "exonStop","GeneName"), stringsAsFactors = F) %>%
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
ggplot(nonROHExonLocation, aes(x = exonStart, y = nonROH)) +
        geom_point() + 
        facet_wrap(~ GeneName, scales = "free") +
        theme_bw() +
        theme(axis.text.x = element_blank())
