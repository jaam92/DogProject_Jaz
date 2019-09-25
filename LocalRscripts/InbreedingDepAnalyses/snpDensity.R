#Load Libraries
library(tidyverse)
library(GenomicRanges)
library(ggplot2)

####Read files in
genes = read.delim("~/DogProject_Jaz/LocalRscripts/CaseControlROH/EnsemblGenes_CanFam3.1.bed", stringsAsFactors = F)
gene_names = read.table("~/DogProject_Jaz/LocalRscripts/CaseControlROH/EnsemblGenes_CanFam3.1_geneNames.txt", stringsAsFactors = F)
bim = read.table("~/DogProject_Jaz/LocalRscripts/PCA_Unrelateds/MergedFile_CornellCanineFitak_allIndivs.bim", stringsAsFactors = F)

####Make file with genes of interest
####Filter to only chromosomes 1-38 
####keep only longest transcript (need to use distinct too bc some transcripts have equal length and we only want to keep one entry)

#gene set for Ensembl
GeneSet = genes %>% 
  mutate(chrom = gsub("chr", "", chrom), 
         AbbrevName = gene_names$V2[match(name, gene_names$V1)], 
         transcript_length = transcriptionEnd - transcriptionStart) %>% 
  group_by(AbbrevName) %>% 
  filter(transcript_length == max(transcript_length) & as.numeric(chrom) <= 38) %>%
  distinct(AbbrevName,.keep_all= TRUE)  %>% 
  select(name, AbbrevName, chrom, transcriptionStart, transcriptionEnd) %>% 
  as.data.frame() %>%
  mutate(TS_50upstream = transcriptionStart - 50000,
         TS_50downstream = transcriptionEnd + 50000)

#Break genes up into (1Kb) steps
dataList = list()
for(i in 1:dim(GeneSet)[1]){
  GeneStepSize = seq.int(from = GeneSet[i,6], to = GeneSet[i,7], by = 1000) %>%
    as.data.frame() %>%
    dplyr::rename("endPos" = ".") %>%
    mutate(startPos = endPos - 1000,
           GeneName = GeneSet[i,2],
           chrom = GeneSet[i,3]) 
  dataList[[i]] = GeneStepSize
}
newGeneSet = bind_rows(dataList)

#Make GRanges object to find overlaps
posRanges = with(bim, GRanges(V1, IRanges(start=V4, end = V4), startPos = V4))
geneRanges = with(newGeneSet, GRanges(chrom, IRanges(start=startPos, end = endPos), id = GeneName, startPos = startPos))

#Find Positions Overlapping Genes
SNP_Gene_Overlaps = findOverlaps(query = posRanges, subject = geneRanges,  type = "within")

#Get info from overlaps and make a data frame
my_query = queryHits(SNP_Gene_Overlaps)
my_subject = subjectHits(SNP_Gene_Overlaps)

snpDensityPerGene = data.frame(GeneName = geneRanges[my_subject]$id, startPos = geneRanges[my_subject]$startPos, snpStart = posRanges[my_query]$startPos) 

x = merge.data.frame(newGeneSet, snpDensityPerGene, by = c("GeneName", "startPos"), all = T) %>%
  mutate(snpPresent = ifelse(is.na(snpStart), "0","1")) %>%
  group_by(GeneName, startPos) %>%
  mutate(countSNPs = sum(is.numeric(snpPresent))) %>%
  ungroup()
