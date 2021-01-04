#Load Libraries
library(tidyverse)
library(GenomicRanges)
library(ggplot2)

####Read files in
genes = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/EnsemblGenes_CanFam3.1.bed")
gene_names = read.table("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/EnsemblGenes_CanFam3.1_geneNames.txt")
gene_names_HGNC = read.table("~/Documents/HIVProject/PopBranchStat/hgnc_symbol_pairs.txt")
ROH = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/ROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt")
goldenRetrievers = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/splitPhenotypeFile/IncludeMixedBreeds/lymphoma_golden_retriever.txt")
lymphomaGenes = read.table("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/Skibola2010_LymphomaAssocGenes.txt")

####Make file with genes of interest
####Filter to only chromosomes 1-38 
####keep only longest transcript (need to use distinct too bc some transcripts have equal length and we only want to keep one entry)

#gene set for Ensembl
GeneSet = genes %>% 
  mutate(chrom = gsub("chr", "", chrom), AbbrevName = gene_names$V2[match(name, gene_names$V1)], transcript_length = transcriptionEnd - transcriptionStart) %>% 
  group_by(AbbrevName) %>% 
  filter(transcript_length == max(transcript_length) & as.numeric(chrom) <= 38) %>%
  distinct(AbbrevName,.keep_all= TRUE)  %>% 
  select(name, AbbrevName, chrom, transcriptionStart, transcriptionEnd) %>% 
  as.data.frame()

#ROH
GR_ROH = ROH %>% 
  filter(INDV %in% goldenRetrievers$dogID)

#Genes associated with lymphoma and mitral valve degeneration
lymphomaGenes_HGNC = lymphomaGenes %>% 
  mutate(HGNC_name = gene_names_HGNC$V2[match(V1,gene_names_HGNC$V1)]) 

#Overlap ROH with genes
gr0 = with(GeneSet, GRanges(chrom, IRanges(start=transcriptionStart, end = transcriptionEnd)))
gr1 = with(GR_ROH, GRanges(CHROM, IRanges(start=AUTO_START, end = AUTO_END)))
compeleteOverlaps_GR = findOverlaps(query = gr0, subject = gr1, type = "within")#complete overlap
ROH_compOverlaps_LymphGR = data.frame(GeneSet[queryHits(compeleteOverlaps_GR),], GR_ROH[subjectHits(compeleteOverlaps_GR),]) %>% 
  mutate(transcriptLength = as.numeric(transcriptionEnd) - as.numeric(transcriptionStart))

####Calculations Lymphoma
GR_AverageCountROHwithinGenes = ROH_compOverlaps_LymphGR %>% 
  group_by(INDV) %>% 
  tally() %>% #total number overlaps
  mutate(status = goldenRetrievers$lymphoma[match(INDV, goldenRetrievers$dogID)]) %>% #add case control label
  group_by(status) %>% 
  summarise(meanROHinGeneCount = mean(n)) #average count genes within ROH for case control

GR_AverageCountROHwithinLymphomaAssocGenes = ROH_compOverlaps_LymphGR %>% 
  filter(AbbrevName %in% lymphomaGenes_HGNC$HGNC_name) %>% #filter to lymphoma associated genes
  group_by(INDV) %>% 
  tally() %>% #total number overlaps
  mutate(status = goldenRetrievers$lymphoma[match(INDV, goldenRetrievers$dogID)]) %>% #add status 
  group_by(status) %>% 
  summarise(meanROHinGeneCount = mean(n))  #average count lymphoma associated genes in ROH for case control

ROHBurdenwithinLymphomaAssociatedGenes = ROH_compOverlaps_LymphGR %>% 
  filter(AbbrevName %in% lymphomaGenes_HGNC$HGNC_name) %>% 
  group_by(INDV) %>% 
  summarise(genesCoveredByROH = sum(as.numeric(transcriptLength))) %>% 
  mutate(status = goldenRetrievers$lymphoma[match(INDV, goldenRetrievers$dogID)],
         plottingID = ifelse(status == 2, "Case", "Control"))
  
####Run Permutation Test
PermutationTest = function(dataFrame, caseIndicator, controlIndicator, numberPerms){
  ROHBurden = dataFrame$genesCoveredByROH
  group = as.character(dataFrame$status)
  testStat = function(w, g) mean(as.numeric(w[g == caseIndicator])) - mean(as.numeric(w[g == controlIndicator])) #test whether there is a difference in median amount of genome in ROH btwn cases and controls
  observedStat = testStat(ROHBurden, group)
  permutations = sapply(1 : numberPerms, function(i) testStat(ROHBurden, sample(group)))
  pvalue = mean(permutations > observedStat) #pvalue 
  pvalue[pvalue == 0] = 1/numberPerms+1 #generate pvalue when there are no perms greater than observed
  bins = seq(1:numberPerms)
  dfPerms = cbind.data.frame(bins,permutations)
  return(list(dfPerms, observedStat, pvalue))
}

#Set seed and run permutation
set.seed(175)
GoldenRetrievers_Perms = PermutationTest(ROHBurdenwithinLymphomaAssociatedGenes, "2", "1", 10000)
GoldenRetrievers_Perms[[3]]

ggplot(ROHBurdenwithinLymphomaAssociatedGenes, aes(x=plottingID, y=genesCoveredByROH, group=plottingID)) + 
  geom_boxplot() + 
  scale_x_discrete(limits= c("Control", "Case")) + 
  theme_bw() 

ggplot(GoldenRetrievers_Perms[[1]], aes(GoldenRetrievers_Perms[[1]]$permutations)) + 
  geom_histogram(binwidth = 30, breaks=seq(-6.0e5, 6.0e5, by =10000),col="coral2", fill="white") + 
  geom_vline(xintercept = GoldenRetrievers_Perms[[2]], col="purple") +
  theme_bw() + 
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        plot.title=element_text(size=26, face = "bold"), 
        axis.title=element_text(size=24), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20)) +
  labs(x= "Permutation Score", 
       y= "Count",
       title = "Lymphoma Associated Genes in ROH Case vs Control (Golden Retriever)")


