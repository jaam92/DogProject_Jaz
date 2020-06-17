#library
library(tidyverse)
library(gridExtra)

#fxn to negate in
`%nin%` = Negate(`%in%`)

#load files
ccrs = read_delim("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/VettingResults/ccrs.autosomes.v2.20180420.bed", delim = "\t")

t10_ccr = ccrs %>%
  filter(ccr_pct >= 90.0) %>%
  distinct(gene, .keep_all = T) #Load ccrs

allGenes = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/VettingResults/EnsemblGenes_CanFam3.1_SingleTranscript_tsStart_tsEnd_perGene.bed", stringsAsFactors = F) #Load all genes

nonROHexons = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/VettingResults/ExonRegion_NonOverlapsROH_cornellData.bed", stringsAsFactors = F) %>%
  distinct(V4, .keep_all = T) #Load file for exons without ROHs add hgnc info

#Fisher's exact test and hypergeomtric LRT
#make contingency table
annots = allGenes %>%
  select(V4) %>%
  mutate(isCCR = ifelse(V4 %in% t10_ccr$gene, as.numeric(1), as.numeric(0)),
         inROH = ifelse(V4 %nin% nonROHexons$V4, as.numeric(1), as.numeric(0))) %>%
  count(isCCR, inROH) %>%
  ungroup()

annots

#Fishers-Exact
counts = c(23,9896,4,5047)
contingencyTable = matrix(counts, nrow =2, ncol = 2)
test = fisher.test(contingencyTable)
sprintf("p-value from fisher's exact test %f", test$p.value)

#Likelihood with hypergeometric
ROH_size = annots %>% filter(inROH == 1) %>% summarise(n = sum(n)) %>% as.numeric()
NonROH_size = annots %>% filter(inROH == 0) %>% summarise(n = sum(n)) %>% as.numeric()
ROH_CCR_count = annots %>% filter(inROH == 1 & isCCR == 1) %>% select(n) %>% as.numeric()
NonROH_CCR_count = annots %>% filter(inROH == 0 & isCCR == 1) %>% select(n) %>% as.numeric()

#Loglikelihood Statistic
p = seq(0.01, 0.99, by = 0.001) #vector of p-values
NonROH = dbinom(NonROH_CCR_count, size = NonROH_size, prob = p, log = TRUE)
ROH = dbinom(ROH_CCR_count, size = ROH_size, prob = p, log = TRUE)

#Restricted Model log likelihood
restrictedTotal = NonROH + ROH
LLrestricted = max(restrictedTotal) #add the two vectors then identify a maximum because we are forcing maximum to be at the same p-value for both models

#p-value at the position where the max loglikelihood is
p[which(restrictedTotal == max(restrictedTotal))]

#Loglikelihood for Full Model 
LLfull = max(NonROH) + max(ROH) #add the maximums because we are assuming there are two different p-values (nonROH and ROH)

#p-value at the position where the max loglikelihood for nonROH and ROH are
p[which(NonROH == max(NonROH))]
p[which(ROH == max(ROH))]

#Likelihood ratio Statistic of the two models
LRStatistic = 2*(LLfull - LLrestricted)

#pvalue of the likelihood ratio under chi-square distribution with 1 df
p_of_LRStatistic = 1-pchisq(LRStatistic,  1)
sprintf("p-value from LLR test exact test %f", p_of_LRStatistic)

#If I sample 27 genes what is the probability that 23 are ccrs in the top 10%
set.seed(2020)
numReps = 100000
histReps = replicate(numReps, sum(sample(allGenes$V4, 27, replace = F) %in% t10_ccr$gene))
resampPval = sum(histReps >= 23)/numReps #p-value 
sprintf("p-value for resampling is %f", resampPval)

#plot results
enrichmentNonROH = ggplot() +
  geom_histogram(aes(x=histReps), bins=25, binwidth = 0.5) +
  geom_vline(aes(xintercept = 23), colour="blue") +
  labs(x="Number of CCR Genes", y="Count Replicates") +
  annotate("text", x=25, y=15000, label= paste0("p = ",resampPval), size = 14) +
  theme_bw() + 
  theme(plot.title=element_text(size=18, face = "bold", hjust=0.5),
        axis.title=element_text(size=16),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))

#plot with contingency table
FlipCounts = c(9896,23,5047,4)
FlipContingencyTable = matrix(FCounts, nrow =2, ncol = 2)
plotContTable = tableGrob(FlipContingencyTable, rows = c("ROH", "non-ROH"), cols =c("CCR", "non-CCR"), theme = ttheme(rownames.style = colnames_style(color = "black",face = "bold",size = 12, fill = "grey80",linewidth = 1,linecolor = "white")))
enrichmentNonROH + 
  annotation_custom(plotContTable, xmin = 7, xmax = 10, ymin = 12000, ymax = 15000)
#Now am I more likely to see an ROH overlapping a CCR in the bottom 20% or top 10%? We do this by first making sure none of the genes in the bottom 20% overlap any genes in the top 10% so if I annotate a gene as being a top CCR that is a priority over a bottom ccr then we downsample top ccrs to match bottom ccrs

b20_ccr = ccrs %>%
  filter(ccr_pct <= 20.0 & gene %nin% t10_ccr$gene) %>% 
  distinct(gene, .keep_all = T)

#label genes that rohs fall in the top 10% and bottom 10% of ccrs
rohsInGenes = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/VettingResults/FullGeneRegion_OverlapsROH_cornellData.bed", stringsAsFactors = F) %>%
  mutate(t10_CCR = ifelse(GeneName %in% t10_ccr$gene, as.numeric(1), as.numeric(0)),
         b20_CCR = ifelse(GeneName %in% b20_ccr$gene, as.numeric(1), as.numeric(0))) 

#Downsample the top 10% rohs and count the number of ROH overlapping those genes
downsampleROHOverlaps = function(sampSize){
  numOverlaps = rohsInGenes %>%
    filter(t10_CCR == "1") %>%
    sample_n(sampSize, replace = F) %>%
    summarise_at(c("NumROHOverlap"), sum, na.rm=TRUE) %>%
    as.numeric()    
  return(numOverlaps)
}

#repeat the downsampling 
overlaps_b20 = rohsInGenes %>%
  filter(b20_CCR == "1") %>%
  summarise_at(c("NumROHOverlap"), sum, na.rm=TRUE) %>%
  as.numeric()
set.seed(48)
histReps = replicate(numReps,downsampleROHOverlaps(3579))
resampPval = sum(histReps >= overlaps_b20)/numReps #p-value 
sprintf("p-value for resampling is %f", resampPval)

#plot the results
downsampleT10 = ggplot() +
  geom_histogram(aes(x=histReps), bins = 50) +
  geom_vline(aes(xintercept = overlaps_b20), colour="blue") +
  labs(x="Number of ROHs overlapping CCR Genes", y="Count Replicates") +
  annotate("text", x=3335500, y=450, label = paste0("p = ", resampPval),  size = 14) +
  theme_bw() + 
  theme(plot.title=element_text(size=18, face = "bold", hjust=0.5),
        axis.title=element_text(size=16),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))

###Incorporating Human Results from Mooney et al. 2018 AJHG
humanNonROHexons = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/VettingResults/ExonRegion_NonOverlapsROH_Mooney.bed", stringsAsFactors = F) %>%
  distinct(V4, .keep_all = T)

humanExons = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/VettingResults/protein_coding_genes_hg19_HGNC.bed", stringsAsFactors = F)

#make contingency table
HumAnnots = allGenes %>%
  select(V4) %>%
  mutate(isCCR = ifelse(V4 %in% t10_ccr$gene, as.numeric(1), as.numeric(0)),
         inROH = ifelse(V4 %nin% humanNonROHexons$V4, as.numeric(1), as.numeric(0))) %>%
  count(isCCR, inROH) %>%
  ungroup()

HumAnnots

#If I sample 27 genes what is the probability that 3981 are ccrs in the top 10%
set.seed(2021)
numReps = 100000
histReps = replicate(numReps, sum(sample(humanExons$V4, 3981, replace = F) %in% t10_ccr$gene))
resampPval = sum(histReps >= 2752)/numReps #p-value 
sprintf("p-value for resampling is %f", resampPval)

humData = ggplot() +
  geom_histogram(aes(x=histReps), bins=50, binwidth = 1) +
  geom_vline(aes(xintercept = 2752), colour="blue") +
  labs(x="Number of CCR Genes", y="Count Replicates") +
  annotate("text", x=2500, y=1200, label = deparse(bquote(~p < 1~e^-05)),  size = 14, parse=T) +    
  theme_bw() + 
  theme(plot.title=element_text(size=18, face = "bold", hjust=0.5),
        axis.title=element_text(size=16),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))
print(humData)


