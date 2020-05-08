#Load files for hgnc symbols and ccr data 
hgnc_symbols = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/VettingResults/hgnc_symbol_pairs.txt", col.names = c("previous", "current"), stringsAsFactors = F)

ccrs = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/VettingResults/ccrs.autosomes.90orhigher.v2.20180420.bed", stringsAsFactors = F, check.names = F) %>%
  distinct(gene, .keep_all = T)

#Load all genes
allGenes = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/ForAbi_EnsemblGenes_CanFam3.1_SingleTranscript.bed", stringsAsFactors = F) %>%
  distinct(V4, .keep_all = T) %>%
  left_join(hgnc_symbols, by = c("V4"="previous"))

#Load file for exons without ROHs add hgnc info
nonROHexons = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/VettingResults/ExonRegion_NonOverlapsROH_cornellData.bed", stringsAsFactors = F) %>% 
  left_join(hgnc_symbols, by = c("V4"="previous")) %>% #merge with hgnc symbols
  distinct(V4, .keep_all = T) #keep one copy of each gene entry

#If I sample 27 genes what is the probability that 23 are ccrs
set.seed(2020)
numReps = 500000
histReps = replicate(numReps, sum(sample(allGenes$current, 27, replace = F) %in% ccrs$gene))
resampPval = sum(histReps >= 23)/numReps #p-value 
sprintf("p-value for resampling is %f", resampPval)

#plot results
ggplot() +
  geom_histogram(aes(x=histReps), bins=20, binwidth = 0.5) +
  geom_vline(aes(xintercept = 23), colour="blue") +
  labs(x="Number of CCR Genes", y="Count Replicates") +
  annotate("text", x=24, y=55000, label= paste0("p-value=",resampPval)) +    
  theme_bw() + 
  theme(plot.title=element_text(size=18, face = "bold", hjust=0.5),
        axis.title=element_text(size=16),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))

#Fisher's exact test and hypergeomtric LRT
#make contingency table
annots = allGenes %>%
  select(current) %>%
  mutate(isCCR = ifelse(current %in% ccrs$gene, as.numeric(1), as.numeric(0)),
         inROH = ifelse(current %in% nonROHexons$current, as.numeric(0), as.numeric(1))) %>%
  count(isCCR, inROH) %>%
  ungroup()

annots

#Fishers-Exact
counts = c(23,9935,4,5079)
contingencyTable = matrix(counts, nrow =2, ncol = 2)
test = fisher.test(contingencyTable)
sprintf("p-value from fisher's exact test %f", test$p.value)

#Likelihood with hypergeometric
ROH_size = annots %>% filter(inROH == 1) %>% summarise(n = sum(n))
NonROH_size = annots %>% filter(inROH == 0) %>% summarise(n = sum(n))
ROH_CCR_count = annots %>% filter(inROH == 1 & isCCR == 1) %>% select(n)
NonROH_CCR_count = annots %>% filter(inROH == 0 & isCCR == 1) %>% select(n)

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



