#Load library 
library(tidyverse)

#fxn to negate in
`%nin%` = Negate(`%in%`)

#fxn to downsample the top 10% rohs and count the number of ROH overlapping those genes
downsample_t10 = function(sampSize){
  numOverlaps = hetWindowed %>%
    filter(t10_CCR == "1") %>%
    sample_n(sampSize, replace = F) %>%
    summarise_at(c("heterozygosity"), sum, na.rm=TRUE) %>%
    as.numeric()    
  return(numOverlaps)
}

#Load ccr file and pull ccrs of interest
ccrs = read_delim("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/VettingResults/ccrs.autosomes.v2.20180420.bed", delim = "\t")

t10_ccr = ccrs %>%
  filter(ccr_pct >= 90.0) %>%
  distinct(gene, .keep_all = T) 

b20_ccr = ccrs %>%
  filter(ccr_pct <= 20.0 & gene %nin% t10_ccr$gene) %>% 
  distinct(gene, .keep_all = T)

#Load windowed heterozygosity and annotate
hetWindowed = read_delim("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/VettingResults/MidpointGeneCoord_heterozygosity250KbWindowEachSide.bed", delim = "\t", col_names = c("chrom","midpointS", "midpointE","geneName","chromBim", "posS","posE","heterozygosity")) %>%
  group_by(geneName) %>%
  summarise_at(.vars= vars(heterozygosity), .funs = mean) %>%
  ungroup() %>%
  mutate(t10_CCR = ifelse(geneName %in% t10_ccr$gene, as.numeric(1), as.numeric(0)),
         b20_CCR = ifelse(geneName %in% b20_ccr$gene, as.numeric(1), as.numeric(0)))

#repeat the downsampling 
count_b20 = sum(ifelse(hetWindowed$b20_CCR == 1, as.numeric(1), as.numeric(0)))
overlaps_b20 = hetWindowed %>%
  filter(b20_CCR == "1") %>%
  summarise_at(c("heterozygosity"), sum, na.rm=TRUE) %>%
  as.numeric()
set.seed(2020)
numReps = 100000
histReps = replicate(numReps,downsample_t10(count_b20))
resampPval = sum(histReps >= overlaps_b20)/numReps #p-value 
sprintf("p-value for resampling is %f", resampPval)

#plotting
histogramHetResamp = ggplot() +
  geom_histogram(aes(x=histReps), bins = 50) +
  geom_vline(aes(xintercept = overlaps_b20), colour="blue") +
  labs(x="Distribution of heterozygosity from midpoint of top 10% CCR Genes", y="Count Replicates") +
  theme_bw() + 
  theme(plot.title=element_text(size=18, face = "bold", hjust=0.5),
        axis.title=element_text(size=16),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))

distHet = ggplot(hetWindowed, aes(x=heterozygosity)) + 
  geom_histogram(bins=25) + 
  labs(x = "Heterozygosity (250Kb window)", title = "Distribution of heterozygosity from midpoint of gene") +
  theme_bw() + 
  theme(plot.title=element_text(size=18, face = "bold", hjust=0.5),
        axis.title=element_text(size=16),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))
