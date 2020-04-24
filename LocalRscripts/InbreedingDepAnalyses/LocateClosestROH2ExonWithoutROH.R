#Load library
library(tidyverse)

#fxn to pull rows of interest by finding a keyword (exon name) and returning 1 row above and below it
findRowsOfInterest <- function(dataFrame, keyword){
  rn <- outer(which(dataFrame$INDV == keyword), -1:1, `+`) %>% 
    as.vector() %>% 
    unique() %>% 
    Filter(function(x) x[1 <= x & x <= nrow(dataFrame)], .)
  filteredDF = dataFrame[rn,] %>%
    group_by(CHROM) %>% #this will make sure we keep all on same chrom in case the next closest value is on a different chrom
    mutate(diff = AUTO_END - lag(AUTO_END, default = first(AUTO_END)))
  return(filteredDF)
}

#Load files
exons = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/ExonRegion_NonOverlapsROH_cornellData.bed") %>% 
  #filter(V4 == "ANKH"| V4 =="FYTTD1"| V4 =="PRMT2" ) %>% #pull exons of interest
  mutate(chrom = as.numeric(gsub("chr", "", V1))) %>% 
  select(V2,V3,V4,chrom) %>% 
  rename(AUTO_START = V2,AUTO_END = V3,INDV = V4,CHROM = chrom) 

lengths = exons %>% 
  mutate(AUTO_LEN = AUTO_END - AUTO_START) %>% group_by(INDV) %>% 
  summarise(lengthGene = sum(AUTO_LEN))

rohs = read.delim(file = "~/Documents/DogProject_Jaz/LocalRscripts/ROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFile_Cornell_allChroms_vcfToolsROH_rmROHlessThan50snps_HaywardDataOnly.txt")

#merge and arrange by start and end
mergedDF = rohs %>%  
  select(AUTO_START,AUTO_END,INDV,CHROM) %>% 
  rbind.data.frame(exons) %>% 
  arrange(CHROM, AUTO_END)

#Grab rows of interest
ANKH = findRowsOfInterest(mergedDF, "ANKH")
PRMT2 = findRowsOfInterest(mergedDF, "PRMT2")
FYTTD1 = findRowsOfInterest(mergedDF, "FYTTD1")
