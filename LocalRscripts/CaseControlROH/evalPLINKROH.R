library(IRanges)
library(data.table)
library(dplyr)
library(tidyr)

#Load Files
IW_pheno = read.delim("~/Documents/DogProject_Jaz/ROH/IrishWolfhounds_Epilepsy.txt")
ROH = read.delim("~/Documents/DogProject_Jaz/ROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt")

#separate ROH by case and control
IW_control = IW_pheno %>% filter(epilepsy_irishWolfhounds == 1)
IW_case = IW_pheno %>% filter(epilepsy_irishWolfhounds == 2)

#Make union of ROH
DT = as.data.table(ROH)

## find interval for each chromosome
DT[,group := { 
  ir =  IRanges(AUTO_START, AUTO_END);
  subjectHits(findOverlaps(ir, reduce(ir)))
},by=CHROM]

## Now I group by group and chrom 
UnionROH = DT[, list(AUTO_START=min(AUTO_START),AUTO_END=max(AUTO_END),INDV=list(INDV),CHROM=unique(CHROM)),
              by=list(group,CHROM)]

#Unnest split ROH
names(UnionROH)[2] = "chrom"
sepINDVconsensus = UnionROH %>% select(CHROM, AUTO_START,AUTO_END,INDV) %>% mutate(INDV= strsplit(as.character(INDV), ",")) %>% unnest(INDV) 

#Find case control
IW_control_ROH = sepINDVconsensus %>% filter(INDV %in% IW_control$dogID) %>% group_by(AUTO_START,AUTO_END) %>% count() %>% rename(countControl = n) %>% as.data.frame() 

IW_case_ROH = ROH %>% filter(FID %in% IW_case$dogID) %>% group_by(POOL) %>% count() %>% rename(countCase = n) %>% as.data.frame() 

