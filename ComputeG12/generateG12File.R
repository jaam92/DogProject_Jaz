#Load Libraries
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyverse)

#Load files with population size greater than 50 
setwd("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/ComputeG12")
pops = c("boxer", "cocker_spaniel", "german_shepherd_dog", "golden_retriever", "grayWolf_Europe", "grayWolf_NorthAmerica", "labrador_retriever", "maltese", "mix", "newfoundland", "poodle", "rottweiler", "village_dog_peru", "yorkshire_terrier")

for (i in pops){
 
  indivs = read.table(file = paste("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/IndividualFiles/UnrelatedByBreed/" ,i, ".txt", sep = ""))

  ids = indivs$V1
 
  for (chrom in 1:38){
  	vcf = fread(file=paste("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/splitVCF/MergedFile_CornellCanineFitak_UnrelatedsOnly_chr", chrom, ".vcf", sep=""),check.names=F) #remove the null column that get's added
    
    #Subset populations of interest and recode
    
    vcfCols = colnames(vcf)
    siteInfo = c(2,4,5) #POS,REF,ALT
    ColsToKeep = append(siteInfo,which(vcfCols %in% ids)) #keep info and individuals of interest
    popSpecificDF = vcf[, ColsToKeep, with = F] #subset based on vector of colnames 
    numCols = length(ColsToKeep) #current column number
    recodedVCF = popSpecificDF %>% 
      mutate_at(.vars=4:numCols, funs(case_when(. == "0/0" ~ REF, 
      . == "1/1" ~ ALT,
      . == "./." ~ "N",
      . == "0/1" ~ ".",
      . == "1/0" ~ "."))) %>%
      select(-c(REF, ALT)) 
    rm(popSpecificDF)

    #Remove fixed sites
    
    end = dim(recodedVCF)[2]
    numSamps = end - 1 
    keep = apply(recodedVCF[2:end], 1, function(x) length(unique(x[!is.na(x)])) > 1)

    #Uncomment to remove sites based on amount of missing data (max missing data is 10%)

    threshold = ceiling(numSamps*.10)
    intermidiateDF = recodedVCF[keep,] 
    intermidiateDF$NumMissing = rowSums(intermidiateDF == "N")
    FinalDF = intermidiateDF %>%
      filter(NumMissing < threshold) %>%
      select(-c(NumMissing))
    
    write.table(FinalDF, file = paste(i,"_n", numSamps,"_Chr", chrom ,"_inputG12haps.txt", sep = ""), sep = "," , quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  
}
    

