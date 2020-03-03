#Load Libraries 
library(tidyverse)
library(tictoc)

#Read autosome length file in 
autosome = read.delim(file = "~/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/chromosomeLengths.txt", check.names = F, stringsAsFactors = F, sep = " ") 

#Fxn for rounding this rounds down 
round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}

#Fxn to turn each chromosome into a vector of 0 or number aka nonROH/ROH and then permute the location of ROH with samp
permuteROH = function(dataFrame, chromosomeNumber){
    #add up all the blocks between ROH and then add end of final ROH to end of chromosome 
    totalNonROH = sum(dataFrame$distBetween) + (autosome[chromosomeNumber,2] - tail(dataFrame, n=1)$AUTO_END)
    #calculate the total number of 100 bp nonROH blocks by rounding totalNonROH to the nearest 100 then dividing by 100
    total_nonROHBlocks = round_any(x = totalNonROH, accuracy = 100) / 100
    #make a vector of the chromosome where 0 is each 100bp nonROH block and ROHs are represented as their lengths
    recodeChrom = c(rep(0, total_nonROHBlocks), dataFrame$AUTO_LEN)
    #permute rohs with sample  
      #recode with rle format and make a data frame
      #convert data frame to chromosome coords
        #multiply number of 0 by 100 since we spaced every 100 bp else leave alone
        #start is where previous element ends and the end is cummulative sum minus 1 because our default on lag is 0-based
      #output roh only
    permuted = data.frame(unclass(rle(sample(recodeChrom))))  %>%
      mutate(lengths = ifelse(values == "0", lengths*100, values),
             start = lag(lengths, default = 0),
             end = cumsum(lengths) - 1) %>% 
      filter(values != "0") %>%
      mutate(chrom = chromosomeNumber) %>%
      select(chrom, start, end) 
    return(permuted)

}

#Permute ROH for all individuals and on all chroms
tic("total time all indivs and chroms")
for(i in 1:38){
  print(i)
  
  mylist = read.delim(file = "~/DogProject_Jaz/LocalRscripts/ROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt", stringsAsFactors = F) %>%
    filter(CHROM == i) %>%
    arrange(INDV, AUTO_START) %>%
    select(INDV, CHROM, AUTO_START, AUTO_END, AUTO_LEN) %>%
    group_by(INDV, CHROM) %>% 
    mutate(distBetween = AUTO_START - lag(AUTO_END, default = 0)) %>%
    ungroup() %>%
    group_by(INDV) %>%
    group_map(~ permuteROH(dataFrame = .x, chromosomeNumber = i))
  
  lapply(mylist, function(x) write.table(x, "~/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/test.bed", quote = F, sep = "\t", row.names = F, col.names = F, append = T))
}


toc()
