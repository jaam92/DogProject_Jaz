#Load Libraries
library(tidyverse)

#read file in and provide a start position
chroms = read.delim("~/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/chromosomeLengths.txt", stringsAsFactors = F, sep = " ") 

#Break chroms up into (50Kb) steps
dataList = list()
for(i in 1:dim(chroms)[1]){
  ChromStepSize = seq.int(from = 0, to = chroms[i,2], by = 50000) %>%
    as.data.frame() %>%
    dplyr::rename("endPos" = ".") %>%
    mutate(startPos = endPos - 50000,
           chrom = chroms[i,1]) 
  dataList[[i]] = ChromStepSize
}
newGeneSet = bind_rows(dataList) %>%
  select(chrom, startPos, endPos)

#write.table(newGeneSet, file = "~/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/Autosome_50KbWindows.bed", col.names = F, row.names = F, quote = F, sep = " ", )