#Load Libraries
library(tidyverse)

#read file in and provide a start position
chroms = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/chromosomeLengths.txt", stringsAsFactors = F, sep = " ") %>%
  mutate(startPOS = as.integer(0))

#Break chroms up into (50Kb) steps
dataList = list()
for(i in 1:dim(chroms)[1]){
  ChromStepSize = seq.int(from = chroms[i,3], to = chroms[i,2], by = 5e6) %>%
    as.data.frame() %>%
    dplyr::rename("endPos" = ".") %>%
    mutate(startPos = endPos - 5e6,
           chrom = chroms[i,1]) %>%
    filter(startPos >= 0)
  dataList[[i]] = ChromStepSize
}
newGeneSet = bind_rows(dataList) %>%
  select(chrom, startPos, endPos) %>%
  mutate(startPos = as.integer(startPos),
         endPos = as.integer(endPos))

write.table(newGeneSet, file = "~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/Autosome_5MbWindows.bed", col.names = F, row.names = F, quote = F, sep = "\t")
