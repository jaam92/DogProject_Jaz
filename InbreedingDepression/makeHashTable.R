#Load Libraries 
#suppressPackageStartupMessages(library(tidyverse))

#Read autosome length and roh files in 
#autosome = read.delim(file = "../chromosomeLengths.txt", check.names = F, stringsAsFactors = F, sep = " ")  %>%
#  mutate(end = cumsum(as.numeric(LENGTH)),
#         start = lag(end, default = 0),
#         check = end - start) 


#write to files
#for(i in 1:38){
#  hash = cbind.data.frame(c(0:autosome[i,2]), c(autosome[i,4]:autosome[i,3])) %>%
#    mutate(chrom = paste0("chr",i)) 
  
#  write.table(hash, file = paste0("hashTableChrom", i, ".txt"), quote = F, sep = "\t", row.names = F, col.names = F)
#}


#concatenate files 
#system("for f in {1..38}; do cat hashTableChrom$f.txt >> hashTableChrom_allChroms.txt;done")

#gzip files
#system("gzip hashTableChrom_allChroms.txt")
