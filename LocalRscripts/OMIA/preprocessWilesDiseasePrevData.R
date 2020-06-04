#####Set Directory
setwd("~/Documents/DogProject_Jaz/LocalRscripts/OMIA")

#####Load Libraries
library(tidyverse)

####Read input file with prevalence data from Wiles et. al 2017 Canine Evolutionary Genomics

preprocess = read.delim("WilesData_DiseasePrev_2017CEG.txt", check.names = F, stringsAsFactors = F)

####Process data frame
  #first address breed names so they match data
  #then process data frame
    #map function will pull the prevalence point estimate
    #then add trait names back
colnames(preprocess) = colnames(preprocess) %>% 
  str_replace(" \\(.*\\)", "") %>%
  str_to_lower() %>% 
  str_replace_all(fixed(" "), "_")

colnames(preprocess)[1] = "trait" #give the first column a name
colnames(preprocess)[18] = "flatcoated_retriever" #rename the flat-coated retriever
colnames(preprocess)[30] = "poodle_miniature" #rename the first instance of poodle as mini the second instance is standard

postProcess = preprocess[,-1] %>% 
  map_dfr(~ parse_number(.x,"[[:digit:]]+")) %>%
  mutate(trait = preprocess$trait) %>%
  select(trait, everything()) %>%
  as.data.frame()

#Save processed data
#write.table(postProcess, file = "~/Documents/DogProject_Jaz/LocalRscripts/OMIA/AKC_DiseasePrev_PointEstWiles2017.txt", quote = F, row.names = F, sep = "t")
