#####Set Directory
setwd("~/Documents/DogProject_Jaz/LocalRscripts/OMIA")

#####Load Libraries
library(tidyverse)

#####Load Files and modify
omiaGenes = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/OMIA/omia-genes_v2.txt", fill = TRUE)

#FinalCausalVars = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/OMIA/causalVars_OMIA.txt", fill = TRUE) %>% 
FinalCausalVars = read.csv("C:/Users/Jazlyn/Downloads/variants.csv", stringsAsFactors = F) %>% 
  mutate(Breed.s.= strsplit(as.character(Breed.s.), ",")) %>% 
  unnest(Breed.s.)  %>% #split the breeds that are comma delimited into separate rows with same info
  select(OMIA.ID.s., Breed.s., Gene) %>%
  rename(OMIA_ID = OMIA.ID.s., Breed = Breed.s.) %>%
  mutate(NCBI_gene_ID = omiaGenes$ncbi_gene_id[match(Gene, omiaGenes$gene_symbol)], #add ncbi gene id
         OMIA_ID = gsub(" ", "_", OMIA_ID), #remove spaces in OMIA ID
         Breed = tolower(Breed), #lower case
         Breed = gsub("(?<=\\w)\\s(?=\\w)", "_", Breed, perl = TRUE), #add underscore between words
         Breed = gsub("-", "_", Breed), #remove dashes
         Breed = gsub("^\\s+", "",Breed), #remove leading white space,
         Breed = gsub("\\s*\\([^\\)]+\\)","", Breed), #remove names in parenthesis
         Breed = gsub("chinese_crested_dog", "chinese_crested",Breed),
         Breed = gsub("chinese_shar_pei", "chinese_sharpei",Breed),
         Breed = gsub("brittany_spaniel", "brittany", Breed),
         Breed = gsub("french_bull_dog", "bulldog_french", Breed),
         Breed = gsub("bull_mastiff", "bullmastiff", Breed),
         Breed = gsub("curly_coated_retriever", "curlycoated_retriever", Breed),
         Breed = gsub("greater_swiss_mountain", "greater_swiss_mountain_dog", Breed),
         Breed = gsub("japanese_chin", "japanese_chin_dog", Breed),
         Breed = gsub("standard_poodle", "poodle", Breed),
         Breed = gsub("italian_spinone", "spinone_italiano", Breed),
         Breed = gsub("wirehaired_fox_terrier", "wire_fox_terrier", Breed),
         Breed = gsub("mexican_hairless_dog", "xoloitzcuintli", Breed),
         Breed = gsub("peruvian_hairless_dog", "inca_hairless", Breed))

#write new file
#write.table(FinalCausalVars, file = "~/Documents/DogProject_Jaz/LocalRscripts/OMIA/processedCausalVarsOMIA.txt", quote = F, row.names = F, col.names = T, sep = "\t")
