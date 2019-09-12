#####Set Directory
setwd("~/Documents/DogProject_Jaz/LocalRscripts/OMIA")

#####Load Libraries
library(dplyr)

#####Load Files
omiaGenes = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/OMIA/omia-genes_v2.txt", fill = TRUE)
causalVars = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/OMIA/causalVars_OMIA.txt", fill = TRUE)

#####Modify Files
######Reformat OMIA files
causalVars$NCBI_gene_ID = omiaGenes$ncbi_gene_id[match(causalVars$Gene, omiaGenes$gene_symbol)] #add ncbi gene id

sepBreedCausalVars = causalVars %>% mutate(Breed.s.= strsplit(as.character(Breed.s.), ",")) %>% unnest(Breed.s.) #split the breeds that are comma delimited into separate rows with same info

######Rename breeds
names(sepBreedCausalVars)[19] = "Breed"
sepBreedCausalVars$Breed = tolower(sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub(" ", "_", sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("-", "_", sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("^\\_","",sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("chinese_crested_dog", "chinese_crested",sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("chinese_shar_pei", "chinese_sharpei",sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("brittany_spaniel", "brittany", sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("french_bull_dog", "bulldog_french", sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("bull_mastiff", "bullmastiff", sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("curly_coated_retriever", "curlycoated_retriever", sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("greater_swiss_mountain", "greater_swiss_mountain_dog", sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("japanese_chin", "japanese_chin_dog", sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("standard_poodle", "poodle", sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("italian_spinone", "spinone_italiano", sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("wirehaired_fox_terrier", "wire_fox_terrier", sepBreedCausalVars$Breed)

#Pull the info I want
FinalCausalVars = sepBreedCausalVars %>% select("OMIA.ID.s.", "NCBI_gene_ID", "Gene","Breed") %>% rename(OMIA_ID = OMIA.ID.s.)

#write new file
#write.table(FinalCausalVars, file = "processedCausalVarsOMIA.txt", quote = F, row.names = F, col.names = T, sep = "\t")
