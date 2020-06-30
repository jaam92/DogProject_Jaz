####This script will take in omia gene names and omia causal variant files and preprocess them for us to pull data from

#Set Directory and Load Libraries
setwd("~/Documents/DogProject_Jaz/LocalRscripts/OMIA")
library(tidyverse)

#fxn to negate in
`%nin%` = Negate(`%in%`)

#Load Files 
#oldData = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/OMIA/causalVars_OMIA.txt", fill = TRUE) 
omiaGenes = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/OMIA/omia-genes_v2.txt", fill = TRUE, stringsAsFactors = F)

omiaData = read_csv("~/Documents/DogProject_Jaz/LocalRscripts/OMIA/variants.csv") %>% 
  mutate(`Breed(s)` = strsplit(as.character(`Breed(s)`), ",")) %>% 
  unnest(`Breed(s)`) #split the breeds that are comma delimited into separate rows with same info

#ID traits that do not affect fitness
nonFitnessTraits = c("Blue eyes", "Bob tail", "Brown", "Coat colour dilution", "Coat colour, agouti", "Coat colour, black-and-tan", "Coat colour, dominant black", "Coat colour, white spotting", "Coat colour, white spotting, KIT-related", "Fawn or sable", "Furnishings (moustache and eyebrows)", "Grizzle", "Improper coat", "Long hair", "Muted, undefined, diluted-brownish hue", "No Merle pattern - diluted-brownish hue", "No Merle pattern - solid coat", "Minimal Merle, areas deleted to white, tweed", "Classic Merle","Recessive black", "Red/yellow coat", "Screw tail", "White coat colour", "White or cream", "Harlequin", "Liver")

#Fitness related traits
FinalCausalVars = omiaData %>%
  filter(`Variant Phenotype` %nin% nonFitnessTraits) %>% #remove nonfitness related traits
  select(`OMIA ID(s)`, `Breed(s)`, Gene, `Variant Phenotype`) %>%
  rename(OMIA_ID = `OMIA ID(s)`, Breed = `Breed(s)`, Phenotype = `Variant Phenotype`) %>%
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

#nonFitness related 
FinalCausalVars_nonFitness = omiaData %>%
  filter(`Variant Phenotype` %in% nonFitnessTraits) %>% #keep only nonfitness related traits
  select(`OMIA ID(s)`, `Breed(s)`, Gene, `Variant Phenotype`) %>%
  rename(OMIA_ID = `OMIA ID(s)`, Breed = `Breed(s)`, Phenotype = `Variant Phenotype`) %>%
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
#write.table(FinalCausalVars_nonFitness, file = "~/Documents/DogProject_Jaz/LocalRscripts/OMIA/processedCausalVarsOMIA_nonFitnessRelated.txt", quote = F, row.names = F, col.names = T, sep = "\t" )
