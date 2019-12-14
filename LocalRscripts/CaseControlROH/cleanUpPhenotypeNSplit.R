#Libraries
library(tidyverse)

#Function to make a list of data frames for a given phenotype of interest
balancePhenoPerBreed = function(phenoColName){
  #browser()
  ListOfFinalDFs = list() #make an empty list to hold my list of data frames
  phenoCol = enquo(phenoColName) #grab the phenotype column
  
  #Make a data frame with breed num case controls and indicator col for sampling
  dfList = phenotypes %>% 
  select(dogID, breed, !!phenoCol) %>% 
  filter(breed != "mix") %>% #remove mixed breed dogs
  na.omit() %>%
  group_by(breed,!!phenoCol) %>%
  count() %>%
  ungroup() %>%
  group_by(breed) %>%
  group_map(~ {
    
  .x %>%
      pivot_wider(names_from  = PSVA, values_from = n )%>%
      dplyr::rename(control=`1`, cases=`2`) %>%
      mutate(downSamp=case_when(
        control > cases ~ as.integer(0), 
        control == cases ~ as.integer(0),
        cases > control ~ as.integer(control)
    ))

}, keep=TRUE) #fxn to identify whether there may be more cases than controls and if there are output the sample size of controls 
  
  #Make a final data frame with balanced case-control if there were more cases than. Otherwise output the original data frame
  for(i in 1:length(dfList)){
    
    if(dfList[[i]]$downSamp != 0){
      finalDF = phenotypes %>%
        select(dogID, breed, !!phenoCol) %>% 
        filter(breed == dfList[[i]]$breed) %>% #remove mixed breed dogs
        na.omit() %>%
        group_by(!!phenoCol) %>%
        sample_n(size = dfList[[i]]$downSamp) #downsample cases to match controls
    }else{
      finalDF = phenotypes %>%
        select(dogID, breed, !!phenoCol) %>% 
        filter(breed == dfList[[i]]$breed) %>% 
        na.omit()}
    ListOfFinalDFs[[i]] = finalDF
  }
  
  
  return(ListOfFinalDFs)
}


#Load file
  popmapDryad = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/breeds_dryad.txt")
  phenotypes = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/phenotypes.txt")  %>% 
    mutate(PSVA = ifelse(is.na(PSVA), PSVA_yorkshireTerriers, PSVA), 
           MCT = ifelse(is.na(MCT), MCT_labradorRetrievers, MCT), 
           lymphoma = ifelse(is.na(lymphoma), lymphoma_goldenRetrievers, lymphoma),
           breed = popmapDryad$breed[match(dogID, popmapDryad$dogID)])


#Make data frames for each phenotype of interest and output new phenotype files
#Elbow Dysplasia
ED_allBreeds = balancePhenoPerBreed("ED")

#Collagen Disorder
CLLD_allBreeds = balancePhenoPerBreed("CLLD")

#Epilepsy Irish Wolfhounds
IrishWolfhounds = balancePhenoPerBreed("epilepsy_irishWolfhounds")

#Lymphoma all Breeds
#sample equal number of cases and controls
Lymphoma_allBreeds = balancePhenoPerBreed("lymphoma")

#Crohn's in Boxer and Bulldog
GC_allBreeds = balancePhenoPerBreed("GC_boxers_bulldogs")

#Mast Cell Tumor all Breeds 
MCT_allBreeds = balancePhenoPerBreed("MCT") 

#PSVA in all breeds
#Add some more from the controls in Yorkshire Terrier data
PSVA_allBreeds = balancePhenoPerBreed("PSVA")

#Mitral Valve data
MitralValve = balancePhenoPerBreed("MVD") 
