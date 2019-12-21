#Libraries
library(tidyverse)

#Function to make a list of data frames for a given phenotype of interest
balancePhenoPerBreed = function(phenoColName){
  ListOfFinalDFs = list() #make an empty list to hold my list of data frames
  phenoCol = enquo(phenoColName) #grab the phenotype column
  
  #Make a data frame with breed num case controls and indicator col for sampling
  #This fxn will identify whether there may be more cases than controls for each breed (list is per breed)
  #If there are more cases than controls downsample and output an equal number than match the sample size of controls 
  #Only output data if there are at least 10 cases and controls
  dfList = phenotypes %>% 
  select(dogID, breed, !!phenoCol) %>%
  filter(breed != "mix") %>% #remove mixed breed dogs
  na.omit() %>%
  dplyr::rename(trait = !!phenoCol) %>% #rename phenotype col makes things easier when pivotting dataframe 
  group_by(breed,trait) %>%
  count() %>%
  ungroup() %>%
  group_by(breed) %>%
  mutate(caseNcontrol = n()) %>% #count up number of occurences of each breed in data frame should be 2 (cases and controls) 
  filter(caseNcontrol == "2") %>% #make sure there are cases and controls
  group_map(~ {
    
  .x %>%
      pivot_wider(names_from = trait, values_from = n )%>%
      dplyr::rename(control=`1`, case=`2`) %>%
      mutate(downSamp=case_when(
        control > case ~ as.integer(0), 
        control == case ~ as.integer(0),
        case > control ~ as.integer(control))) %>%
      filter(control >= 10 & case >= 10) #keep only if there are at least 10 cases and controls

}, keep=TRUE)
  
  #Remove the empty dataframes from the list of dataframes
  dfList = dfList[map(dfList, function(x) dim(x)[1]) > 0]
  
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
        na.omit()
    }
    
    #make list of data frames per trait
    ListOfFinalDFs[[i]] = finalDF[,1:3] #remove weird extra col
    
    #write the final data frame to file
    breed = unique(finalDF$breed)
    write.table(finalDF, file=paste0("~/DogProject_Jaz/LocalRscripts/CaseControlROH/splitPhenotypeFile/", phenoColName, "_", breed, ".txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
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
set.seed(303)
#Elbow Dysplasia
ED_allBreeds = balancePhenoPerBreed("ED")

#Collagen Disorder
CLLD_allBreeds = balancePhenoPerBreed("CLLD")

#Epilepsy Irish Wolfhounds
IrishWolfhounds_allBreeds = balancePhenoPerBreed("epilepsy_irishWolfhounds")

#Lymphoma all Breeds
#sample equal number of cases and controls
Lymphoma_allBreeds = balancePhenoPerBreed("lymphoma")

#Crohn's in Boxer and Bulldog
GC_allBreeds = balancePhenoPerBreed("GC_boxers_bulldogs")

#Mast Cell Tumor all Breeds 
MCT_allBreeds = balancePhenoPerBreed("MCT") 

#PSVA in all breeds
PSVA_allBreeds = balancePhenoPerBreed("PSVA")

#Mitral Valve data
MitralValve_allBreeds = balancePhenoPerBreed("MVD") 