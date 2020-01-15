#Libraries
library(tidyverse)

#Functions to make a list of data frames for a given phenotype of interest for each breed and over all breeds
#all Breeds fxn
balancePheno = function(phenoColName){
  df = data.frame() #make an empty list to hold my list of data frames
  phenoCol = enquo(phenoColName) #grab the phenotype column
  
  #Make a data frame with breed num case controls and indicator col for sampling
  #This fxn will identify whether there may be more cases than controls for each breed (list is per breed)
  #If there are more cases than controls downsample and output an equal number than match the sample size of controls 
  #Only output data if there are at least 10 cases and controls
  df = phenotypes %>% 
    select(dogID, breed, !!phenoCol) %>%
    filter(dogID %in% indivs$dogID & breed != "mix") %>% #remove mixed breed dogs
    na.omit() %>%
    dplyr::rename(trait = !!phenoCol) %>% #rename phenotype col makes things easier when pivotting dataframe 
    group_by(trait) %>%
    count() %>%
    pivot_wider(names_from = trait, values_from = n) %>%
    dplyr::rename(control=`1`, case=`2`) %>%
    mutate(downSamp=case_when(
      control > case ~ as.integer(0), 
      control == case ~ as.integer(0),
      case > control ~ as.integer(control))) %>%
    filter(control >= 10 & case >= 10)
  
  #Downsample if needed and write to output file
  if(df$downSamp != 0){
    finalDF = phenotypes %>%
      select(dogID, breed, !!phenoCol) %>% 
      dplyr::rename(trait = !!phenoCol) %>% #rename column to get group_by to work
      filter(dogID %in% indivs$dogID & breed != "mix") %>% #use only dogs that pass QC and aren't mixed breed
      na.omit() %>%
      group_by(trait) %>%
      sample_n(size = df$downSamp) #downsample cases to match controls
    names(finalDF)[3] = phenoColName #put original column name back
  }else{
    finalDF = phenotypes %>%
      select(dogID, breed, !!phenoCol) %>% 
      filter(dogID %in% indivs$dogID & breed != "mix") %>% #use only dogs that pass QC and aren't mixed breed
      na.omit()
  }
  
  #write the final data frame to file
  write.table(finalDF, file=paste0("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/splitPhenotypeFile/", phenoColName, ".txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
  
  return(finalDF)
  }
  

#per breed fxn
balancePhenoPerBreed = function(phenoColName){
  dfList = list() #make an empty list to hold my list of data frames
  phenoCol = enquo(phenoColName) #grab the phenotype column
  
  #Make a data frame with breed num case controls and indicator col for sampling
  #This fxn will identify whether there may be more cases than controls for each breed (list is per breed)
  #If there are more cases than controls downsample and output an equal number than match the sample size of controls 
  #Only output data if there are at least 10 cases and controls
  dfList = phenotypes %>% 
  select(dogID, breed, !!phenoCol) %>%
  filter(dogID %in% indivs$dogID & breed != "mix") %>% #remove mixed breed dogs
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
  
  #Write each component of list to separate file and downsample if needed
  for(i in 1:length(dfList)){
    phenoCol = enquo(phenoColName)
    if(dfList[[i]]$downSamp != 0){
      finalDF = phenotypes %>%
        select(dogID, breed, !!phenoCol) %>% 
        dplyr::rename(trait = !!phenoCol) %>% #rename column to get group_by to work
        filter(dogID %in% indivs$dogID & breed == dfList[[i]]$breed) %>% #remove mixed breed dogs
        na.omit() %>%
        group_by(trait) %>%
        sample_n(size = dfList[[i]]$downSamp) #downsample cases to match controls
      names(finalDF)[3] = phenoColName #put original column name back
    }else{
      finalDF = phenotypes %>%
        select(dogID, breed, !!phenoCol) %>% 
        filter(dogID %in% indivs$dogID & breed == dfList[[i]]$breed) %>% 
        na.omit()
    }
    
    #write the final data frame to file
    breed = unique(finalDF$breed)
    write.table(finalDF, file=paste0("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/splitPhenotypeFile/", phenoColName, "_", breed, ".txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
    
  }

  return(dfList)
}

#Load files
popmapDryad = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/breeds_dryad.txt")
phenotypes = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/phenotypes.txt")  %>% 
    mutate(PSVA = ifelse(is.na(PSVA), PSVA_yorkshireTerriers, PSVA), 
           MCT = ifelse(is.na(MCT), MCT_labradorRetrievers, MCT), 
           lymphoma = ifelse(is.na(lymphoma), lymphoma_goldenRetrievers, lymphoma),
           breed = popmapDryad$breed[match(dogID, popmapDryad$dogID)])
indivs = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/Individuals_allBreeds_mergedFitakCornell.txt", col.names = c("dogID", "breed"))

#Make data frames for each phenotype of interest and output new phenotype files
set.seed(303)
#Elbow Dysplasia
ED_allBreeds = balancePheno("ED")
ED_perBreed = balancePhenoPerBreed("ED")

#Collagen Disorder
CLLD_allBreeds = balancePheno("CLLD")
CLLD_perBreed = balancePhenoPerBreed("CLLD")

#Epilepsy Irish Wolfhounds
IrishWolfhounds_allBreeds = balancePheno("epilepsy_irishWolfhounds")

#Lymphoma all Breeds
Lymphoma_allBreeds = balancePheno("lymphoma")
Lymphoma_perBreed = balancePhenoPerBreed("lymphoma")

#Crohn's in Boxer and Bulldog
GC_allBreeds = balancePheno("GC_boxers_bulldogs")
GC_perBreed = balancePhenoPerBreed("GC_boxers_bulldogs")

#Mast Cell Tumor all Breeds 
MCT_allBreeds = balancePheno("MCT") 
MCT_perBreed = balancePhenoPerBreed("MCT") 

#PSVA in all breeds
PSVA_allBreeds = balancePheno("PSVA")
PSVA_perBreed = balancePhenoPerBreed("PSVA")

#Mitral Valve data
MitralValve_allBreeds = balancePheno("MVD") 
MitralValve_perBreed = balancePhenoPerBreed("MVD") 
