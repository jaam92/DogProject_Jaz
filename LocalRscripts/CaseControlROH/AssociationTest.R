#Load libraries
library(tidyverse)

#Define fxns
#Function to pull phenotype data
phenoKinship = function(phenoColName){
  phenoCol = enquo(phenoColName)
  dogsOfInterest = phenotypes %>% 
    select(dogID, !!phenoCol) %>% 
    na.omit() 
  #keep only those rows and columns that correspond to dogs of interest
  kinshipMat = pcRelateMat %>%
    as.matrix()
  kinshipMat = kinshipMat[, colnames(kinshipMat) %in% dogsOfInterest$dogID]
  kinshipMat = kinshipMat[ rownames(kinshipMat) %in% dogsOfInterest$dogID,]
  return(kinshipMat)
}
#Run ROH score test
#function to inverse normal transform the score
qn <- function(exp_vector) {
  result = qnorm(rank(exp_vector)/(length(exp_vector)+1))
  return(result)
}
#function to run the association
runAssociation = function(phenoKinshipMatrix, phenotypeColName){
  ROHLoad = colnames(phenoKinshipMatrix) %>% 
    as.data.frame() %>% 
    rename("dogID" = ".") %>%
    mutate(ROHLoad = rohs$totalROH[match(dogID, rohs$INDV)],
           QNSCORE = qn(as.numeric(ROHLoad)), 
           pheno = phenotypes[,phenotypeColName][match(dogID, phenotypes$dogID)],
           status = ifelse(pheno == "2", as.numeric(1), as.numeric(0)))
  #test whether there are more controls than cases, if not balance the data
  if(table(ROHLoad$pheno)[1] > table(ROHLoad$pheno)[2]){
    finalROHLoad = ROHLoad
  }else{
    ###start here tomorrow
    finalROHLoad = 
  }
  out.logistic<-glmmkin(status~QNSCORE,data=finalROHLoad,kins=phenoKinshipMatrix*2,id = "dogID",family=binomial(link="logit"))
  
  LogisticCoef = out.logistic$coef[2]
  LogisticStat = out.logistic$coef[2]/sqrt(out.logistic$cov[2,2])
  LogisticPval = ifelse(out.logistic$coef[2]/sqrt(out.logistic$cov[2,2])<0, pnorm(out.logistic$coef[2]/sqrt(out.logistic$cov[2,2]),lower=TRUE)*2,pnorm(out.logistic$coef[2]/sqrt(out.logistic$cov[2,2]),lower=FALSE)*2)
  LogRegOutput = cbind.data.frame(LogisticCoef,LogisticStat,LogisticPval)
  return(LogRegOutput)
}

#Load kinship matrix, ROHs, and phenotype data
pcRelateMat <- readRDS("~/DogProject_Jaz/LocalRscripts/CaseControlROH/pcRelateMatrix_Unrelateds.rds")
rohs = read.delim(file = "~/DogProject_Jaz/LocalRscripts/ROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt") %>%
  group_by(INDV) %>%
  summarise(totalROH = sum(as.numeric(AUTO_LEN)))
phenotypes = read.delim("~/DogProject_Jaz/LocalRscripts/BreedCladeInfo/phenotypes.txt")  %>% 
  mutate(PSVA = ifelse(is.na(PSVA), PSVA_yorkshireTerriers, PSVA), 
         MCT = ifelse(is.na(MCT), MCT_labradorRetrievers, MCT), 
         lymphoma = ifelse(is.na(lymphoma), lymphoma_goldenRetrievers, lymphoma))

#Make data frames for each phenotype of interest and run association test
#Hip dysplasia
CHD_allBreeds = phenoKinship("CHD")

#Elbow Dysplasia
ED_allBreeds = phenoKinship("ED")

#Collagen Disorder
CLLD_allBreeds = phenoKinship("CLLD")

#Epilepsy Irish Wolfhounds
IrishWolfhounds = phenoKinship("epilepsy_irishWolfhounds")

#Lymphoma Golden Retrievers 
GoldenRetrievers = phenoKinship("lymphoma_goldenRetrievers")
runAssociation(GoldenRetrievers, "lymphoma_goldenRetrievers")

#Lymphoma all Breeds
#sample equal number of cases and controls
Lymphoma_allBreeds = phenoKinship("lymphoma")
runAssociation(Lymphoma_allBreeds, "lymphoma")

#Crohn's in Boxer
Boxers = phenoKinship("GC_boxers_bulldogs")

#Crohn's in Boxer and Bulldog
GC_allBreeds = phenoKinship("GC_boxers_bulldogs")

#Mast Cell Tumor in Labrador Retrievers 
#sample equal number of cases and controls
LabradorRetrievers = phenoKinship("MCT_labradorRetrievers") 

#Mast Cell Tumor all Breeds 
MCT_allBreeds = phenoKinship("MCT") 

#PSVA in Yorkshire Terriers
YorkshireTerriers = phenoKinship("PSVA_yorkshireTerriers")

#PSVA in all breeds
#Add some more from the controls in Yorkshire Terrier data
PSVA_allBreeds = phenoKinship("PSVA")

#Mitral Valve data
MitralValve = phenoKinship("MVD") 
