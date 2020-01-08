#Load libraries
library(tidyverse)
library(GMMAT)

#Define fxns

#Function to inverse normal transform the score
qn = function(exp_vector) {
  result = qnorm(rank(exp_vector)/(length(exp_vector)+1))
  return(result)
}

#Load kinship matrix, ROHs, and phenotype data
pcRelateMat = readRDS("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/pcRelateMatrix_Unrelateds.rds")

rohs = read.delim(file = "~/Documents/DogProject_Jaz/LocalRscripts/ROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt") %>%
  group_by(INDV) %>%
  summarise(totalROH = sum(as.numeric(AUTO_LEN))) 

fnames = list.files(path="~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/splitPhenotypeFile", pattern = "\\.txt$")

#Generate data frame
##create columns with fileName, population, and compute pi 
df = rbindlist(sapply(fnames, read.delim, simplify = FALSE), use.names = TRUE, idcol = "FileName")

phenotypes = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/splitPhenotypeFile/CLLD_german_shepherd_dog.txt")

names(phenotypes)[3] = "trait"
#keep only those rows and columns that correspond to dogs of interest
kinshipMat = pcRelateMat %>%
  as.matrix()
kinshipMat = kinshipMat[, colnames(kinshipMat) %in% phenotypes$dogID]
kinshipMat = kinshipMat[ rownames(kinshipMat) %in% phenotypes$dogID,]

#Grab the ROH information
ROHLoad = colnames(kinshipMat) %>% 
  as.data.frame() %>% 
  rename("dogID" = ".") %>%
  mutate(ROHLoad = rohs$totalROH[match(dogID, rohs$INDV)],
         QNSCORE = qn(as.numeric(ROHLoad)), 
         pheno = phenotypes$trait[match(dogID, phenotypes$dogID)],
         status = ifelse(pheno == "2", as.numeric(1), as.numeric(0)))

out.logistic = glmmkin(status~QNSCORE, data=ROHLoad, kins=kinshipMat*2, id = "dogID", family=binomial(link="logit"))

LogisticCoef = out.logistic$coef[2]
LogisticStat = out.logistic$coef[2]/sqrt(out.logistic$cov[2,2])
LogisticPval = ifelse(out.logistic$coef[2]/sqrt(out.logistic$cov[2,2])<0, pnorm(out.logistic$coef[2]/sqrt(out.logistic$cov[2,2]),lower=TRUE)*2,pnorm(out.logistic$coef[2]/sqrt(out.logistic$cov[2,2]),lower=FALSE)*2)
LogRegOutput = cbind.data.frame(LogisticCoef,LogisticStat,LogisticPval,fname[i])



#Make data frames for each phenotype of interest and run association test
#Elbow Dysplasia
ED_allBreeds = phenoKinship("ED")
runAssociation(ED_allBreeds, "ED")

#Collagen Disorder
CLLD_allBreeds = phenoKinship("CLLD")
runAssociation(CLLD_allBreeds, "CLLD")

#Epilepsy Irish Wolfhounds
IrishWolfhounds = phenoKinship("epilepsy_irishWolfhounds")
runAssociation(IrishWolfhounds, "epilepsy_irishWolfhounds")

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
runAssociation(GC_allBreeds, "GC_boxers_bulldogs")

#Mast Cell Tumor in Labrador Retrievers 
#sample equal number of cases and controls
LabradorRetrievers = phenoKinship("MCT_labradorRetrievers") 
runAssociation(LabradorRetrievers, "MCT_labradorRetrievers")

#Mast Cell Tumor all Breeds 
MCT_allBreeds = phenoKinship("MCT") 
runAssociation(MCT_allBreeds, "MCT")

#PSVA in Yorkshire Terriers
YorkshireTerriers = phenoKinship("PSVA_yorkshireTerriers")
runAssociation(YorkshireTerriers, "PSVA_yorkshireTerriers")

#PSVA in all breeds
#Add some more from the controls in Yorkshire Terrier data
PSVA_allBreeds = phenoKinship("PSVA")
runAssociation(PSVA_allBreeds, "PSVA")

#Mitral Valve data
MitralValve = phenoKinship("MVD") 
runAssociation(MitralValve, "MVD")
