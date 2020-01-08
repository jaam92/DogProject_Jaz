#Load libraries
library(tidyverse)
library(GMMAT)

#Function to inverse normal transform the score
qn = function(exp_vector) {
  result = qnorm(rank(exp_vector)/(length(exp_vector)+1))
  return(result)
}

#Load kinship matrix, ROHs, and phenotype data filename
pcRelateMat = readRDS("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/pcRelateMatrix_Unrelateds.rds")

rohs = read.delim(file = "~/Documents/DogProject_Jaz/LocalRscripts/ROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt") %>%
  group_by(INDV) %>%
  summarise(totalROH = sum(as.numeric(AUTO_LEN))) 

fnames = paste0("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/splitPhenotypeFile/", list.files(path="~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/splitPhenotypeFile", pattern = "\\.txt$")) #paste the path in front of the filename

#Loop through all files and run association test
AssociationTestResults = data.frame() #make data frame to store results
for (i in seq_along(fnames)){
  phenotypes = read.delim(fnames[i]) #read file in
  names(phenotypes)[3] = "trait" #rename third column to trait
  trait = gsub(".*[/]([^.]+)[.txt].*", "\\1", fnames[i]) #keep track of trait and breed tested
  
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
           status = gsub("1", "0", pheno), 
           status = gsub("2", "1", status),
           status = as.numeric(status)) #recode case-control status to 0 and 1
  
  #Count number of cases and controls
  CaseControlCount = ROHLoad %>% 
    count(status) %>% 
    pivot_wider(names_from = status, values_from = n) %>% 
    rename("Case" = "1", "Control" = "0")
  
  #Run the glmm as a logisitic, compute coeffcients and  p-value  
  out.logistic = glmmkin(status~QNSCORE, data=ROHLoad, kins=kinshipMat*2, id = "dogID", family=binomial(link="logit"))
  LogisticCoef = out.logistic$coef[2]
  LogisticStat = out.logistic$coef[2]/sqrt(out.logistic$cov[2,2])
  LogisticPval = ifelse(out.logistic$coef[2]/sqrt(out.logistic$cov[2,2])<0, pnorm(out.logistic$coef[2]/sqrt(out.logistic$cov[2,2]),lower=TRUE)*2,pnorm(out.logistic$coef[2]/sqrt(out.logistic$cov[2,2]),lower=FALSE)*2)
  
  #Save Association Test output case-control count
  LogRegOutput = cbind.data.frame(trait, CaseControlCount, LogisticCoef,LogisticStat,LogisticPval)
  AssociationTestResults = rbind.data.frame(AssociationTestResults, LogRegOutput)
}





