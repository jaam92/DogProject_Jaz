####You will want to replace all the file names and paths here 

#Load libraries
library(tidyverse)
library(GMMAT)

#Function to inverse normal transform the ROH score
qn = function(exp_vector) {
  result = qnorm(rank(exp_vector)/(length(exp_vector)+1))
  return(result)
}

#Load kinship matrix, ROHs, and phenotype data filename
pcRelateMat = readRDS("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/rohGRM_allIndivs.rds")

#Load ROH filter to only those greater than or equal to 2MB
rohs = read.delim(file = "~/Documents/DogProject_Jaz/LocalRscripts/ROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFile_Cornell_allChroms_vcfToolsROH_rmROHlessThan50snps_HaywardDataOnly.txt") %>%
  filter(AUTO_LEN >= 2e+06) %>%
  group_by(INDV) %>%
  summarise(totalROH = sum(as.numeric(AUTO_LEN))) #aggregate roh for each individual to get total amount of the genome within a ROH per individual

fnames = paste0("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/splitPhenotypeFile/IncludeMixedBreeds/", list.files(path="~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/splitPhenotypeFile/IncludeMixedBreeds/", pattern = "\\.txt$")) #grab all the file names from a specific folder and paste the path in front of the filename

#Loop through all files and run association test
AssociationTestResults = data.frame() #make data frame to store results
for (i in seq_along(fnames)){
  phenotypes = read.delim(fnames[i]) #read file in
  names(phenotypes)[3] = "trait" #rename third column to trait
  trait = gsub(".*[/]([^.]+)[.txt].*", "\\1", fnames[i]) #keep track of trait and breed tested
  
  #keep only those rows and columns of kinship matrix that correspond to dogs of interest
  kinshipMat = pcRelateMat %>%
    as.matrix()
  kinshipMat = kinshipMat[, colnames(kinshipMat) %in% phenotypes$dogID]
  kinshipMat = kinshipMat[ rownames(kinshipMat) %in% phenotypes$dogID,]
  
  #Grab the ROH information
  ROHLoad = colnames(kinshipMat) %>% 
    as.data.frame() %>% 
    rename("dogID" = ".") %>%
    mutate(ROHLoad = rohs$totalROH[match(dogID, rohs$INDV)],
           QNSCORE = qn(as.numeric(ROHLoad)), #inverse quantile normalize ROH for each trait
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
  LogBeta = out.logistic$coef[2]
  OddsRatio = exp(LogBeta)
  LogBetaStandardized = out.logistic$coef[2]/sqrt(out.logistic$cov[2,2]) #standardized
  LogisticPval = ifelse(out.logistic$coef[2]/sqrt(out.logistic$cov[2,2])<0, pnorm(out.logistic$coef[2]/sqrt(out.logistic$cov[2,2]),lower=TRUE)*2,pnorm(out.logistic$coef[2]/sqrt(out.logistic$cov[2,2]),lower=FALSE)*2)
  confintUpper = out.logistic$coefficients + 1.96*sqrt(diag(out.logistic$cov)) #diagonal of cov matrix gives variance need SE so take the sqrt
  confintLower = out.logistic$coefficients - 1.96*sqrt(diag(out.logistic$cov))
  
  #Save Association Test output case-control count and the confidence interval of my fixed effect beta
  LogRegOutput = cbind.data.frame(trait, CaseControlCount, OddsRatio, LogBeta, LogBetaStandardized, confintUpper[2], confintLower[2], LogisticPval) %>%
    dplyr::rename(lowerBound = `confintLower[2]`, upperBound = `confintUpper[2]`)
  LogRegOutput$significant = ifelse(LogRegOutput$LogisticPval <= 0.05, "yes","no") #attach indicator 
  AssociationTestResults = rbind.data.frame(AssociationTestResults, LogRegOutput)
}

#write to file
#write.table(AssociationTestResults, file = "~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/AssociationTestResults.txt", col.names = T, row.names = F, quote = F, sep = "\t")

#regular expressions to clean up trait name and make first letter uppercase 
AssociationTestResults$trait = gsub("(?<=^|_)([a-z])", "\\U\\1", AssociationTestResults$trait, perl=TRUE)
AssociationTestResults$trait = gsub("_Dog", "_dog", AssociationTestResults$trait, perl = TRUE)

#Plot the results
ggplot(AssociationTestResults, aes(x=gsub("_", " ", trait), y = LogBeta, colour= significant)) +
  geom_hline(yintercept = 0) + 
  geom_errorbar(aes(ymin=lowerBound, ymax=upperBound), colour="gray40", width=.2) + 
  geom_point() + 
  coord_flip() +  
  scale_colour_manual(values = c("yes"= "red", "no"="black")) + 
  labs(x = "Trait", y="Effect Size (log-odds)") +
  theme_bw() + 
  theme(axis.text.x = element_text( hjust= 0.5, vjust=1, size=20), 
        axis.text.y = element_text(size =20), 
        plot.title=element_text(size =24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=24), 
        legend.position = "none") 
  




