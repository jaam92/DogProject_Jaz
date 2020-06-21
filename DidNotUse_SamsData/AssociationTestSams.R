#Load libraries
library(data.table)
library(tidyverse)
library(reshape2)
library(GMMAT)

#####Reformat case control data to match hayward and pull only diseases with at least 10 cases also select only 300 controls
#Load Files
#cases = read.delim("~/Documents/DogProject_Jaz/DidNotUse_SamsData/at_risk_id_key.tsv", stringsAsFactors = F) %>%
#  separate_rows(disorder_name, sep = ",") %>%
#  mutate(disorder_name = gsub("^\\s+|\\s+$", "", disorder_name),
#         disorder_name = gsub("[()]", "", disorder_name),
#         disorder_name = gsub(";", "_aka_", disorder_name)) 

#set.seed(505)
#controls = read.delim("~/Documents/DogProject_Jaz/DidNotUse_SamsData/breed_id_key.tsv", stringsAsFactors = F) %>%
#  sample_n(300) %>% #sample 300 controls
#  select(IID)

#count number of cases 
#countCases = cases %>%
#  group_by(disorder_name) %>%
#  count() %>%
#  ungroup() %>%
#  filter(n >= 10)

#Reformat data to match hayward data
#reformatTraits = cases %>%
#  filter(disorder_name %in% countCases$disorder_name) #keep only traits with at least 10 cases
#  select(IID, disorder_name) %>%
#  mutate(status = "2",
#         disorder_name = gsub(" ", "_", disorder_name)) %>%
#  pivot_wider(names_from = disorder_name, values_from = status)

#Add columns with diseases and negative status for controls
#controls[setdiff(names(reformatTraits), names(controls))] <- "1"

#This fxn will balance cases and controls if need and write out diseases to seperate files
balancePhenoSams = function(phenoColName){
  df = data.frame() #make an empty list to hold my list of data frames
  phenoCol = enquo(phenoColName) #grab the phenotype column
  
  #Make a data frame with breed num case controls and indicator col for sampling
  #This fxn will identify whether there may be more cases than controls for each breed (list is per breed)
  #If there are more cases than controls downsample and output an equal number than match the sample size of controls 
  #Only output data if there are at least 10 cases and controls
  df = phenotypesSamsData %>% 
    select(dogID, !!phenoCol) %>%
    na.omit(!!phenoCol) %>%
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
    finalDF = phenotypesSamsData %>%
      select(dogID, !!phenoCol) %>% 
      dplyr::rename(trait = !!phenoCol) %>% #rename column to get group_by to work
      na.omit(!!phenoCol) %>%
      group_by(trait) %>%
      sample_n(size = df$downSamp) #downsample cases to match controls
    names(finalDF)[3] = phenoColName #put original column name back
  }else{
    finalDF = phenotypesSamsData %>%
      select(dogID, !!phenoCol) %>% 
      na.omit()
  }
  
  #write the final data frame to file
  #write.table(finalDF, file=paste0("~/Documents/DogProject_Jaz/DidNotUse_SamsData/splitPheno/", phenoColName, ".txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
  
  return(finalDF)
}

#Make data frames for each phenotype of interest and output new phenotype files
#set.seed(301)
#Canine_Multifocal_Retinopathy = balancePhenoSams("Canine_Multifocal_Retinopathy")
#Collie_Eye_Anomaly_aka_Choroidal_Hypoplasia = balancePhenoSams("Collie_Eye_Anomaly_aka_Choroidal_Hypoplasia")
#Degenerative_Myelopathy = balancePhenoSams("Degenerative_Myelopathy")
#Exercise_Induced_Collapse  = balancePhenoSams("Exercise-Induced_Collapse")    
#Factor_VII_Deficiency = balancePhenoSams("Factor_VII_Deficiency") 
#Hyperuricosuria_and_Hyperuricemia_or_Urolithiasis = balancePhenoSams("Hyperuricosuria_and_Hyperuricemia_or_Urolithiasis")  
#Ichthyosis = balancePhenoSams("Ichthyosis")
#Progressive_Retinal_Atrophy_PRA = balancePhenoSams("Progressive_Retinal_Atrophy_PRA")                     
#Von_Willebrand_Disease_TypeI = balancePhenoSams("Von_Willebrand_Disease_TypeI")

#####Making a kinship matrix based on ROH
#pcRelateMat = readRDS("~/Documents/DogProject_Jaz/DidNotUse_SamsData/pcRelateMatrix_allIndivs_SamsData.rds") %>% 
#  as.matrix()
#grmROHs = read_delim("~/Documents/DogProject_Jaz/DidNotUse_SamsData/ROHGRM_Sams.ped", col_names = F, delim = "\t") %>%
#  mutate(X1 = gsub(".bed","", X1, perl = T),
#         X2 = gsub(".bed","", X2, perl = T))
#allCompsGRM = grmROHs %>%
#  select(X2, X1, X3) %>% #flip comparisons
#  rename("X1" = "X2", "X2"= "X1")
#fullROHGRM = rbind.data.frame(allCompsGRM, grmROHs) %>%
#  rename("V1" = "X1", "V2"= "X2", "V3"="X3")#add original comps
#rm(grmROHs,allCompsGRM) #delete intermidiate df

#Compute bounded ROH burden for each individual
#ROHperIndiv = rohs %>%
#  group_by(INDV) %>%
#  summarise(ROHBurden = as.numeric(sum(AUTO_LEN))) %>% 
#  ungroup() %>%
#  mutate(ROHBurdenNorm = as.numeric((ROHBurden-min(ROHBurden))/(max(ROHBurden)-min(ROHBurden)))) %>%
#  select(INDV, ROHBurdenNorm) 

#longData = melt(pcRelateMat) %>%
#  mutate(V1 = as.character(Var1),
#         V2 = as.character(Var2),
#         kinshipNorm = as.numeric((value-min(value))/(max(value)-min(value))), 
#         rohNorm = 1-abs(ROHperIndiv$ROHBurdenNorm[match(V1, ROHperIndiv$INDV)] - ROHperIndiv$ROHBurdenNorm[match(V2, ROHperIndiv$INDV)])) %>%
#  select(V1, V2, kinshipNorm, rohNorm)

#Add kinship based on shared ROH
#left join the data with count of total bp of shared ROH
#replace NA with approximate length of genome so comparing self to self will be 1
#bound sharing of genome between 0 and 1
#last add breed info
#mergedDF = longData %>%
#  select(-c(kinshipNorm)) %>%
#  left_join(fullROHGRM) %>%
#  mutate(V3 = ifelse(is.na(V3), as.numeric(3e9), V3), 
#         rohGRMNorm = (V3-min(V3))/(max(V3)-min(V3))) %>%
#  select(-c(V3)) 

#convert to distance matrices and mantel test
#finalgrmROHMatrix = mergedDF %>%
#  select(V1, V2, rohGRMNorm) %>%
#  acast(V1~V2, value.var="rohGRMNorm")

#saveRDS(finalgrmROHMatrix, "~/Documents/DogProject_Jaz/DidNotUse_SamsData/rohGRM_allIndivs_SamsData.rds")




#####Run Association Test
#Function to inverse normal transform the score
qn = function(exp_vector) {
  result = qnorm(rank(exp_vector)/(length(exp_vector)+1))
  return(result)
}

#Load kinship matrix, ROHs, and phenotype data filename
pcRelateMat = readRDS("~/Documents/DogProject_Jaz/DidNotUse_SamsData/rohGRM_allIndivs_SamsData.rds")

rohs = read.delim(file = "~/Documents/DogProject_Jaz/DidNotUse_SamsData/TrueROH_propCoveredwithin1SDMean_SamsData_allChroms_vcfTools_rmROHlessThan50snps.txt", stringsAsFactors = F) %>%
  group_by(INDV) %>%
  summarise(totalROH = sum(as.numeric(AUTO_LEN))) 

#paste the path in front of the filename
fnames = paste0("~/Documents/DogProject_Jaz/DidNotUse_SamsData/splitPheno/", list.files(path="~/Documents/DogProject_Jaz/DidNotUse_SamsData/splitPheno/", pattern = "\\.txt$")) 

#Loop through all files and run association test and plot differences in ROH burden
AssociationTestResults = data.frame() #make data frame to store results
pdf(file ="~/Documents/DogProject_Jaz/DidNotUse_SamsData/ROHBurden_CaseControl.pdf", height = 8, width = 10) #open pdf to save differences in ROH burden

for (i in seq_along(fnames)){
  phenotypes = read.delim(fnames[i]) #read file in
  names(phenotypes)[2] = "trait" #rename third column to trait
  trait = gsub("_", " ", gsub(".*[/]([^.]+)[.txt].*", "\\1", fnames[i])) #keep track of trait and breed tested
  
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
  
  #Plot ROH burdens
  plotROHBurden = ggplot(ROHLoad, aes(x=pheno, y=ROHLoad/10^6, group=pheno, colour=pheno)) +
    geom_boxplot() +
    geom_point() +
    theme_bw() +
    labs(x ="Case-Control Status", y="ROH Burden (Mb)", title = paste(trait)) +
    scale_x_continuous(breaks = c(1,2), labels = c("1"="Controls", "2"="Cases")) +
    scale_color_gradient(low = "blue", high = "red") +
    theme(axis.text.x = element_text(size = 18), 
          axis.text.y = element_text(size = 18), 
          plot.title= element_text(size=24, face = "bold", hjust=0.5), 
          axis.title= element_text(size=20),
          strip.text = element_text(size = 14),
          legend.position = "none") 
  
  #put in pdf
  #print(plotROHBurden)
}
dev.off() #close pdf

#Plot the output
ggplot(AssociationTestResults, aes(x=gsub("_", " ", trait), y = LogBeta, colour= significant)) +
  geom_hline(yintercept = 0) + 
  geom_errorbar(aes(ymin=lowerBound, ymax=upperBound), colour="gray40", width=.2) + 
  geom_point() + 
  coord_flip() +  
  scale_colour_manual(values = c("yes"= "red", "no"="black")) + 
  labs(x = "Trait", y="Effect Size") +
  theme_bw() + 
  theme(axis.text.x = element_text( hjust= 0.5, vjust=1, size=20), 
        axis.text.y = element_text(size =20), 
        plot.title=element_text(size =24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=24), 
        legend.position = "none") 





