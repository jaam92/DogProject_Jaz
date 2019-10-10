#Load Libraries
library(tidyverse)
library(data.table)

#Load files
popmapDryad = read.delim("~/DogProject_Jaz/LocalRscripts/BreedCladeInfo/breeds_dryad.txt")
popmapMerge = read.delim("~/DogProject_Jaz/LocalRscripts/BreedCladeInfo/BreedAndCladeInfo_mergedFitakCornell.txt") %>%
  mutate(breed = gsub("large_munsterlander","munsterlander_large", breed))
phenotypes = read.delim("~/DogProject_Jaz/LocalRscripts/BreedCladeInfo/phenotypes.txt")
unrelateds = read.table("~/DogProject_Jaz/LocalRscripts/PCA_Unrelateds/Hayward2016_ogFiles/UnrelatedIndividuals_allBreeds.txt")
ROH = read.delim("~/DogProject_Jaz/LocalRscripts/CaseControlROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt")
IBD = fread("~/DogProject_Jaz/LocalRscripts/IBDSegs/IBDSeq/MergedFitakCornell_allChroms_Haplotypes_IBDSeq.ibd")

#Define fxns
#Function to pull phenotype data
phenoData = function(phenoColName){
  phenoCol = enquo(phenoColName)
  phenotypes_unrelateds %>% 
    select("dogID", "breed", !!phenoCol, "totalLenIBDMb", "totalLenROHMb") %>% 
    rename("status" = phenoColName ) %>%
    na.omit()
}

#Function to plot phenodata boxplots
plotPhenoBoxPlot = function(dataFrame){
  ggplot(dataFrame, aes(x=status, y= totalLenROHMb, group=status)) + 
    geom_boxplot() + 
    scale_x_discrete( limits= c("1", "2")) + 
    labs(y = "Length of Genome\nin ROH(Mb)", x="Case-Control Status") + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 20), 
          axis.text.y = element_text(size = 10), 
          plot.title=element_text(size=26, face = "bold", hjust=0.5), 
          axis.title=element_text(size=20),
          legend.title=element_text(size=20), 
          legend.text=element_text(size=18)) 
}

#Function to run permutation test on case-control data
#case-control status should be supplied in quotes
PermutationTest = function(dataFrame, caseIndicator, controlIndicator, numberPerms){
  ROHBurden = dataFrame$totalLenROHMb
  group = as.character(dataFrame$status)
  testStat = function(w, g) mean(as.numeric(w[g == caseIndicator])) - mean(as.numeric(w[g == controlIndicator])) #test whether there is a difference in median amount of genome in ROH btwn cases and controls
  observedStat = testStat(ROHBurden, group)
  permutations = sapply(1 : numberPerms, function(i) testStat(ROHBurden, sample(group)))
  pvalue = mean(permutations > observedStat) #pvalue 
  pvalue[pvalue == 0] = 1/numberPerms+1 #generate pvalue when there are no perms greater than observed
  bins = seq(1:numberPerms)
  dfPerms = cbind.data.frame(bins,permutations)
  return(list(dfPerms, observedStat, pvalue))
}

#Fxn to plot permutation results
plotPerms = function(dataFrame, plotTitle, breakStart, breakEnd){
  ggplot(dataFrame[[1]], aes(dataFrame[[1]]$permutations)) + 
    geom_histogram(binwidth = 30, breaks=seq(breakStart, breakEnd, by =5),
                   col="coral2", fill="white") + 
    geom_vline(xintercept = dataFrame[[2]], col="purple") +
    labs(y = "Count", x="Permutation Score", title = paste(plotTitle)) + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 20), 
          axis.text.y = element_text(size = 10), 
          plot.title=element_text(size=26, face = "bold", hjust=0.5), 
          axis.title=element_text(size=20),
          legend.title=element_text(size=20), 
          legend.text=element_text(size=18))
}

#Grab phenotypes for unrelateds
phenotypes_unrelateds = phenotypes %>% 
  mutate(PSVA = ifelse(is.na(PSVA), PSVA_yorkshireTerriers, PSVA), 
         MCT = ifelse(is.na(MCT), MCT_labradorRetrievers, MCT), 
         lymphoma = ifelse(is.na(lymphoma), lymphoma_goldenRetrievers, lymphoma),
         Unrelated = ifelse(dogID %in% unrelateds$V1, 1, 0), 
         breed = popmapDryad$breed[match(dogID, popmapDryad$dogID)]) %>% 
  filter(Unrelated == 1) 

#Aggregate ROH Segments for unrelateds
aggregateROH = ROH %>% 
  filter(INDV %in% phenotypes_unrelateds$dogID) %>% 
  select(AUTO_LEN, INDV) %>% 
  group_by(INDV) %>% 
  summarise(totalLenROH = sum(AUTO_LEN)) %>% 
  mutate(totalLenROHMb = totalLenROH/10^6)  

rm(ROH)#delete data frame

#Generate an IBD data frame with scores per individual
#Only using IBD segments within breed and greater than 4Mb 
#chose 4Mb rather than 3Mb because that's the cutoff used to reliably call segments from array data with IBDSeq
#Remove IBD segments that are not shared within the same breed of dog/wolf
aggregateIBD = IBD %>% 
  filter(V1 %in% phenotypes_unrelateds$dogID) %>% 
  mutate(Breed1 = popmapMerge$breed[match(V1, popmapMerge$dogID)], 
         Breed2 = popmapMerge$breed[match(V3, popmapMerge$dogID)],
         segLen = as.numeric(V7) - as.numeric(V6)) %>% 
  filter(segLen >= 4e6 & Breed1 == Breed2) %>% 
  group_by(V1) %>% 
  summarise(totalLenIBD = sum(segLen)) %>% 
  mutate(totalLenIBDMb = totalLenIBD/10^6) %>% 
  rename("INDV" = V1) 

rm(IBD) #delete data frame

#Add ROH and IBD Scores
phenotypes_unrelateds$totalLenIBDMb = aggregateIBD$totalLenIBDMb[match(phenotypes_unrelateds$dogID, aggregateIBD$INDV)]
phenotypes_unrelateds$totalLenROHMb = aggregateROH$totalLenROHMb[match(phenotypes_unrelateds$dogID, aggregateROH$INDV)]

#Make data frames for each phenotype of interest
#Hip dysplasia
CHD_allBreeds = phenoData("CHD")

#Elbow Dysplasia
ED_allBreeds = phenoData("ED")
plotPhenoBoxPlot(ED_allBreeds)

#Collagen Disorder
CLLD_allBreeds = phenoData("CLLD")
plotPhenoBoxPlot(CLLD_allBreeds)

#Epilepsy Irish Wolfhounds
IrishWolfhounds = phenoData("epilepsy_irishWolfhounds")
plotPhenoBoxPlot(IrishWolfhounds)

#Lymphoma Golden Retrievers 
GoldenRetrievers = phenoData("lymphoma_goldenRetrievers")
plotPhenoBoxPlot(GoldenRetrievers)

#Lymphoma all Breeds
#sample equal number of cases and controls
set.seed(10)
Lymphoma_allBreeds = phenoData("lymphoma") %>% 
  group_by(status) %>%  
  sample_n(74)
plotPhenoBoxPlot(Lymphoma_allBreeds)

#Crohn's in Boxer
Boxers = phenoData("GC_boxers_bulldogs") %>%
  filter(breed == "boxer")
plotPhenoBoxPlot(Boxers)

#Crohn's in Boxer and Bulldog
GC_allBreeds = phenoData("GC_boxers_bulldogs")
plotPhenoBoxPlot(GC_allBreeds)

#Mast Cell Tumor in Labrador Retrievers 
#sample equal number of cases and controls
set.seed(11)
LabradorRetrievers = phenoData("MCT_labradorRetrievers") %>% 
  group_by(status) %>%  
  sample_n(44)
plotPhenoBoxPlot(LabradorRetrievers)

#Mast Cell Tumor all Breeds 
#sample equal number of cases and controls
set.seed(12)
MCT_allBreeds = phenoData("MCT") %>% 
  group_by(status) %>%  
  sample_n(55)
plotPhenoBoxPlot(MCT_allBreeds) 

#PSVA in Yorkshire Terriers
YorkshireTerriers = phenoData("PSVA_yorkshireTerriers")
plotPhenoBoxPlot(YorkshireTerriers)

#PSVA in all breeds
#Add some more from the controls in Yorkshire Terrier data
set.seed(13)
PSVA_allBreeds = phenotypes_unrelateds %>% 
  select("dogID", "breed", "PSVA","PSVA_yorkshireTerriers","totalLenIBDMb", "totalLenROHMb") %>% 
  rename("status"="PSVA") %>% 
  mutate(status = ifelse(is.na(status), PSVA_yorkshireTerriers, status)) %>%
  select(-one_of("PSVA_yorkshireTerriers")) %>% 
  na.omit() %>% 
  group_by(status) %>% 
  sample_n(65)

plotPhenoBoxPlot(PSVA_allBreeds)

#Plot Mitral Valve data
set.seed(14)
MitralValve = phenoData("MVD") %>% 
  group_by(status) %>% 
  sample_n(73)
plotPhenoBoxPlot(MitralValve) 

#Check whether individuals that are cases for more than one trait have more ROH
#only traits where we have multiple diagnoses are ED and CLLD
#N=154 (14 cases and 140 controls)
multipleDiagnosis_ED_CLLD_controls = phenotypes_unrelateds %>% 
  select("dogID", "breed", "ED", "CLLD", "totalLenROHMb") %>% 
  filter(ED == 1 & CLLD == 1) %>%
  mutate(status = "1") %>%
  na.omit() 

multipleDiagnosis_ED_CLLD_cases = phenotypes_unrelateds %>% 
  select("dogID", "breed", "ED", "CLLD", "totalLenROHMb") %>% 
  filter(ED == 2 & CLLD == 2) %>%
  mutate(status = "2") %>%
  na.omit()

multipleDiagnosis_ED_CLLD = rbind.data.frame(multipleDiagnosis_ED_CLLD_controls, multipleDiagnosis_ED_CLLD_cases)

#Run Permutation Test to test whether there is a significant difference in ROH burden btwn cases and controls
set.seed(800)
ED_Perms = PermutationTest(ED_allBreeds, "2", "1", 10000)
CLLD_Perms = PermutationTest(CLLD_allBreeds, "2", "1", 10000)
IrishWolfhound_Perms = PermutationTest(IrishWolfhounds, "2", "1", 10000)
YorkshireTerriers_Perms = PermutationTest(YorkshireTerriers, "2", "1", 10000)
PSVA_allBreeds_Perms = PermutationTest(PSVA_allBreeds, "2", "1", 10000)
LabradorRetrievers_Perms = PermutationTest(LabradorRetrievers, "2", "1", 10000)
MCT_allBreeds_Perms = PermutationTest(MCT_allBreeds, "2", "1", 10000 )
Boxers_Perms = PermutationTest(Boxers, "2", "1", 10000)
GC_allBreeds_Perms = PermutationTest(GC_allBreeds, "2", "1", 10000)
GoldenRetrievers_Perms = PermutationTest(GoldenRetrievers, "2", "1", 10000)
Lymphoma_allBreeds_Perms = PermutationTest(Lymphoma_allBreeds, "2", "1", 10000 )
MitralValve_Perms = PermutationTest(MitralValve, "2", "1", 10000)
multipleDiagnosis_ED_CLLD_Perms = PermutationTest(multipleDiagnosis_ED_CLLD, "2", "1", 10000)

#Make data frame of pvalues and print results
dataframePvaluesPerms = t(cbind.data.frame(IrishWolfhound_Perms[[3]], YorkshireTerriers_Perms[[3]], PSVA_allBreeds_Perms[[3]], LabradorRetrievers_Perms[[3]], MCT_allBreeds_Perms[[3]], GC_allBreeds_Perms[[3]], Boxers_Perms[[3]], GoldenRetrievers_Perms[[3]], Lymphoma_allBreeds_Perms[[3]], MitralValve_Perms[[3]], multipleDiagnosis_ED_CLLD_Perms[[3]]))

cat(sprintf("permutation test results"), sep = "\n")
dataframePvaluesPerms 

#Plot permutations
plotPerms(GoldenRetrievers_Perms, "ROH Burden in Golden Retriever Lymphoma Case vs Control", -150, 200)
plotPerms(Lymphoma_allBreeds_Perms, "ROH Burden in Lymphoma Case vs Control", -150, 150)
