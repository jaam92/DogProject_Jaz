#Load Libraries
library(tidyverse)
library(data.table)
library(ggplot2)

#Load files
popmapDryad = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/breeds_dryad.txt")
popmapMerge = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/BreedAndCladeInfo_mergedFitakCornell.txt")
phenotypes = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/phenotypes.txt")
unrelateds = read.table("~/Documents/DogProject_Jaz/LocalRscripts/PCA_Unrelateds/Hayward2016_ogFiles/UnrelatedIndividuals_allBreeds.txt")
ROH = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt")
IBD = fread("~/Documents/DogProject_Jaz/LocalRscripts/IBDSegs/IBDSeq/MergedFitakCornell_allChroms_Haplotypes_IBDSeq.ibd")

#update files
phenotypes_unrelateds = phenotypes %>% mutate(PSVA = ifelse(is.na(PSVA), PSVA_yorkshireTerriers, PSVA), MCT = ifelse(is.na(MCT), MCT_labradorRetrievers, MCT), lymphoma = ifelse(is.na(lymphoma), lymphoma_goldenRetrievers, lymphoma), Unrelated = ifelse(dogID %in% unrelateds$V1, 1, 0), breed = popmapDryad$breed[match(dogID, popmapDryad$dogID)]) %>% filter(Unrelated == 1)

popmapMerge$breed = gsub("large_munsterlander","munsterlander_large", popmapMerge$breed)

#Aggregate ROH Segments for unrelateds
aggregateROH = ROH %>% select(AUTO_LEN, INDV) %>% group_by(INDV) %>% summarise(totalLenROH = sum(AUTO_LEN)) %>% mutate(totalLenROHMb = totalLenROH/10^6) %>% filter(INDV %in% phenotypes_unrelateds$dogID) %>% as.data.frame()
rm(ROH)

#Aggregate IBD Segments
aggregateIBD = IBD %>% filter(V1 %in% phenotypes_unrelateds$dogID) %>% mutate(Breed1 = popmapMerge$breed[match(V1, popmapMerge$dogID)], Breed2 = popmapMerge$breed[match(V3, popmapMerge$dogID)]) %>% filter(Breed1 == Breed2) %>% mutate(segLen = as.numeric(V7) - as.numeric(V6)) %>% group_by(V1) %>% summarise(totalLenIBD = sum(segLen)) %>% mutate(totalLenIBDMb = totalLenIBD/10^6) %>% plyr::rename(c("V1" = "INDV")) %>% as.data.frame()


#Phenotypes with more into
phenotypes_unrelateds$totalLenIBDMb = aggregateIBD$totalLenIBDMb[match(phenotypes_unrelateds$dogID, aggregateIBD$INDV)]
phenotypes_unrelateds$totalLenROHMb = aggregateROH$totalLenROHMb[match(phenotypes_unrelateds$dogID, aggregateROH$INDV)]

#Hip dysplasia
CHD_allBreeds = phenotypes_unrelateds %>% select("dogID", "breed", "CHD", "totalLenIBDMb", "totalLenROHMb") %>% na.omit()

ggplot(CHD_allBreeds, aes(x=totalLenROHMb, y=CHD)) + geom_point() + geom_smooth(method = "lm") + theme_bw() #plot looks like a cloud I doubt there is any correlation

#Elbow Dysplasia
ED_allBreeds = phenotypes_unrelateds %>% select("dogID", "breed", "ED", "totalLenIBDMb", "totalLenROHMb") %>% plyr::rename(c( "ED" = "status")) %>% na.omit()

ggplot(ED_allBreeds, aes(x=status, y= totalLenROHMb, group=status)) + geom_boxplot() + scale_x_discrete( limits= c("1", "2")) + theme_bw() 

#Collagen Disorder
CLLD_allBreeds = phenotypes_unrelateds %>% select("dogID", "breed", "CLLD", "totalLenIBDMb", "totalLenROHMb") %>% plyr::rename(c( "CLLD" = "status")) %>% na.omit()

ggplot(CLLD_allBreeds, aes(x=status, y= totalLenROHMb, group=status)) + geom_boxplot() + scale_x_discrete( limits= c("1", "2")) + theme_bw() 


#Epilepsy Irish Wolfhounds
IrishWolfhounds = phenotypes_unrelateds[grep("wolf", phenotypes_unrelateds$breed),] %>% select("dogID", "breed", "epilepsy_irishWolfhounds", "totalLenIBDMb", "totalLenROHMb") %>% plyr::rename(c( "epilepsy_irishWolfhounds" = "status")) %>% na.omit()

ggplot(IrishWolfhounds, aes(x=status, y= totalLenROHMb, group=status)) + geom_boxplot() + scale_x_discrete( limits= c("1", "2")) + theme_bw() 

#Lymphoma Golden Retrievers 
GoldenRetrievers = phenotypes_unrelateds[grep("golden", phenotypes_unrelateds$breed),] %>% select("dogID", "breed", "lymphoma_goldenRetrievers", "totalLenIBDMb", "totalLenROHMb") %>% plyr::rename(c( "lymphoma_goldenRetrievers" = "status")) %>% na.omit()

ggplot(GoldenRetrievers, aes(x=status, y= totalLenROHMb, group=status)) + geom_boxplot() + scale_x_discrete( limits= c("1", "2")) + theme_bw() 

#Lymphoma all Breeds
#sample equal number of cases and controls
Lymphoma_allBreeds = phenotypes_unrelateds %>% select("dogID", "breed", "lymphoma", "totalLenIBDMb", "totalLenROHMb") %>% plyr::rename(c( "lymphoma" = "status")) %>% na.omit() %>% group_by(status) %>%  sample_n(74)
  
ggplot(Lymphoma_allBreeds, aes(x=status, y= totalLenROHMb, group=status)) + geom_boxplot() + scale_x_discrete( limits= c("1", "2")) + theme_bw() 

#Crohn's in Boxer
Boxers = phenotypes_unrelateds[grep("boxer", phenotypes_unrelateds$breed),] %>% select("dogID", "breed", "GC_boxers_bulldogs", "totalLenIBDMb", "totalLenROHMb") %>% plyr::rename(c( "GC_boxers_bulldogs" = "status")) %>% na.omit()


ggplot(Boxers, aes(x=status, y= totalLenROHMb, group=status)) + geom_boxplot() + scale_x_discrete( limits= c("1", "2")) + theme_bw()

#Crohn's in Boxer and Bulldog
GC_allBreeds = phenotypes_unrelateds %>% select("dogID", "breed", "GC_boxers_bulldogs", "totalLenIBDMb", "totalLenROHMb") %>% plyr::rename(c( "GC_boxers_bulldogs" = "status")) %>% na.omit()

ggplot(GC_allBreeds, aes(x=status, y= totalLenROHMb, group=status)) + geom_boxplot() + scale_x_discrete( limits= c("1", "2")) + theme_bw()

#Mast Cell Tumor in Labrador Retrievers 
#sample equal number of cases and controls
LabradorRetrievers = phenotypes_unrelateds[grep("labrador", phenotypes_unrelateds$breed),] %>% select("dogID", "breed", "MCT_labradorRetrievers", "totalLenIBDMb", "totalLenROHMb") %>%  plyr::rename(c( "MCT_labradorRetrievers" = "status")) %>% na.omit() %>% group_by(status) %>%  sample_n(44)

ggplot(LabradorRetrievers, aes(x=status, y= totalLenROHMb, group=status)) + geom_boxplot() + scale_x_discrete( limits= c("1", "2")) + theme_bw() 

#Mast Cell Tumor all Breeds 
#sample equal number of cases and controls
MCT_allBreeds = phenotypes_unrelateds %>% select("dogID", "breed", "MCT", "totalLenIBDMb", "totalLenROHMb") %>%  plyr::rename(c( "MCT" = "status")) %>% na.omit() %>% group_by(status) %>%  sample_n(55)

ggplot(MCT_allBreeds, aes(x=status, y= totalLenROHMb, group=status)) + geom_boxplot() + scale_x_discrete( limits= c("1", "2")) + theme_bw() 

#PSVA in Yorkshire Terriers
YorkshireTerriers = phenotypes_unrelateds[grep("yorkshire", phenotypes_unrelateds$breed),] %>% select("dogID", "breed", "PSVA_yorkshireTerriers", "totalLenIBDMb", "totalLenROHMb") %>% plyr::rename(c( "PSVA_yorkshireTerriers" = "status")) %>% na.omit()

ggplot(YorkshireTerriers, aes(x=status, y= totalLenROHMb, group=status)) + geom_boxplot() + scale_x_discrete( limits= c("1", "2")) + theme_bw()

#PSVA in all breeds
#Add some more from the controls in Yorkshire Terrier data
PSVA_allBreeds = phenotypes_unrelateds %>% select("dogID", "breed", "PSVA","PSVA_yorkshireTerriers","totalLenIBDMb", "totalLenROHMb") %>% plyr::rename(c( "PSVA" = "status")) %>% mutate(status = ifelse(is.na(status), PSVA_yorkshireTerriers, status)) %>% select(-one_of("PSVA_yorkshireTerriers")) %>% na.omit() %>% group_by(status) %>% sample_n(65)

ggplot(PSVA_allBreeds, aes(x=status, y= totalLenROHMb, group=status)) + geom_boxplot() + scale_x_discrete( limits= c("1", "2")) + theme_bw()

#Plot Mitral Valve data
MitralValve = phenotypes_unrelateds %>% select("dogID", "breed", "MVD", "totalLenIBDMb", "totalLenROHMb") %>% plyr::rename(c("MVD" = "status")) %>% na.omit() %>% group_by(status) %>% sample_n(73)

ggplot(MitralValve, aes(x=status, y= totalLenROHMb, group=status)) + geom_boxplot() + scale_x_discrete( limits= c("1", "2")) + theme_bw() 

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

dataframePvaluesPerms = t(cbind.data.frame(IrishWolfhound_Perms[[3]], YorkshireTerriers_Perms[[3]], PSVA_allBreeds_Perms[[3]], LabradorRetrievers_Perms[[3]], MCT_allBreeds_Perms[[3]], GC_allBreeds_Perms[[3]], Boxers_Perms[[3]], GoldenRetrievers_Perms[[3]], Lymphoma_allBreeds_Perms[[3]], MitralValve_Perms[[3]], multipleDiagnosis_ED_CLLD_Perms[[3]]))

dataframePvaluesPerms

#Plot significant stuff

ggplot(GoldenRetrievers_Perms[[1]], aes(GoldenRetrievers_Perms[[1]]$permutations)) + geom_histogram(binwidth = 30, breaks=seq(-150, 200, by =5),col="coral2", fill="white") + geom_vline(xintercept = GoldenRetrievers_Perms[[2]], col="purple") + theme_bw() + theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size  = 20), plot.title=element_text(size=26, face = "bold"), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) + labs(x= "Permutation Score", y= "Count") + ggtitle("ROH Burden in Golden Retriever Lymphoma Case vs Control")

ggplot(Lymphoma_allBreeds_Perms[[1]], aes(Lymphoma_allBreeds_Perms[[1]]$permutations)) + geom_histogram(binwidth = 30, breaks=seq(-150, 150, by =5),col="coral2", fill="white") + geom_vline(xintercept = Lymphoma_allBreeds_Perms[[2]], col="purple") + theme_bw() + theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size  = 20), plot.title=element_text(size=26, face = "bold"), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) + labs(x= "Permutation Score", y= "Count") + ggtitle("ROH Burden in Lymphoma Case vs Control")

ggplot(MitralValve_Perms[[1]], aes(MitralValve_Perms[[1]]$permutations)) + geom_histogram(binwidth = 30, breaks=seq(-250, 250, by =5) ,col="coral2", fill="white") + geom_vline(xintercept = MitralValve_Perms[[2]], col="purple") + theme_bw() + theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size  = 20), plot.title=element_text(size=26, face = "bold"), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) + labs(x= "Permutation Score", y= "Count") + ggtitle("ROH Burden in Mitral Valve Degeneration Case vs Control")
