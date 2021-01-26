#Load libraries
library(tidyverse)
library(data.table)
library(mgsub)

#Load files
setwd("~/Documents/DogProject_Jaz/LocalRscripts/ROH/ROHOverlapIBDperChrom/")
sharedROHIBD = read.delim("FinalSharedROHwithinIBDSegs_allChroms.bed", col.names = c("chrom","ROHShareStart", "ROHShareEnd", "INDV1", "INDV2","chrom","IBDStart", "IBDEnd", "INDV1", "INDV2","LenOverlap"), fill = NA) %>% 
  select(-c("INDV1.1", "INDV2.1","chrom.1")) %>%
  mutate(LengthIBDSeg = as.numeric(IBDEnd) - as.numeric(IBDStart),
         NotCoveredROH = as.numeric(LengthIBDSeg) - as.numeric(LenOverlap))
popmapDryad = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/breeds_dryad.txt")
phenotypes = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/phenotypes.txt")
unrelateds = read.table("~/Documents/DogProject_Jaz/LocalRscripts/PCA_Unrelateds/UnrelatedIndividuals_allBreeds_mergedFitakCornell.txt")

fnames = paste0("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/splitPhenotypeFile/IncludeMixedBreeds/", 
                list.files(path="~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/splitPhenotypeFile/IncludeMixedBreeds/", 
                           pattern = "[^_]") %>%
                  str_subset(., "_", negate = TRUE)) #remove breed specific files

####Function for computing normalized score


####Looking at ROH sharing across traits
#Loop through all traits and run comparisons within a trait
allTraits = list()
allTraits_Pval = list()

for (i in seq_along(fnames)) {
  trait = gsub(".*[/]([^.]+)[.txt].*", "\\1", fnames[i]) #keep track of trait and breed tested
  traitDF = read.delim(fnames[i])
  names(traitDF)[3] = "status"
  
  CountCaseControlperBreed = sharedROHIBD %>%
    filter(INDV1 %in% traitDF$dogID & INDV2 %in% traitDF$dogID ) %>%
    mutate(status_INDV1 = traitDF$status[match(INDV1, traitDF$dogID)],
           status_INDV2 = traitDF$status[match(INDV2, traitDF$dogID)],
           status = ifelse(status_INDV1 == status_INDV2, "same", "all"),
           status = ifelse(status_INDV1 == 1 & status == "same", "control", status),
           status = ifelse(status_INDV1 == 2 & status == "same", "case", status),
           Breed1 = popmapDryad$breed[match(INDV1, popmapDryad$dogID)],
           Breed2 = popmapDryad$breed[match(INDV2, popmapDryad$dogID)]) %>%
    distinct(status, status_INDV1, status_INDV2, Breed1, Breed2) %>%
    group_by(status, Breed1, Breed2) %>% 
    count() %>%
    mutate(numIndivs = n*2, #there are two individuals per comparison
           normConstant = (choose(2*as.numeric(numIndivs), 2)) - as.numeric(numIndivs)) %>%
    select(-c(n))
  
  caseControl = sharedROHIBD %>%
    filter(INDV1 %in% traitDF$dogID & INDV2 %in% traitDF$dogID ) %>%
    mutate(status_INDV1 = traitDF$status[match(INDV1, traitDF$dogID)],
           status_INDV2 = traitDF$status[match(INDV2, traitDF$dogID)],
           status = ifelse(status_INDV1 == status_INDV2, "same", "all"),
           status = ifelse(status_INDV1 == 1 & status == "same", "control", status),
           status = ifelse(status_INDV1 == 2 & status == "same", "case", status),
           Breed1 = popmapDryad$breed[match(INDV1, popmapDryad$dogID)],
           Breed2 = popmapDryad$breed[match(INDV2, popmapDryad$dogID)]) %>%
    group_by(status, Breed1, Breed2) %>% 
    summarise(GroupScore = sum(as.numeric(LenOverlap))) %>% #ROH sharing within IBD seg
    left_join(CountCaseControlperBreed) %>% #join on status, breed1, breed2
    mutate(NormGroupScorePerMb = (GroupScore/normConstant)/10^6,
           FinalStatus = case_when(
             (Breed1 == Breed2 & status == "control") ~ "within breed control", 
             (Breed1 != Breed2 & status == "control") ~ "between breed control",
             (Breed1 == Breed2 & status == "case") ~ "within breed case", 
             (Breed1 != Breed2 & status == "case") ~ "between breed case",
             (Breed1 == Breed2 & status == "all") ~ "within breed case vs control", 
             (Breed1 != Breed2 & status == "all") ~ "between breed case vs control"),
           trait = trait)
  
  countCaseSummaryTable = 
    caseControl %>%  
    group_by(FinalStatus) %>%
    summarise(count = n(),
              median = median(NormGroupScorePerMb),
              max = max(NormGroupScorePerMb),
              min = min(NormGroupScorePerMb))
  
  #run a pairwise wilcoxon-test on all comparisons
  Pvals_temp = pairwise.wilcox.test(caseControl$NormGroupScorePerMb, caseControl$FinalStatus,p.adjust.method = "BH")$p.value
  Pvals_df = data.frame(expand.grid(dimnames(Pvals_temp)),array(Pvals_temp)) %>%
    mutate(Trait = trait) %>%
    rename("BH_adj_pvalue" = "array.Pvals_temp.") %>%
    mutate(numBreedsVar1 = countCaseSummaryTable$count[match(Var1, countCaseSummaryTable$FinalStatus)],
           numBreedsVar2 = countCaseSummaryTable$count[match(Var2, countCaseSummaryTable$FinalStatus)],
           medianShareVar1 = countCaseSummaryTable$median[match(Var1, countCaseSummaryTable$FinalStatus)],
           medianShareVar2 = countCaseSummaryTable$median[match(Var2, countCaseSummaryTable$FinalStatus)],
           minShareVar1 = countCaseSummaryTable$min[match(Var1, countCaseSummaryTable$FinalStatus)],
           minShareVar2 = countCaseSummaryTable$min[match(Var2, countCaseSummaryTable$FinalStatus)],
           maxShareVar1 = countCaseSummaryTable$max[match(Var1, countCaseSummaryTable$FinalStatus)],
           maxShareVar2 = countCaseSummaryTable$max[match(Var2, countCaseSummaryTable$FinalStatus)])

  allTraits[[i]] = caseControl
  allTraits_Pval[[i]] = Pvals_df
  
  print(trait)
}

names(allTraits) = gsub(".*[/]([^.]+)[.txt].*", "\\1", fnames)
names(allTraits_Pval) = gsub(".*[/]([^.]+)[.txt].*", "\\1", fnames)


####Look at IBD Sharing across traits
allTraits_IBD = list()
allTraits_IBD_Pval = list()

for (i in seq_along(fnames)) {
  trait = gsub(".*[/]([^.]+)[.txt].*", "\\1", fnames[i]) #keep track of trait and breed tested
  traitDF = read.delim(fnames[i])
  names(traitDF)[3] = "status"
  
  CountCaseControlperBreed = sharedROHIBD %>%
    filter(INDV1 %in% traitDF$dogID & INDV2 %in% traitDF$dogID ) %>%
    mutate(status_INDV1 = traitDF$status[match(INDV1, traitDF$dogID)],
           status_INDV2 = traitDF$status[match(INDV2, traitDF$dogID)],
           status = ifelse(status_INDV1 == status_INDV2, "same", "all"),
           status = ifelse(status_INDV1 == 1 & status == "same", "control", status),
           status = ifelse(status_INDV1 == 2 & status == "same", "case", status),
           Breed1 = popmapDryad$breed[match(INDV1, popmapDryad$dogID)],
           Breed2 = popmapDryad$breed[match(INDV2, popmapDryad$dogID)]) %>%
    distinct(status, status_INDV1, status_INDV2, Breed1, Breed2) %>%
    group_by(status, Breed1, Breed2) %>% 
    count() %>%
    mutate(numIndivs = n*2, #there are two individuals per comparison
           normConstant = (choose(2*as.numeric(numIndivs), 2)) - as.numeric(numIndivs)) %>%
    select(-c(n))
  
  caseControl = sharedROHIBD %>%
    filter(INDV1 %in% traitDF$dogID & INDV2 %in% traitDF$dogID ) %>%
    mutate(status_INDV1 = traitDF$status[match(INDV1, traitDF$dogID)],
           status_INDV2 = traitDF$status[match(INDV2, traitDF$dogID)],
           status = ifelse(status_INDV1 == status_INDV2, "same", "all"),
           status = ifelse(status_INDV1 == 1 & status == "same", "control", status),
           status = ifelse(status_INDV1 == 2 & status == "same", "case", status),
           Breed1 = popmapDryad$breed[match(INDV1, popmapDryad$dogID)],
           Breed2 = popmapDryad$breed[match(INDV2, popmapDryad$dogID)]) %>%
    group_by(status, Breed1, Breed2) %>% 
    summarise(GroupScore = sum(as.numeric(NotCoveredROH))) %>% #IBD sharing outside of a ROH
    left_join(CountCaseControlperBreed) %>% #join on status, breed1, breed2
    mutate(NormGroupScorePerMb = (GroupScore/normConstant)/10^6,
           FinalStatus = case_when(
             (Breed1 == Breed2 & status == "control") ~ "within breed control", 
             (Breed1 != Breed2 & status == "control") ~ "between breed control",
             (Breed1 == Breed2 & status == "case") ~ "within breed case", 
             (Breed1 != Breed2 & status == "case") ~ "between breed case",
             (Breed1 == Breed2 & status == "all") ~ "within breed case vs control", 
             (Breed1 != Breed2 & status == "all") ~ "between breed case vs control"),
           trait = trait)
  
  countCaseSummaryTable = 
    caseControl %>%  
    group_by(FinalStatus) %>%
    summarise(count = n(),
              median = median(NormGroupScorePerMb),
              max = max(NormGroupScorePerMb),
              min = min(NormGroupScorePerMb))
  
  #run a pairwise wilcoxon-test on all comparisons
  Pvals_temp = pairwise.wilcox.test(caseControl$NormGroupScorePerMb, caseControl$FinalStatus,p.adjust.method = "BH")$p.value
  Pvals_df = data.frame(expand.grid(dimnames(Pvals_temp)),array(Pvals_temp)) %>%
    mutate(Trait = trait) %>%
    rename("BH_adj_pvalue" = "array.Pvals_temp.") %>%
    mutate(numBreedsVar1 = countCaseSummaryTable$count[match(Var1, countCaseSummaryTable$FinalStatus)],
           numBreedsVar2 = countCaseSummaryTable$count[match(Var2, countCaseSummaryTable$FinalStatus)],
           medianShareVar1 = countCaseSummaryTable$median[match(Var1, countCaseSummaryTable$FinalStatus)],
           medianShareVar2 = countCaseSummaryTable$median[match(Var2, countCaseSummaryTable$FinalStatus)],
           minShareVar1 = countCaseSummaryTable$min[match(Var1, countCaseSummaryTable$FinalStatus)],
           minShareVar2 = countCaseSummaryTable$min[match(Var2, countCaseSummaryTable$FinalStatus)],
           maxShareVar1 = countCaseSummaryTable$max[match(Var1, countCaseSummaryTable$FinalStatus)],
           maxShareVar2 = countCaseSummaryTable$max[match(Var2, countCaseSummaryTable$FinalStatus)])
  
  allTraits_IBD[[i]] = caseControl
  allTraits_IBD_Pval[[i]] = Pvals_df
  
  print(trait)
}

names(allTraits_IBD) = gsub(".*[/]([^.]+)[.txt].*", "\\1", fnames)
names(allTraits_IBD_Pval) = gsub(".*[/]([^.]+)[.txt].*", "\\1", fnames)

####Bind the lists for ROH and IBD sharing across traits
ROHSharing = bind_rows(allTraits)
ROHSharing_pvals = bind_rows(allTraits_Pval)
ROHSharing_pvals$Type = "ROH and IBD"
IBDSharing = bind_rows(allTraits_IBD)
IBDSharing_pvals = bind_rows(allTraits_IBD_Pval)
IBDSharing_pvals$Type = "IBD"

###aggregate across traits ROH
ROH_allTraits = bind_rows(allTraits)
ROH_allTraitsSummary = ROH_allTraits %>% 
  group_by(FinalStatus) %>%
  summarise(count = n(),
            median = median(GroupScore),
            max = max(GroupScore),
            min = min(GroupScore))

Pvals_temp = pairwise.wilcox.test(ROH_allTraits$GroupScore, ROH_allTraits$FinalStatus,p.adjust.method = "BH")$p.value

Pvals_ROH_allTraits = data.frame(expand.grid(dimnames(Pvals_temp)),array(Pvals_temp)) %>%
  rename("BH_adj_pvalue" = "array.Pvals_temp.")  %>% 
  mutate(Trait = "all traits",
         numBreedsVar1 = ROH_allTraitsSummary$count[match(Var1, ROH_allTraitsSummary$FinalStatus)],
         numBreedsVar2 = ROH_allTraitsSummary$count[match(Var2, ROH_allTraitsSummary$FinalStatus)],
         medianShareVar1 = ROH_allTraitsSummary$median[match(Var1, ROH_allTraitsSummary$FinalStatus)],
         medianShareVar2 = ROH_allTraitsSummary$median[match(Var2, ROH_allTraitsSummary$FinalStatus)],
         minShareVar1 = ROH_allTraitsSummary$min[match(Var1, ROH_allTraitsSummary$FinalStatus)],
         minShareVar2 = ROH_allTraitsSummary$min[match(Var2, ROH_allTraitsSummary$FinalStatus)],
         maxShareVar1 = ROH_allTraitsSummary$max[match(Var1, ROH_allTraitsSummary$FinalStatus)],
         maxShareVar2 = ROH_allTraitsSummary$max[match(Var2, ROH_allTraitsSummary$FinalStatus)],
         Type = "ROH and IBD")

###aggregate across traits ROH
IBD_allTraits = bind_rows(allTraits_IBD)
IBD_allTraitsSummary = IBD_allTraits %>% 
  group_by(FinalStatus) %>%
  summarise(count = n(),
            median = median(GroupScore),
            max = max(GroupScore),
            min = min(GroupScore))

Pvals_temp = pairwise.wilcox.test(IBD_allTraits$GroupScore, IBD_allTraits$FinalStatus,p.adjust.method = "BH")$p.value

Pvals_IBD_allTraits = data.frame(expand.grid(dimnames(Pvals_temp)),array(Pvals_temp)) %>%
  rename("BH_adj_pvalue" = "array.Pvals_temp.")  %>% 
  mutate(Trait = "all traits",
         numBreedsVar1 = IBD_allTraitsSummary$count[match(Var1, IBD_allTraitsSummary$FinalStatus)],
         numBreedsVar2 = IBD_allTraitsSummary$count[match(Var2, IBD_allTraitsSummary$FinalStatus)],
         medianShareVar1 = IBD_allTraitsSummary$median[match(Var1, IBD_allTraitsSummary$FinalStatus)],
         medianShareVar2 = IBD_allTraitsSummary$median[match(Var2, IBD_allTraitsSummary$FinalStatus)],
         minShareVar1 = IBD_allTraitsSummary$min[match(Var1, IBD_allTraitsSummary$FinalStatus)],
         minShareVar2 = IBD_allTraitsSummary$min[match(Var2, IBD_allTraitsSummary$FinalStatus)],
         maxShareVar1 = IBD_allTraitsSummary$max[match(Var1, IBD_allTraitsSummary$FinalStatus)],
         maxShareVar2 = IBD_allTraitsSummary$max[match(Var2, IBD_allTraitsSummary$FinalStatus)],
         Type = "IBD")

####Aggregate all the results
merged = rbind.data.frame(ROHSharing_pvals, Pvals_ROH_allTraits ,IBDSharing_pvals, Pvals_IBD_allTraits) %>% 
  mutate_if(is.numeric, round, digits=3)

#write.table(merged, "Pvalues_sharedROHwithIBD_sepByTrait_all.txt", row.names = F, col.names = T, sep = "\t", quote = F)

withinBreedOnly = merged %>% 
  filter(Var1 == "within breed control" & Var2 == "within breed case") %>%
  mutate(RangeVar1 = paste0("[", minShareVar1, "-", maxShareVar1, "]"), 
         RangeVar2 = paste0("[", minShareVar2, "-", maxShareVar2, "]")) %>%
  select(-c(minShareVar1, maxShareVar1,minShareVar2, maxShareVar2))

#write.table(withinBreedOnly, "Pvalues_sharedROHwithIBD_sepByTrait_withinBreeds.txt", row.names = F, col.names = T, sep = "\t", quote = F)