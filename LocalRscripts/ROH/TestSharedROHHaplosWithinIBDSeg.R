#Load libraries
library(tidyverse)
library(data.table)
library(mgsub)

#Empty lists
allTraits_ROH = list()
allTraits_ROH_Pval = list()
allTraits_IBD = list()
allTraits_IBD_Pval = list()

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


####Looking at ROH sharing across traits
for (i in seq_along(fnames)) {
  trait = gsub(".*[/]([^.]+)[.txt].*", "\\1", fnames[i]) #keep track of trait and breed tested
  traitDF = read.delim(fnames[i])
  names(traitDF)[3] = "status"
  
  caseControlCountsBreed = traitDF %>% 
    filter(dogID %in% sharedROHIBD$INDV1 | dogID %in% sharedROHIBD$INDV2) %>%
    mutate(status = ifelse(status == 1, "control", "case")) %>%
    rename("Breed1" = "breed") %>%
    group_by(Breed1, status) %>%
    count() 
  
  caseControlCounts = traitDF %>% 
    filter(dogID %in% sharedROHIBD$INDV1 | dogID %in% sharedROHIBD$INDV2) %>%
    mutate(status = ifelse(status == 1, "within breed control", "within breed case")) %>% #rename now to make easier to attach to pvalues data frame
    rename("Breed1" = "breed") %>%
    group_by(status) %>%
    count()  
  
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
    summarise(GroupScore = sum(as.numeric(LenOverlap))) %>% #sum up shared ROH within IBD for group
    ungroup() %>%
    mutate(FinalStatus = case_when(
      (Breed1 == Breed2 & status == "control") ~ "within breed control", 
      (Breed1 != Breed2 & status == "control") ~ "between breed control",
      (Breed1 == Breed2 & status == "case") ~ "within breed case", 
      (Breed1 != Breed2 & status == "case") ~ "between breed case",
      (Breed1 == Breed2 & status == "all") ~ "within breed control vs case", 
      (Breed1 != Breed2 & status == "all") ~ "between breed control vs case"),
      trait = trait) %>%
    filter(FinalStatus == "within breed control" | FinalStatus == "within breed case") %>%
    left_join(caseControlCountsBreed) %>% #join on status, breed1, breed2
    filter(n > 1) %>% #number of individuals in a breed has to be greater than 1 (it should be but just in case)
    mutate(normConstant = (choose(2*as.numeric(n), 2)) - as.numeric(n),
           NormGroupScorePerMb = (GroupScore/normConstant)/10^6)
  
  
  SummaryTable = caseControl %>%  
    group_by(FinalStatus) %>%
    summarise(numBreeds = n(),
              median = median(NormGroupScorePerMb),
              max = max(NormGroupScorePerMb),
              min = min(NormGroupScorePerMb))
  
  #run a pairwise wilcoxon-test on all comparisons
  Pvals_temp = pairwise.wilcox.test(caseControl$NormGroupScorePerMb, caseControl$FinalStatus, p.adjust.method="none")$p.value
  Pvals_df = data.frame(expand.grid(dimnames(Pvals_temp)),array(Pvals_temp)) %>%
    mutate(Trait = trait) %>%
    rename("pvalue" = "array.Pvals_temp.") %>%
    mutate(numIndivsVar1 = caseControlCounts$n[match(Var1, caseControlCounts$status)],
           numBreedsVar1 = SummaryTable$numBreeds[match(Var1, SummaryTable$FinalStatus)],
           medianShareVar1 = SummaryTable$median[match(Var1, SummaryTable$FinalStatus)],
           minShareVar1 = SummaryTable$min[match(Var1, SummaryTable$FinalStatus)],
           maxShareVar2 = SummaryTable$max[match(Var2, SummaryTable$FinalStatus)],
           numIndivsVar2 = caseControlCounts$n[match(Var2, caseControlCounts$status)],
           numBreedsVar2 = SummaryTable$numBreeds[match(Var2, SummaryTable$FinalStatus)],
           medianShareVar2 = SummaryTable$median[match(Var2, SummaryTable$FinalStatus)],
           minShareVar2 = SummaryTable$min[match(Var2, SummaryTable$FinalStatus)],
           maxShareVar1 = SummaryTable$max[match(Var1, SummaryTable$FinalStatus)]) 
  
  allTraits_ROH[[i]] = caseControl
  allTraits_ROH_Pval[[i]] = Pvals_df
  
  print(trait)
}

####Looking at IBD sharing across traits
for (i in seq_along(fnames)) {
  trait = gsub(".*[/]([^.]+)[.txt].*", "\\1", fnames[i]) #keep track of trait and breed tested
  traitDF = read.delim(fnames[i])
  names(traitDF)[3] = "status"
  
  caseControlCountsBreed = traitDF %>% 
    filter(dogID %in% sharedROHIBD$INDV1 | dogID %in% sharedROHIBD$INDV2) %>%
    mutate(status = ifelse(status == 1, "control", "case")) %>%
    rename("Breed1" = "breed") %>%
    group_by(Breed1, status) %>%
    count() 
  
  caseControlCounts = traitDF %>% 
    filter(dogID %in% sharedROHIBD$INDV1 | dogID %in% sharedROHIBD$INDV2) %>%
    mutate(status = ifelse(status == 1, "within breed control", "within breed case")) %>% #rename now to make easier to attach to pvalues data frame
    rename("Breed1" = "breed") %>%
    group_by(status) %>%
    count()  
  
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
    summarise(GroupScore = sum(as.numeric(NotCoveredROH))) %>% #sum up IBD not covered by ROH for group
    ungroup() %>%
    mutate(FinalStatus = case_when(
      (Breed1 == Breed2 & status == "control") ~ "within breed control", 
      (Breed1 != Breed2 & status == "control") ~ "between breed control",
      (Breed1 == Breed2 & status == "case") ~ "within breed case", 
      (Breed1 != Breed2 & status == "case") ~ "between breed case",
      (Breed1 == Breed2 & status == "all") ~ "within breed control vs case", 
      (Breed1 != Breed2 & status == "all") ~ "between breed control vs case"),
      trait = trait) %>%
    filter(FinalStatus == "within breed control" | FinalStatus == "within breed case") %>%
    left_join(caseControlCountsBreed) %>% #join on status, breed1, breed2
    filter(n > 1) %>% #number of individuals in a breed has to be greater than 1 (it should be but just in case)
    mutate(normConstant = (choose(2*as.numeric(n), 2)) - as.numeric(n),
           NormGroupScorePerMb = (GroupScore/normConstant)/10^6)
  
  
  SummaryTable = caseControl %>%  
    group_by(FinalStatus) %>%
    summarise(numBreeds = n(),
              median = median(NormGroupScorePerMb),
              max = max(NormGroupScorePerMb),
              min = min(NormGroupScorePerMb))
  
  #run a pairwise wilcoxon-test on all comparisons
  Pvals_temp = pairwise.wilcox.test(caseControl$NormGroupScorePerMb, caseControl$FinalStatus, p.adjust.method="none")$p.value
  Pvals_df = data.frame(expand.grid(dimnames(Pvals_temp)),array(Pvals_temp)) %>%
    mutate(Trait = trait) %>%
    rename("pvalue" = "array.Pvals_temp.") %>%
    mutate(numIndivsVar1 = caseControlCounts$n[match(Var1, caseControlCounts$status)],
           numBreedsVar1 = SummaryTable$numBreeds[match(Var1, SummaryTable$FinalStatus)],
           medianShareVar1 = SummaryTable$median[match(Var1, SummaryTable$FinalStatus)],
           minShareVar1 = SummaryTable$min[match(Var1, SummaryTable$FinalStatus)],
           maxShareVar2 = SummaryTable$max[match(Var2, SummaryTable$FinalStatus)],
           numIndivsVar2 = caseControlCounts$n[match(Var2, caseControlCounts$status)],
           numBreedsVar2 = SummaryTable$numBreeds[match(Var2, SummaryTable$FinalStatus)],
           medianShareVar2 = SummaryTable$median[match(Var2, SummaryTable$FinalStatus)],
           minShareVar2 = SummaryTable$min[match(Var2, SummaryTable$FinalStatus)],
           maxShareVar1 = SummaryTable$max[match(Var1, SummaryTable$FinalStatus)]) 
  
  allTraits_IBD[[i]] = caseControl
  allTraits_IBD_Pval[[i]] = Pvals_df
  
  print(trait)
}

####Name lists
names(allTraits_ROH) = gsub(".*[/]([^.]+)[.txt].*", "\\1", fnames)
names(allTraits_ROH_Pval) = gsub(".*[/]([^.]+)[.txt].*", "\\1", fnames) 
names(allTraits_IBD) = gsub(".*[/]([^.]+)[.txt].*", "\\1", fnames)
names(allTraits_IBD_Pval) = gsub(".*[/]([^.]+)[.txt].*", "\\1", fnames)

####Bind the lists for ROH and IBD sharing across traits
ROHSharing = bind_rows(allTraits_ROH)
ROHSharing_pvals = bind_rows(allTraits_ROH_Pval)
ROHSharing_pvals$Type = "ROH and IBD"
IBDSharing = bind_rows(allTraits_IBD)
IBDSharing_pvals = bind_rows(allTraits_IBD_Pval)
IBDSharing_pvals$Type = "IBD"

###aggregate across traits ROH
ROH_allTraits = bind_rows(allTraits_ROH)
ROH_allTraitsSummary = ROH_allTraits %>% 
  group_by(FinalStatus) %>%
  summarise(count = n(),
            median = median(NormGroupScorePerMb),
            max = max(NormGroupScorePerMb),
            min = min(NormGroupScorePerMb),
            numIndivs = sum(n))

Pvals_temp = pairwise.wilcox.test(ROH_allTraits$NormGroupScorePerMb, ROH_allTraits$FinalStatus, p.adjust.method="none")$p.value

Pvals_ROH_allTraits = data.frame(expand.grid(dimnames(Pvals_temp)),array(Pvals_temp)) %>%
  rename("pvalue" = "array.Pvals_temp.")  %>% 
  mutate(Trait = "All traits",
         numIndivsVar1 = ROH_allTraitsSummary$numIndivs[match(Var1, ROH_allTraitsSummary$FinalStatus)],
         numBreedsVar1 = ROH_allTraitsSummary$count[match(Var1, ROH_allTraitsSummary$FinalStatus)],
         medianShareVar1 = ROH_allTraitsSummary$median[match(Var1, ROH_allTraitsSummary$FinalStatus)],
         minShareVar1 = ROH_allTraitsSummary$min[match(Var1, ROH_allTraitsSummary$FinalStatus)],
         maxShareVar1 = ROH_allTraitsSummary$max[match(Var1, ROH_allTraitsSummary$FinalStatus)],
         numIndivsVar2 = ROH_allTraitsSummary$numIndivs[match(Var2, ROH_allTraitsSummary$FinalStatus)],
         numBreedsVar2 = ROH_allTraitsSummary$count[match(Var2, ROH_allTraitsSummary$FinalStatus)],
         medianShareVar2 = ROH_allTraitsSummary$median[match(Var2, ROH_allTraitsSummary$FinalStatus)],
         minShareVar2 = ROH_allTraitsSummary$min[match(Var2, ROH_allTraitsSummary$FinalStatus)],
         maxShareVar2 = ROH_allTraitsSummary$max[match(Var2, ROH_allTraitsSummary$FinalStatus)],
         Type = "ROH and IBD")

###aggregate across traits IBD
IBD_allTraits = bind_rows(allTraits_IBD)

IBD_allTraitsSummary = IBD_allTraits %>% 
  group_by(FinalStatus) %>%
  summarise(count = n(),
            median = median(NormGroupScorePerMb),
            max = max(NormGroupScorePerMb),
            min = min(NormGroupScorePerMb),
            numIndivs = sum(n))

Pvals_temp = pairwise.wilcox.test(IBD_allTraits$NormGroupScorePerMb, IBD_allTraits$FinalStatus, p.adjust.method="none")$p.value

Pvals_IBD_allTraits = data.frame(expand.grid(dimnames(Pvals_temp)),array(Pvals_temp)) %>%
  rename("pvalue" = "array.Pvals_temp.")  %>% 
  mutate(Trait = "All traits",
         numIndivsVar1 = ROH_allTraitsSummary$numIndivs[match(Var1, ROH_allTraitsSummary$FinalStatus)],
         numBreedsVar1 = IBD_allTraitsSummary$count[match(Var1, IBD_allTraitsSummary$FinalStatus)],
         medianShareVar1 = IBD_allTraitsSummary$median[match(Var1, IBD_allTraitsSummary$FinalStatus)],
         minShareVar1 = IBD_allTraitsSummary$min[match(Var1, IBD_allTraitsSummary$FinalStatus)],
         maxShareVar1 = IBD_allTraitsSummary$max[match(Var1, IBD_allTraitsSummary$FinalStatus)],
         numIndivsVar2 = ROH_allTraitsSummary$numIndivs[match(Var2, ROH_allTraitsSummary$FinalStatus)],
         numBreedsVar2 = IBD_allTraitsSummary$count[match(Var2, IBD_allTraitsSummary$FinalStatus)],
         medianShareVar2 = IBD_allTraitsSummary$median[match(Var2, IBD_allTraitsSummary$FinalStatus)],
         minShareVar2 = IBD_allTraitsSummary$min[match(Var2, IBD_allTraitsSummary$FinalStatus)],
         maxShareVar2 = IBD_allTraitsSummary$max[match(Var2, IBD_allTraitsSummary$FinalStatus)],
         Type = "IBD")

####Aggregate all the results
merged = rbind.data.frame(ROHSharing_pvals, Pvals_ROH_allTraits ,IBDSharing_pvals, Pvals_IBD_allTraits) %>% 
  mutate_if(is.numeric, round, digits=3) %>%
  mutate(RangeVar1 = paste0("[", minShareVar1, "-", maxShareVar1, "]"), 
         RangeVar2 = paste0("[", minShareVar2, "-", maxShareVar2, "]")) %>%
  select(-c(minShareVar1, maxShareVar1,minShareVar2, maxShareVar2))

#write.table(merged, "Pvalues_sharedROHwithIBD_sepByTrait_withinBreeds.txt", row.names = F, col.names = T, sep = "\t", quote = F)
