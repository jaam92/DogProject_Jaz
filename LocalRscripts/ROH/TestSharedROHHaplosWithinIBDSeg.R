#Load libraries
library(tidyverse)
library(data.table)
library(mgsub)
library(cowplot)
library(ggpubr)

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
  
  caseControl = sharedROHIBD %>%
    filter(INDV1 %in% traitDF$dogID & INDV2 %in% traitDF$dogID ) %>%
    mutate(status_INDV1 = traitDF$status[match(INDV1, traitDF$dogID)],
           status_INDV2 = traitDF$status[match(INDV2, traitDF$dogID)],
           status = ifelse(status_INDV1 == status_INDV2, "same", "all"),
           status = ifelse(status_INDV1 == 1 & status == "same", "control", status),
           status = ifelse(status_INDV1 == 2 & status == "same", "case", status),
           Breed1 = popmapDryad$breed[match(INDV1, popmapDryad$dogID)],
           Breed2 = popmapDryad$breed[match(INDV2, popmapDryad$dogID)],
           key = paste(pmin(Breed1, Breed2), pmax(Breed1, Breed2), sep = "-")) %>%
    group_by(status, key) %>%
    summarise(GroupScore = sum(as.numeric(LenOverlap))) %>% #sum up shared ROH within IBD for group
    ungroup() %>%
    separate(col = key, into = c("Breed1", "Breed2"), sep = "-") %>% #get breed info back
    mutate(FinalStatus = case_when(
      (Breed1 == Breed2 & status == "control") ~ "within breed control", 
      (Breed1 != Breed2 & status == "control") ~ "between breed control",
      (Breed1 == Breed2 & status == "case") ~ "within breed case", 
      (Breed1 != Breed2 & status == "case") ~ "between breed case",
      (Breed1 == Breed2 & status == "all") ~ "within breed control vs case", 
      (Breed1 != Breed2 & status == "all") ~ "between breed control vs case"),
      trait = trait) %>%
    filter(status != "all") %>% 
    mutate(sampSizeB1 = caseControlCountsBreed$n[match(Breed1, caseControlCountsBreed$Breed1)],
           sampSizeB2 = caseControlCountsBreed$n[match(Breed2, caseControlCountsBreed$Breed1)],
           finalSamp = ifelse(Breed1==Breed2, sampSizeB1, sampSizeB1+sampSizeB2)) %>% #if same breed don't need to add sample sizes together
    filter(finalSamp > 1) %>% #number of individuals in a breed has to be greater than 1 (it should be but just in case)
    mutate(normConstant = (choose(2*as.numeric(finalSamp), 2)) - as.numeric(finalSamp),
           NormGroupScorePerMb = (GroupScore/normConstant)/10^6)
  
  SummaryTable = caseControl %>%  
    group_by(FinalStatus) %>%
    summarise(numBreeds = n(),
              numSamps = sum(finalSamp),
              median = median(NormGroupScorePerMb),
              max = max(NormGroupScorePerMb),
              min = min(NormGroupScorePerMb))
  
  #run a pairwise wilcoxon-test on all comparisons
  Pvals_temp = pairwise.wilcox.test(caseControl$NormGroupScorePerMb, caseControl$FinalStatus, p.adjust.method="none")$p.value
  Pvals_df = data.frame(expand.grid(dimnames(Pvals_temp)),array(Pvals_temp)) %>%
    mutate(Trait = trait) %>%
    rename("pvalue" = "array.Pvals_temp.") %>%
    mutate(numIndivsVar1 = SummaryTable$numSamps[match(Var1, SummaryTable$FinalStatus)],
           numBreedsVar1 = SummaryTable$numBreeds[match(Var1, SummaryTable$FinalStatus)],
           medianShareVar1 = SummaryTable$median[match(Var1, SummaryTable$FinalStatus)],
           minShareVar1 = SummaryTable$min[match(Var1, SummaryTable$FinalStatus)],
           maxShareVar2 = SummaryTable$max[match(Var2, SummaryTable$FinalStatus)],
           numIndivsVar2 = SummaryTable$numSamps[match(Var2, SummaryTable$FinalStatus)],
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
  
  caseControl = sharedROHIBD %>%
    filter(INDV1 %in% traitDF$dogID & INDV2 %in% traitDF$dogID ) %>%
    mutate(status_INDV1 = traitDF$status[match(INDV1, traitDF$dogID)],
           status_INDV2 = traitDF$status[match(INDV2, traitDF$dogID)],
           status = ifelse(status_INDV1 == status_INDV2, "same", "all"),
           status = ifelse(status_INDV1 == 1 & status == "same", "control", status),
           status = ifelse(status_INDV1 == 2 & status == "same", "case", status),
           Breed1 = popmapDryad$breed[match(INDV1, popmapDryad$dogID)],
           Breed2 = popmapDryad$breed[match(INDV2, popmapDryad$dogID)],
           key = paste(pmin(Breed1, Breed2), pmax(Breed1, Breed2), sep = "-")) %>%
    group_by(status, key) %>%
    summarise(GroupScore = sum(as.numeric(NotCoveredROH))) %>% #sum up shared IBD for group
    ungroup() %>%
    separate(col = key, into = c("Breed1", "Breed2"), sep = "-") %>% #get breed info back 
    mutate(FinalStatus = case_when(
      (Breed1 == Breed2 & status == "control") ~ "within breed control", 
      (Breed1 != Breed2 & status == "control") ~ "between breed control",
      (Breed1 == Breed2 & status == "case") ~ "within breed case", 
      (Breed1 != Breed2 & status == "case") ~ "between breed case",
      (Breed1 == Breed2 & status == "all") ~ "within breed control vs case", 
      (Breed1 != Breed2 & status == "all") ~ "between breed control vs case"),
      trait = trait) %>%
    filter(status != "all") %>% 
    mutate(sampSizeB1 = caseControlCountsBreed$n[match(Breed1, caseControlCountsBreed$Breed1)],
           sampSizeB2 = caseControlCountsBreed$n[match(Breed2, caseControlCountsBreed$Breed1)],
           finalSamp = ifelse(Breed1==Breed2, sampSizeB1, sampSizeB1+sampSizeB2)) %>% #if same breed don't need to add sample sizes together
    filter(finalSamp > 1) %>% #number of individuals in a breed has to be greater than 1 (it should be but just in case)
    mutate(normConstant = (choose(2*as.numeric(finalSamp), 2)) - as.numeric(finalSamp),
           NormGroupScorePerMb = (GroupScore/normConstant)/10^6)
  
  SummaryTable = caseControl %>%  
    group_by(FinalStatus) %>%
    summarise(numBreeds = n(),
              numSamps = sum(finalSamp),
              median = median(NormGroupScorePerMb),
              max = max(NormGroupScorePerMb),
              min = min(NormGroupScorePerMb))
  
  #run a pairwise wilcoxon-test on all comparisons
  Pvals_temp = pairwise.wilcox.test(caseControl$NormGroupScorePerMb, caseControl$FinalStatus, p.adjust.method="none")$p.value
  Pvals_df = data.frame(expand.grid(dimnames(Pvals_temp)),array(Pvals_temp)) %>%
    mutate(Trait = trait) %>%
    rename("pvalue" = "array.Pvals_temp.") %>%
    mutate(numIndivsVar1 = SummaryTable$numSamps[match(Var1, SummaryTable$FinalStatus)],
           numBreedsVar1 = SummaryTable$numBreeds[match(Var1, SummaryTable$FinalStatus)],
           medianShareVar1 = SummaryTable$median[match(Var1, SummaryTable$FinalStatus)],
           minShareVar1 = SummaryTable$min[match(Var1, SummaryTable$FinalStatus)],
           maxShareVar2 = SummaryTable$max[match(Var2, SummaryTable$FinalStatus)],
           numIndivsVar2 = SummaryTable$numSamps[match(Var2, SummaryTable$FinalStatus)],
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
ROHSharing_pvals = bind_rows(allTraits_ROH_Pval) %>%
  na.omit() #remove pairwise comparisons with same group
ROHSharing_pvals$Type = "ROH and IBD"
IBDSharing = bind_rows(allTraits_IBD)
IBDSharing_pvals = bind_rows(allTraits_IBD_Pval) %>%
  na.omit() #remove pairwise comparisons with same group
IBDSharing_pvals$Type = "IBD"

###aggregate across traits ROH
ROH_allTraitsSummary = ROHSharing %>% 
  group_by(FinalStatus) %>%
  summarise(count = n(),
            median = median(NormGroupScorePerMb),
            max = max(NormGroupScorePerMb),
            min = min(NormGroupScorePerMb),
            numIndivs = sum(finalSamp))

Pvals_temp = pairwise.wilcox.test(ROHSharing$NormGroupScorePerMb, ROHSharing$FinalStatus, p.adjust.method="none")$p.value

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
IBD_allTraitsSummary = IBDSharing %>% 
  group_by(FinalStatus) %>%
  summarise(count = n(),
            median = median(NormGroupScorePerMb),
            max = max(NormGroupScorePerMb),
            min = min(NormGroupScorePerMb),
            numIndivs = sum(finalSamp))

Pvals_temp = pairwise.wilcox.test(IBDSharing$NormGroupScorePerMb, IBDSharing$FinalStatus, p.adjust.method="none")$p.value

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
  select(-c(minShareVar1, maxShareVar1,minShareVar2, maxShareVar2)) %>%
  na.omit()  #remove comparisons between same group
  
withinBreedComps = merged %>% 
  filter(!str_detect(Var1, "between*") & !str_detect(Var2, "between*"))

#write.table(withinBreedComps, "Pvalues_sharedROHwithIBD_sepByTrait_withinBreeds.txt", row.names = F, col.names = T, sep = "\t", quote = F)

####Plots
cbPalette = c("All traits" = "gray25", "CLLD" = "#D55E00",  "PSVA" = "steelblue", "lymphoma" = "#009E73", "MCT" = "gold3", "MVD" = "mediumpurple4", "GC" = "#CC79A7", "ED" = "#867BCF")

ROH = ggplot(ROHSharing, aes(y=FinalStatus, x=NormGroupScorePerMb)) +
  geom_boxplot(outlier.shape = NA) + #remove outlier points and only use jitter
  geom_jitter(aes(colour=trait), size = 3) +
  labs(x="Pairwise sharing within\nROH and IBD segement (Mb)", y="Comparison group") +
  scale_colour_manual(name = "Trait", values = cbPalette) + 
  theme_bw() +
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        plot.title=element_text(size=24), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))

ROH_inset = ROHSharing %>%
  filter(str_detect(FinalStatus, "between*")) %>%
  ggplot(., aes(y=FinalStatus, x=NormGroupScorePerMb)) +
  geom_boxplot(outlier.shape = NA) + #remove outlier points and only use jitter
  geom_jitter(aes(colour=trait), size = 3) +
  labs(x="Pairwise sharing within\nROH and IBD segement (Mb)", y="Comparison group") +
  scale_colour_manual(name = "Trait", values = cbPalette) + 
  theme_bw() +
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        plot.title=element_blank(), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none")

plotROH = ROH + 
  annotation_custom(grob=ggplotGrob(ROH_inset), xmin = 500, xmax=1500, ymin = "between breed case", ymax="between breed control")

IBD_inset = IBDSharing %>%
  filter(str_detect(FinalStatus, "between*")) %>%
  ggplot(., aes(y=FinalStatus, x=NormGroupScorePerMb)) +
  geom_boxplot(outlier.shape = NA) + #remove outlier points and only use jitter
  geom_jitter(aes(colour=trait), size = 3) +
  labs(x="Pairwise sharing outside of ROH\nwithin IBD segement (Mb)", y="Comparison group") +
  scale_colour_manual(name = "Trait", values = cbPalette) + 
  theme_bw() +
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_text(size  = 20), 
        plot.title=element_blank(), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none")

IBD = ggplot(IBDSharing, aes(y=FinalStatus, x=NormGroupScorePerMb)) +
  geom_boxplot(outlier.shape = NA) + #remove outlier points and only use jitter
  geom_jitter(aes(colour=trait), size = 3) +
  labs(x="Pairwise sharing outside of ROH\nwithin IBD segement (Mb)", y="Comparison group") +
  scale_colour_manual(name = "Trait", values = cbPalette) + 
  theme_bw() +
  theme(axis.text.x = element_text(size  = 20), 
        axis.text.y = element_blank(), 
        plot.title=element_text(size=24), 
        axis.title.x=element_text(size=24),
        axis.title.y=element_blank(),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))


plotIBD = IBD + 
  annotation_custom(grob=ggplotGrob(IBD_inset), xmin = 5000, xmax=15000, ymin = "between breed case", ymax="between breed control")

#pdf(file = "HaplotypeSharingInIBDandROH.pdf", height = 10, width = 24)
print(ggarrange(plotROH, 
          plotIBD, 
          common.legend =TRUE, 
          legend = "right",
          align = "h"))
#dev.off()
