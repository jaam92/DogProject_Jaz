#Load libraries
library(tidyverse)
library(randomcoloR)

#Load files
setwd("~/Documents/DogProject_Jaz/LocalRscripts/ROH")
sharedROHIBD= read.delim("~/Documents/DogProject_Jaz/LocalRscripts/ROH/FinalSharedROHwithinIBDSegs_allChroms.bed", col.names = c("chrom","INDV1", "INDV2","segLen"))
popmapDryad = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/breeds_dryad.txt")
phenotypes = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/phenotypes.txt")
unrelateds = read.table("~/Documents/DogProject_Jaz/LocalRscripts/PCA_Unrelateds/UnrelatedIndividuals_allBreeds_mergedFitakCornell.txt")
fnames = paste0("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/splitPhenotypeFile/IncludeMixedBreeds/", 
                list.files(path="~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/splitPhenotypeFile/IncludeMixedBreeds/", 
                           pattern = "[^_]") %>%
                  str_subset(., "_", negate = TRUE)) #remove breed specific files

#Loop through all traits and run comparisons within a trait
allTraits = list()
allTraits_Pval = list()
allTraits_AvgSharing = list()

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
    summarise(GroupScore = sum(as.numeric(segLen))) %>%
    left_join(CountCaseControlperBreed) %>% #join on status, breed1, breed2
    mutate(NormGroupScore = GroupScore/normConstant,
           FinalStatus = case_when(
             (Breed1 == Breed2 & status == "control") ~ "within breed control", 
             (Breed1 != Breed2 & status == "control") ~ "between breed control",
             (Breed1 == Breed2 & status == "case") ~ "within breed case", 
             (Breed1 != Breed2 & status == "case") ~ "between breed case",
             (Breed1 == Breed2 & status == "all") ~ "within breed case vs control", 
             (Breed1 != Breed2 & status == "all") ~ "between breed case vs control"))
  
  #run a pairwise wilcoxon-test on all comparisons
  Pvals_temp = pairwise.wilcox.test(caseControl$NormGroupScore, caseControl$FinalStatus,p.adjust.method = "BH")$p.value
  Pvals_df = data.frame(expand.grid(dimnames(Pvals_temp)),array(Pvals_temp)) %>%
    mutate(trait = trait) %>%
    rename("BH_adj_pvalue" = "array.Pvals_temp.") %>%
    na.omit(BH_adj_pvalue) %>%
    filter(BH_adj_pvalue < 0.05)
  
  caseControlAvgSharing = caseControl %>%
    group_by(FinalStatus) %>%
    summarise(GroupScore = sum(as.numeric(NormGroupScore))/10^6) %>%
    mutate(trait = trait)

  allTraits[[i]] = caseControl
  allTraits_AvgSharing[[i]] = caseControlAvgSharing
  allTraits_Pval[[i]] = Pvals_df
  
  print(trait)
}

names(allTraits) = gsub(".*[/]([^.]+)[.txt].*", "\\1", fnames)
names(allTraits_AvgSharing) = gsub(".*[/]([^.]+)[.txt].*", "\\1", fnames)
names(allTraits_Pval) = gsub(".*[/]([^.]+)[.txt].*", "\\1", fnames)



#Looking at sharing across traits
ROHSharing = bind_rows(allTraits_AvgSharing)

meanROHSharing = ROHSharing %>%
  group_by(FinalStatus) %>%
  summarise(GroupScore = sum(as.numeric(GroupScore)))

Pvals_temp = pairwise.wilcox.test(ROHSharing$GroupScore, ROHSharing$FinalStatus,p.adjust.method = "BH")$p.value

Pvals_df = data.frame(expand.grid(dimnames(Pvals_temp)),array(Pvals_temp)) %>%
  rename("bonferroni_adj_pvalue" = "array.Pvals_temp.") %>%
  na.omit(bonferroni_adj_pvalues) %>%
  mutate(meanVar1 = mean$GroupScore[match(Var1, mean$FinalStatus)],
         meanVar2 = mean$GroupScore[match(Var2, mean$FinalStatus)])



#expand color palette of choice to hve number of colors equal to number of clades
colourCount_pop = length(unique(ROHSharing$FinalStatus)) 
palette = distinctColorPalette(colourCount_pop)

ggplot(ROHSharing, aes(x=GroupScore, y=FinalStatus, group=FinalStatus, colour=trait)) +
  geom_boxplot(show.legend = FALSE) +
  geom_point() + 
  scale_colour_manual(values = palette, na.value="grey") +
  theme_bw()  + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 10), 
        plot.title=element_text(size=26, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18)) +
  labs(x="Haplotype Sharing within ROH (Mb)", y="Group Comparison")