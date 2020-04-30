#load libraries
library(biomaRt)
library(tidyverse)

#Load files for hgnc symbols, dnds 
hgnc_symbols = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/VettingResults/hgnc_symbol_pairs.txt", col.names = c("previous", "current"), stringsAsFactors = F)

DNDS = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/VettingResults/primates_dn_ds_v69.txt", col.names = c("ensembl_gene_id","dn","ds"), stringsAsFactors = F) %>%
  mutate(dnds = dn/ds)

#Load biomart to get gene names
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=83, GRCh = 37)

genes = getBM(attributes = c("ensembl_gene_id", "transcript_length", "hgnc_symbol"), filters = "ensembl_gene_id", values = DNDS$ensembl_gene_id, mart = ensembl) %>%
  mutate(dnds = DNDS$dnds[match(ensembl_gene_id, DNDS$ensembl_gene_id)],
         data = "ensembl") %>%
  group_by(ensembl_gene_id) %>%
  filter(transcript_length == max(transcript_length)) %>%
  ungroup() %>%
  na.omit() %>%
  select(-c(transcript_length))
#rm(ensembl) #delete the mart

#Load file for exons without ROHs add hgnc and dnds info
exons = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/VettingResults/ExonRegion_NonOverlapsROH_cornellData.bed", stringsAsFactors = F) %>% 
  mutate(data = "no ROH") %>%
  left_join(hgnc_symbols, by = c("V4"="previous")) %>% #merge with hgnc symbols
  mutate(dnds = genes$dnds[match(current, genes$hgnc_symbol)], #add dn/ds
         data = "no ROH") %>%
  distinct(V4, .keep_all = T) %>% #keep one copy of each gene entry
  na.omit() %>%
  select(current,dnds,data)

#Grab morbid map data from OMIM
omimPLINK = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/VettingResults/MorbidMap_rmNonDiseasesProvisionalDiseases.txt") %>%
  separate_rows(geneName, sep = ",") %>% 
  filter(geneName %in% exons$current) %>%
  mutate(dnds = genes$dnds[match(geneName, genes$hgnc_symbol)],
         phenotype_MIM_INFO = str_replace_all(phenotype_MIM_INFO, "[[:digit:]]", ""),
         phenotype_MIM_INFO = str_replace_all(phenotype_MIM_INFO, "[[:punct:]]", " "))

#merge together to plot


