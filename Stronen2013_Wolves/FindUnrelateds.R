#Download packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("gdsfmt")
#biocLite("SNPRelate")

#Load libraries
library(gdsfmt)
library(SNPRelate)
library(dplyr)

# PLINK BED files for unpruned data aka random sample of 100000 snps
bed.fn = ("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/Stronen2013_Wolves/Stronen2013_CanFam3_MergedFitakCornellSharedSites.bed")
bim.fn = ("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/Stronen2013_Wolves/Stronen2013_CanFam3_MergedFitakCornellSharedSites.bim")
fam.fn = ("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/Stronen2013_Wolves/Stronen2013_CanFam3_MergedFitakCornellSharedSites.fam")

#convert
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/Stronen2013_Wolves/Stronen2013_CanFam3_MergedFitakCornellSharedSites.gds")
snpgdsSummary("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/Stronen2013_Wolves/Stronen2013_CanFam3_MergedFitakCornellSharedSites.gds")

#Open files
genofile = snpgdsOpen("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/Stronen2013_Wolves/Stronen2013_CanFam3_MergedFitakCornellSharedSites.gds")
popmap = read.delim("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/Stronen2013_Wolves/IndividualFiles/SampleInfo_Stronen2013_Wolves.txt", fill = T)

#LD prune
snpset = snpgdsLDpruning(genofile, ld.threshold=0.8, maf = 0.05, missing.rate = 0.1, slide.max.n = 50, autosome.only = F)
snpset.id = unlist(snpset)


#IBD and kinship analysis with KING
#Must be done on LD pruned data

#Get sample ids and breed info from gds file
gdsSampIDs = read.gdsn(index.gdsn(genofile, "sample.id"))
family.id = popmap$FID_Col1[match(gdsSampIDs, popmap$IID_Col2)]

#run KING
ibd.robust = snpgdsIBDKING(genofile, sample.id=gdsSampIDs,family.id=family.id, snp.id=snpset.id, num.thread=2, autosome.only = F)
KINGdf = snpgdsIBDSelection(ibd.robust)

#Identify relateds and remove them
RmRelated = unique(subset.data.frame(KINGdf, KINGdf$kinship >= 1/16)$ID1)

#Unrelated
Unrelateds = gdsSampIDs[!(gdsSampIDs %in% RmRelated)]
unrelated.ibd.robust = snpgdsIBDKING(genofile, sample.id=Unrelateds,family.id=popmap$FID_Col1[match(Unrelateds, popmap$IID_Col2)],snp.id=snpset.id, num.thread=2, autosome.only = F)
unrelatedKING = snpgdsIBDSelection(unrelated.ibd.robust)

#check that individuals areat most first cousins and write unrelateds to file
summary(unrelatedKING$kinship)
dfUnrelateds_Full = cbind.data.frame(Unrelateds, popmap$FID_Col1[match(Unrelateds, popmap$IID_Col2)],popmap$Population[match(Unrelateds, popmap$IID_Col2)], popmap$Cluster[match(Unrelateds, popmap$IID_Col2)])
names(dfUnrelateds_Full)[2] = "FID"
names(dfUnrelateds_Full)[3] = "Population"
names(dfUnrelateds_Full)[4] = "Cluster"

write.table(dfUnrelateds_Full, "Stronen2013_WolvesUnrelatedIndividuals_allBreeds.txt", sep = "\t", row.names = F, col.names = F, quote = F)

#Make a larger data frame with population and cluster info
UnrelatedsPerBreed_n30 = dfUnrelateds_Full %>% group_by(Cluster) %>% tally() %>% filter(n>=30) #find breeds with at least 30 unrelateds
UnrelatedsPerBreed_n30

Unrelated_sampsGrEql30 = dfUnrelateds_Full %>% filter(Cluster %in% UnrelatedsPerBreed_n30$Cluster) %>% select(FID,Unrelateds)#output in order of plink files
Unrelated_sampsGrEql30$breed = popmap$Cluster[match(Unrelated_sampsGrEql30$Unrelateds, popmap$IID_Col2)]#Save it as breed because this is colname for Rscript that parses ibd segments expects
write.table(Unrelated_sampsGrEql30, "Stronen2013_WolvesUnrelatedIndividuals_grEql30.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#Close file
snpgdsClose(genofile)
