#Download packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("gdsfmt")
#biocLite("SNPRelate")

#Load libraries
library(gdsfmt)
library(SNPRelate)
library(dplyr)

# PLINK BED files for unpruned data aka random sample of 100000 snps
bed.fn = ("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/MergeFitkakAndCornell/IntermidiateFile/MergedFile_CornellCanineFitak.bed")
bim.fn = ("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/MergeFitkakAndCornell/IntermidiateFile/MergedFile_CornellCanineFitak.bim")
fam.fn = ("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/MergeFitkakAndCornell/IntermidiateFile/MergedFile_CornellCanineFitak.fam")

#convert
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/MergeFitkakAndCornell/IntermidiateFile/MergedFile_CornellCanineFitak.gds")
snpgdsSummary("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/MergeFitkakAndCornell/IntermidiateFile/MergedFile_CornellCanineFitak.gds")

#Open file
genofile = snpgdsOpen("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/MergeFitkakAndCornell/IntermidiateFile/MergedFile_CornellCanineFitak.gds")

#LD prune
snpset = snpgdsLDpruning(genofile, ld.threshold=0.8, maf = 0.05, missing.rate = 0.1, slide.max.n = 50, autosome.only = F)
snpset.id = unlist(snpset)

#IBD and kinship analysis with KING
#Must be done on LD pruned data

#Get sample ids and breed info from gds file
gdsSampIDs = read.gdsn(index.gdsn(genofile, "sample.id"))
matchSamples = as.data.frame(gsub('(.*)-\\w+', '\\1',gdsSampIDs))
names(matchSamples)[1] = "sample.ids"
breed.id = gsub(".*-","",gdsSampIDs)#remove everything before the dash to get breed

#run KING
ibd.robust = snpgdsIBDKING(genofile, sample.id=gdsSampIDs,family.id=breed.id,snp.id=snpset.id, num.thread=2)
KINGdf = snpgdsIBDSelection(ibd.robust)

####Identify duplicates and relateds then remove them
Relateds = KINGdf %>% filter(kinship >= 1/16)
newSampIds  = gdsSampIDs[!(gdsSampIDs %in% unique(Relateds$ID1))]
unrelated.ibd.robust = snpgdsIBDKING(genofile, sample.id=newSampIds,family.id=gsub(".*-","",newSampIds),snp.id=snpset.id, num.thread=2)
unrelatedKING = snpgdsIBDSelection(unrelated.ibd.robust)

#check that individuals are at most first cousins
summary(unrelatedKING$kinship)
dfUnrelateds = cbind.data.frame(gsub('(.*)-\\w+', '\\1', newSampIds), gsub(".*-","",newSampIds))
names(dfUnrelateds)[1] = "Unrelateds"
names(dfUnrelateds)[2] = "breed"
write.table(dfUnrelateds, "UnrelatedIndividuals_allBreeds_mergedFitakCornell.txt", sep = "\t", row.names = F, col.names = F, quote = F)

####Only remove duplicate individuals
PotentialDuplicates = KINGdf[which(KINGdf$kinship > 0.49),]
id1_PotDups = unique(PotentialDuplicates$ID1)
id2_PotDups = unique(PotentialDuplicates$ID2)
mergePotDups = unique(c(id1_PotDups, id2_PotDups))
#manually picked individuals to keep that also matched with dfUnrelateds
keepIndivs = c("Box_LU132-boxer", "Box_LU134-boxer","MB-858-mexWolf","Pdl_GT332-poodle", "Pdl_GT333-poodle", "BoC_GT62-border_collie", "CKC_GT90-cavalier_king_charles_spaniel","CWD_GT107-czechoslovakian_wolf_dog", "EBT_GT147-english_bull_terrier", "GRe_GT218-golden_retriever", "PFZ43A03-irish_wolfhound","Ter_GT386-mix","Ter_GT387-mix","Wlf_LUb2-grayWolf", "WO_SoutheastAK_02zPOWIS-63-grayWolf")
excludeIndivs = mergePotDups[!(mergePotDups %in% keepIndivs)]
notDupsSampIds  = gdsSampIDs[!(gdsSampIDs %in% excludeIndivs)]
rmDups.ibd.robust = snpgdsIBDKING(genofile, sample.id=notDupsSampIds,family.id=gsub(".*-","",notDupsSampIds),snp.id=snpset.id, num.thread=2)
rmDupsKING = snpgdsIBDSelection(rmDups.ibd.robust)
dfRmDups = cbind.data.frame(gsub('(.*)-\\w+', '\\1', notDupsSampIds), gsub(".*-","",notDupsSampIds))
names(dfRmDups)[1] = "dogIDs"
names(dfRmDups)[2] = "breed"
write.table(dfRmDups, "Individuals_allBreeds_mergedFitakCornell.txt", sep = "\t", row.names = F, col.names = F, quote = F)

#ID groups with at least 50 unrelateds
wolfData = read.table("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/IndividualFiles/canine-cluster2.txt")
UnrelatedsPerBreed_n50 = dfUnrelateds %>% group_by(breed) %>% tally() %>% filter(n>=50) #find breeds with at least 50 unrelateds
Unrelated_sampsGrEql50 = dfUnrelateds %>% filter(breed %in% UnrelatedsPerBreed_n50$breed)
Unrelated_sampsGrEql50$WolfType = wolfData$V3[match(Unrelated_sampsGrEql50$Unrelateds, wolfData$V2)]#all the NA are dogs
Unrelated_sampsGrEql50 %>% group_by(WolfType) %>% tally() #only european breeds have at least 50 so use those
#subset out dogs and european wolves
FinalUnrelatedDF = Unrelated_sampsGrEql50 %>% filter(is.na(WolfType) | WolfType == "EURO") %>% select(Unrelateds, breed) %>% as.data.frame()
write.table(FinalUnrelatedDF, "UnrelatedIndividuals_grEql50_MergedFitakCornell.txt", sep = "\t", row.names = F, col.names = F, quote = F)

#Close file
snpgdsClose(genofile)
