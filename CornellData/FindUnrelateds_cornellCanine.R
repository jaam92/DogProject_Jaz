. /u/local/Modules/default/init/modules.sh
module load R/3.4.2

#Load libraries
library(gdsfmt)
library(SNPRelate)
library(dplyr)

# PLINK BED files for unpruned data 
#bed.fn = ("cornell_canine_updatedFID.bed")
#bim.fn = ("cornell_canine_updatedFID.bim")
#fam.fn = ("cornell_canine_updatedFID.fam")

#convert
#snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "cornell_canine_updatedFID.gds")
#snpgdsSummary("cornell_canine_updatedFID.gds")

#Open file
genofile = snpgdsOpen("cornell_canine_updatedFID.gds")

#LD prune
snpset = snpgdsLDpruning(genofile, ld.threshold=0.8, maf = 0.05, missing.rate = 0.1, slide.max.n = 50, autosome.only = F)
snpset.id = unlist(snpset)


#IBD and kinship analysis with KING
#Must be done on LD pruned data

#Get sample ids and breed info from gds file
gdsSampIDs = read.gdsn(index.gdsn(genofile, "sample.id"))
breed.id = gsub(".*-","",gdsSampIDs)#remove everything before the dash to get breed

#run KING
ibd.robust = snpgdsIBDKING(genofile, sample.id=gdsSampIDs,family.id=breed.id,snp.id=snpset.id, num.thread=2)
KINGdf = snpgdsIBDSelection(ibd.robust)

#Identify relateds and remove them
RmRelated = unique(subset.data.frame(KINGdf, KINGdf$kinship >= 1/16)$ID1)

#Unrelated
Unrelateds = gdsSampIDs[!(gdsSampIDs %in% RmRelated)]
unrelated.ibd.robust = snpgdsIBDKING(genofile, sample.id=Unrelateds,family.id=gsub(".*-","",Unrelateds),snp.id=snpset.id, num.thread=2)
unrelatedKING = snpgdsIBDSelection(unrelated.ibd.robust)
#check that individuals areat most first cousins
summary(unrelatedKING$kinship)
dfUnrelateds = cbind.data.frame(Unrelateds, gsub(".*-","",Unrelateds))
names(dfUnrelateds)[2] = "breed"
dfUnrelateds$Unrelateds = gsub("-.*","",dfUnrelateds$Unrelateds) #fix sample id
#write.table(dfUnrelateds, "UnrelatedIndividuals_allBreeds.txt", sep = "\t", row.names = F, col.names = F, quote = F)
UnrelatedsPerBreed_n30 = dfUnrelateds %>% group_by(breed) %>% tally() %>% filter(n>=30) #find breeds with at least 30 unrelateds
Unrelated_sampsGrEql30 = dfUnrelateds %>% filter(breed %in% UnrelatedsPerBreed_n30$breed)
#write.table(Unrelated_sampsGrEql30, "UnrelatedIndividuals_grEql30.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#Close file
snpgdsClose(genofile)


