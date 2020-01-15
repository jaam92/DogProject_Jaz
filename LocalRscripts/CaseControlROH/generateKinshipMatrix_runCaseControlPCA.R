#Load libraries
library(tidyverse)
library(gdsfmt)
library(SNPRelate)
library(GENESIS)
library(GWASTools)

####Make gds file
#bed.fn = ("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/MergedFile_CornellCanineFitak_allIndivs.bed")
#bim.fn = ("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/MergedFile_CornellCanineFitak_allIndivs.bim")
#fam.fn = ("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/MergedFile_CornellCanineFitak_allIndivs.fam")

#convert
#snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/MergedFile_CornellCanineFitak_allIndivs.gds")
#snpgdsSummary("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/MergedFile_CornellCanineFitak_allIndivs.gds")

#open the new gds file
genofile = snpgdsOpen("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/MergedFile_CornellCanineFitak_allIndivs.gds")
sampIds = read.gdsn(index.gdsn(genofile, "sample.id")) #grab sample ids 
famIds = gsub(".*-","",sampIds) #make family ids

#LD prune because data set is small, set r2 = 0.7, maf 5%, and 50 snp window 
snpset = snpgdsLDpruning(genofile, sample.id=sampIds, method="corr", slide.max.n=50, ld.threshold=0.3, maf = 0.05, autosome.only = F) 
snpset.id = unlist(snpset)

###Make a kinship matrix that accounts for known relatedness of samples
#Run king
king = snpgdsIBDKING(genofile, snp.id=snpset.id)
kingMat = king$kinship
colnames(kingMat) <- rownames(kingMat) <- king$sample.id
snpgdsClose(genofile)#Close file

#Re-Open GDS data 
genoFile = GdsGenotypeReader(filename = "~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/MergedFile_CornellCanineFitak_allIndivs.gds")#read in GDS data
genoData = GenotypeData(genoFile)#create a GenotypeData class object

##Partition data into relateds and unrelated (less than first cousins)
set.seed(100)
sampset = pcairPartition(kinobj = kingMat, kin.thresh=2^(-9/2), divobj = kingMat, div.thresh=2^(-9/2)) 
unrelsID = sampset$unrels

#Run PC-AiR
canidpcAir = pcair(genoData, kinobj = kingMat, kin.thresh=2^(-9/2), divobj = kingMat, div.thresh=2^(-9/2), unrel.set = unrelsID, snp.include = snpset.id, autosome.only = F)

#Prep GDS for PC-Relate
genoData = GenotypeBlockIterator(genoData, snpBlock = 20000)#take input from PC-AiR and convert to Genotype block iterator

#Run PC-Relate
canidpcRelate = pcrelate(genoData, pcs = canidpcAir$vectors[,1:2], training.set = canidpcAir$unrels, ibd.probs = FALSE) #use first 2 pcs to correct kinship for population structure (aka ancestry)
pcRelateMat = pcrelateToMatrix(canidpcRelate, scaleKin = 1) #convert pcrelate output to GRM and don't scale kinship 

###Make the kinship matrix
#remove duplicate individuals and remove the breed from the sample id
#only keep individuals that have a phenotype
newRowNames = gsub('(.*)-\\w+', '\\1',rownames(pcRelateMat))
colnames(pcRelateMat) <- rownames(pcRelateMat) <- newRowNames
saveRDS(pcRelateMat, "~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/pcRelateMatrix_allIndivs.rds")

snpgdsClose(genofile)
