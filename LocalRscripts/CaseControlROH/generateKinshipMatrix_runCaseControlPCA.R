#Load libraries
library(tidyverse)
library(gdsfmt)
library(SNPRelate)
library(GENESIS)
library(GWASTools)

####Make gds file
#bed.fn = ("~/DogProject_Jaz/LocalRscripts/CaseControlROH/MergedFile_CornellCanineFitak_UnrelatedsOnly.bed")
#bim.fn = ("~/DogProject_Jaz/LocalRscripts/CaseControlROH/MergedFile_CornellCanineFitak_UnrelatedsOnly.bim")
#fam.fn = ("~/DogProject_Jaz/LocalRscripts/CaseControlROH/MergedFile_CornellCanineFitak_UnrelatedsOnly.fam")

#convert
#snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "~/DogProject_Jaz/LocalRscripts/CaseControlROH/MergedFile_CornellCanineFitak_UnrelatedsOnly.gds")
#snpgdsSummary("~/DogProject_Jaz/LocalRscripts/CaseControlROH/MergedFile_CornellCanineFitak_UnrelatedsOnly.gds")

#open the new gds file
genofile = snpgdsOpen("~/DogProject_Jaz/LocalRscripts/CaseControlROH/MergedFile_CornellCanineFitak_UnrelatedsOnly.gds")
sampIds = read.gdsn(index.gdsn(genofile, "sample.id")) #grab sample ids 
famIds = gsub(".*-","",sampIds) #make family ids

#LD prune because data set is small, set r2 = 0.7, maf 5%, and 50 snp window 
snpset = snpgdsLDpruning(genofile, sample.id=sampIds, method="corr", slide.max.n=50, ld.threshold=0.3, maf = 0.05, autosome.only = F) 
snpset.id = unlist(snpset)

#IBD and kinship analysis with KING
#Must be done on LD pruned data

#Get sample ids and breed info from gds file
gdsSampIDs = read.gdsn(index.gdsn(genofile, "sample.id"))
breed.id = gsub(".*-","",gdsSampIDs)#remove everything before the dash to get breed

#Run PCA on goldens
sampleid_goldens = gdsSampIDs[grep("golden_retriever", gdsSampIDs)]
samp.id.goldens = unlist(sampleid_goldens)

#LD prune
snpset_goldens = snpgdsLDpruning(genofile, sample.id = samp.id.goldens, method="corr", slide.max.n=50, ld.threshold=0.3, maf = 0.05, autosome.only = F)
snpset.id.goldens = unlist(snpset_goldens)

#Run PCA
pca_goldens = snpgdsPCA(genofile, snp.id = snpset.id.goldens, sample.id = samp.id.goldens, autosome.only = F)

#Look at percentage of variance explained by each PC
pc.percent_goldens = pca_goldens$varprop*100
pc_goldens = head(round(pc.percent_goldens,2))
pc_goldens
pca_goldens$sample.id = gsub('(.*)-\\w+', '\\1', pca_goldens$sample.id)

#Make data frame with first two pcs
df_PCA_goldens = data.frame(sample.id = pca_goldens$sample.id, 
                            status = factor(as.character(phenotypes$lymphoma_goldenRetrievers))[match(pca_goldens$sample.id, phenotypes$dogID)], 
                            EV1 = pca_goldens$eigenvect[,1], 
                            EV2 = pca_goldens$eigenvect[,2], 
                            EV3 = pca_goldens$eigenvect[,3], 
                            EV4 = pca_goldens$eigenvect[,4], 
                            stringsAsFactors = FALSE) %>%
  na.omit()

#Plot the pca
ggplot(df_PCA_goldens, aes(y=EV2, x=EV1, colour=status)) +
  geom_point(size=2) +
  labs(y=bquote('PC2' ~'('~.(pc_goldens[2])~'%'~')'), x=bquote('PC1'~'('~.(pc_goldens[1])~'%'~')')) +
  theme_bw() +
  theme(axis.text.x = element_text(size  = 24), 
        axis.text.y = element_text(size  = 24), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=24), 
        legend.text=element_text(size=24))


###Make a kinship matrix that accounts for known relatedness of samples
#Run king
king = snpgdsIBDKING(genofile, snp.id=snpset.id)
kingMat = king$kinship
colnames(kingMat) <- rownames(kingMat) <- king$sample.id
snpgdsClose(genofile)#Close file

#Re-Open GDS data 
genoFile = GdsGenotypeReader(filename = "~/DogProject_Jaz/LocalRscripts/CaseControlROH/MergedFile_CornellCanineFitak_UnrelatedsOnly.gds")#read in GDS data
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
#saveRDS(pcRelateMat, "~/DogProject_Jaz/LocalRscripts/CaseControlROH/pcRelateMatrix_Unrelateds.rds")

snpgdsClose(genofile)
