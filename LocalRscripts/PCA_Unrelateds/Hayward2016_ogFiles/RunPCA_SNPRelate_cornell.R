#Download packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("gdsfmt")
#biocLite("SNPRelate")

#Load libraries
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(randomcoloR)
library(dplyr)

# PLINK BED files for unpruned data aka random sample of 100000 snps
setwd("~/Documents/DogProject_Jaz/PCA_Unrelateds/Hayward2016_ogFiles/")
bed.fn = ("cornell_canine_updatedFID.bed")
bim.fn = ("cornell_canine_updatedFID.bim")
fam.fn = ("cornell_canine_updatedFID.fam")

#convert
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "cornell_canine_updatedFID.gds")
snpgdsSummary("cornell_canine_updatedFID.gds")

#Open file
genofile = snpgdsOpen("cornell_canine_updatedFID.gds")

#LD prune
snpset = snpgdsLDpruning(genofile, ld.threshold=0.8, maf = 0.05, missing.rate = 0.1, slide.max.n = 50, autosome.only = F)
snpset.id = unlist(snpset)

#Run PCA with SNPs at frequecy higher than 5% and missingness less than 10%
pca = snpgdsPCA(genofile, snp.id=snpset.id, autosome.only = F)

#Look at percentage of variance explained by each PC
pc.percent = pca$varprop*100
pc = head(round(pc.percent,2))
pca$sample.id = gsub("-.*","",pca$sample.id) #remove the breed info from sampleID
pc

#Load Sample information
popmap = read.delim("~/Documents/DogProject_Jaz/BreedCladeInfo/breeds_dryad.txt")
popmap$breed = gsub("_dog_.*","",popmap$breed) #turn village dog to just village
sample.id = as.character(popmap$dogID)
population = as.character(popmap$breed)
breedCladeInfo = read.delim("~/Documents/DogProject_Jaz/BreedCladeInfo/BreedAndCladeInfo_mergedFitakCornell.txt")

#Make data frame with first two pcs
df_PCA = data.frame(sample.id = pca$sample.id, population = factor(population)[match(pca$sample.id, sample.id)], EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], EV3 = pca$eigenvect[,3], EV4 = pca$eigenvect[,4],stringsAsFactors = FALSE)

df_PCA$clade = breedCladeInfo$clade[match(df_PCA$population, breedCladeInfo$breed)]
#head(df_PCA)
#tail(df_PCA)

#Evaluate data frame with breed and clade information
#expand color palette of choice to hve number of colors equal to number of clades
colourCount = length(unique(df_PCA$clade))
palette = distinctColorPalette(colourCount)

#plot the nice version
a = ggplot(df_PCA, aes(y=EV2, x=EV1, colour=clade)) + geom_point(size=2) + scale_colour_manual(values = palette) + theme_bw() + labs(y=bquote('PC2' ~'('~.(pc[2])~'%'~')'), x=bquote('PC1'~'('~.(pc[1])~'%'~')'))  + theme(axis.text.x = element_text(size  = 24), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=24), legend.text=element_text(size=24))

b = ggplot(df_PCA, aes(y=EV3, x=EV2, colour=clade)) + geom_point(size=2) + scale_colour_manual(values = palette) + theme_bw() + labs(y=bquote('PC3' ~'('~.(pc[3])~'%'~')'), x=bquote('PC2'~'('~.(pc[2])~'%'~')'))  + theme(axis.text.x = element_text(size  = 24), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=24), legend.text=element_text(size=24))

print(a)
print(b)


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
write.table(Unrelated_sampsGrEql30, "UnrelatedIndividuals_grEql30.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#Close file
snpgdsClose(genofile)

#Plot Unrelateds
KINGdf$BreedID1 = gsub(".*-", "", KINGdf$ID1, perl = T)
KINGdf$BreedID2 = gsub(".*-", "", KINGdf$ID2, perl = T)
KINGdf$SamePop = ifelse(KINGdf$BreedID1 == KINGdf$BreedID2, 1, 0)

#Check maltese against other breeds
checkAcrossBreeds = subset(KINGdf, KINGdf$ID1 == "PFZ44E12-labrador_retriever")
checkAcrossBreeds$Breed = gsub(".*-", "", checkAcrossBreeds$ID2, perl = T)
ggplot(checkAcrossBreeds, aes(x=IBS0,y=kinship, colour=Breed)) + geom_point() + theme_bw() + theme(legend.position="none") + ylim(-0.25,0)

########################## Manichaikul boxes ################
# from the Manichaikul manuscript: cutoffs (see my Manichaikul presentation)
# values: 
a = 1/ (2^(3/2)) # 0.3535534
b = 1/ (2^(5/2)) # 0.1767767
c = 1/ (2^(7/2)) # 0.08838835 
d = 1/ (2^(9/2)) # 0.04419417
unrelatedBox=data.frame(ymin=0,ymax=d,xmin=(1-b),xmax=1)
thirdDegreeBox=data.frame(ymin=d,ymax=c,xmin=1-a,xmax=1-b)
secondDegreeBox=data.frame(ymin=c,ymax=b,xmin=a,xmax=0.365)
fullSibBox=data.frame(ymin=b,ymax=a,xmin=0.1,xmax=0.365)
ParentBox=data.frame(ymin=b,ymax=a,xmin=0,xmax=0.1)
MonozygoticBox=data.frame(ymin=a,ymax=Inf,xmin=0,xmax=0.1)
# these values are from manichaikul 2010 table 1 
# when using geom rect set aes for each thing separately

#Check one breed 
subsetDF = subset.data.frame(KINGdf, BreedID1 == "maltese" & BreedID2 == "maltese")
MalteseVMaltese = ggplot()+
  geom_rect(data=unrelatedBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3)+
  geom_text(aes(x=unrelatedBox$xmin+0.1,y=unrelatedBox$ymin+0.02,label="unrelated"))+
  geom_rect(data=thirdDegreeBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3)+
  geom_text(aes(x=thirdDegreeBox$xmin+0.1,y=thirdDegreeBox$ymin+0.02,label="third degree"))+
  geom_rect(data=secondDegreeBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3)+
  geom_text(aes(x=secondDegreeBox$xmin+0.05,y=secondDegreeBox$ymin+0.02,label="second degree"))+
  geom_rect(data=fullSibBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3)+
  geom_text(aes(x=fullSibBox$xmin+0.1,y=fullSibBox$ymin+0.1,label="full sib"))+
  geom_rect(data=ParentBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3,fill="blue")+
  geom_text(aes(x=ParentBox$xmin+.1,y=ParentBox$ymin+.15,label="parent-offspring"))+
  geom_rect(data=MonozygoticBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3,fill="red")+
  geom_text(aes(x=MonozygoticBox$xmin + 0.05,y=MonozygoticBox$ymin + 0.15,label="monozygotic"))+
  geom_point(data=subsetDF,aes(x=IBS0, y=kinship,color=BreedID1))+
  xlab("Proportion of Zero IBS")+
  ylab("Estimated Kinship Coefficient (KING-robust)")+
  theme_bw()+
  scale_x_continuous(limits = c(0,1))+
  theme(legend.title = element_blank()) +
  ylim(-0.5,0.5)

#Check closely related breeds
subsetDF = subset.data.frame(KINGdf, BreedID1 == "maltese" | BreedID1 == "coton_de_tulear" | BreedID1 == "havanese")
MalteseVCloseBreeds = ggplot()+
  geom_rect(data=unrelatedBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3)+
  geom_text(aes(x=unrelatedBox$xmin+0.1,y=unrelatedBox$ymin+0.02,label="unrelated"))+
  geom_rect(data=thirdDegreeBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3)+
  geom_text(aes(x=thirdDegreeBox$xmin+0.1,y=thirdDegreeBox$ymin+0.02,label="third degree"))+
  geom_rect(data=secondDegreeBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3)+
  geom_text(aes(x=secondDegreeBox$xmin+0.05,y=secondDegreeBox$ymin+0.02,label="second degree"))+
  geom_rect(data=fullSibBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3)+
  geom_text(aes(x=fullSibBox$xmin+0.1,y=fullSibBox$ymin+0.1,label="full sib"))+
  geom_rect(data=ParentBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3,fill="blue")+
  geom_text(aes(x=ParentBox$xmin+.1,y=ParentBox$ymin+.15,label="parent-offspring"))+
  geom_rect(data=MonozygoticBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3,fill="red")+
  geom_text(aes(x=MonozygoticBox$xmin + 0.05,y=MonozygoticBox$ymin + 0.15,label="monozygotic"))+
  geom_point(data=subsetDF,aes(x=IBS0, y=kinship,color=BreedID1))+
  xlab("Proportion of Zero IBS")+
  ylab("Estimated Kinship Coefficient (KING-robust)")+
  theme_bw()+
  scale_x_continuous(limits = c(0,1))+
  theme(legend.title = element_blank()) +
  ylim(-0.5,0.5)

#Check one breed versus a bunch 
subsetDF = subset.data.frame(KINGdf, BreedID1 == "maltese")
MalteseVMany = ggplot()+
  geom_rect(data=unrelatedBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3)+
  geom_text(aes(x=unrelatedBox$xmin+0.1,y=unrelatedBox$ymin+0.02,label="unrelated"))+
  geom_rect(data=thirdDegreeBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3)+
  geom_text(aes(x=thirdDegreeBox$xmin+0.1,y=thirdDegreeBox$ymin+0.02,label="third degree"))+
  geom_rect(data=secondDegreeBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3)+
  geom_text(aes(x=secondDegreeBox$xmin+0.05,y=secondDegreeBox$ymin+0.02,label="second degree"))+
  geom_rect(data=fullSibBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3)+
  geom_text(aes(x=fullSibBox$xmin+0.1,y=fullSibBox$ymin+0.1,label="full sib"))+
  geom_rect(data=ParentBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3,fill="blue")+
  geom_text(aes(x=ParentBox$xmin+.1,y=ParentBox$ymin+.15,label="parent-offspring"))+
  geom_rect(data=MonozygoticBox,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.3,fill="red")+
  geom_text(aes(x=MonozygoticBox$xmin + 0.05,y=MonozygoticBox$ymin + 0.15,label="monozygotic"))+
  geom_point(data=subsetDF,aes(x=IBS0, y=kinship,color=BreedID2))+
  xlab("Proportion of Zero IBS")+
  ylab("Estimated Kinship Coefficient (KING-robust)")+
  theme_bw()+
  scale_x_continuous(limits = c(0,1))+
  theme(legend.title = element_blank(), legend.position = "none") +
  ylim(-0.5,0.5)

plot_grid(MalteseVMaltese, MalteseVCloseBreeds, MalteseVMany, labels = c("A", "B", "C"), ncol = 1, align = 'h')