#Download packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("gdsfmt")
#biocLite("SNPRelate")

#Load libraries
library(GWASTools)
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(cowplot)
library(dplyr)

#Load Sample information
setwd("~/Documents/DogProject_Jaz/PCA_Unrelateds")
popmap = read.delim("~/Documents/DogProject_Jaz/BreedCladeInfo/SamplesUsedInfoWolves.txt", fill = T)
sample.id = as.character(popmap$IID_Col2)
population = as.character(popmap$Population)
cluster = as.character(popmap$Cluster)

# PLINK BED files for unpruned data aka random sample of 100000 snps
setwd("~/Documents/DogProject_Jaz/PCA_Unrelateds/")
#bed.fn = ("WolvesAll_22feb2012.bed")
#bim.fn = ("WolvesAll_22feb2012.bim")
#fam.fn = ("WolvesAll_22feb2012.fam")

# convert
#snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "WolvesAll_22feb2012.gds")

#look at file summary 
#snpgdsSummary("WolvesAll_22feb2012.gds")

#Open file
genofile = snpgdsOpen("WolvesAll_22feb2012.gds")

#LD prune data
snpset = snpgdsLDpruning(genofile, ld.threshold=0.8, maf = 0.05, missing.rate = 0.1, slide.max.n = 50, autosome.only = F)
snpset.id = unlist(snpset)

#Run PCA with SNPs at frequecy higher than 5% and missingness less than 10%
pca = snpgdsPCA(genofile, snp.id=snpset.id, autosome.only = F)

#Look at percentage of variance explained by each PC
pc.percent = pca$varprop*100
pc = head(round(pc.percent,2))
pc


#Make data frame with first two pcs
df_PCA = data.frame(sample.id = pca$sample.id, population = factor(population)[match(pca$sample.id, sample.id)], cluster = factor(cluster)[match(pca$sample.id, sample.id)], EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], EV3 = pca$eigenvect[,3], EV4 = pca$eigenvect[,4],stringsAsFactors = FALSE)

head(df_PCA)
tail(df_PCA)

#plot the nice version
WolvesPC1vPC2 = ggplot(df_PCA, aes(y=EV2, x=EV1, colour=population)) + geom_point(size=2) + theme_bw() + labs(y=bquote('PC2' ~'('~.(pc[2])~'%'~')'), x=bquote('PC1'~'('~.(pc[1])~'%'~')'))  + theme(axis.text.x = element_text(size  = 24), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=24), legend.text=element_text(size=24))

WolvesPC2vPC3 = ggplot(df_PCA, aes(y=EV3, x=EV2, colour=population)) + geom_point(size=2) + theme_bw() + labs(y=bquote('PC3' ~'('~.(pc[3])~'%'~')'), x=bquote('PC2'~'('~.(pc[2])~'%'~')'))  + theme(axis.text.x = element_text(size  = 24), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=24), legend.text=element_text(size=24))

print(WolvesPC1vPC2)
print(WolvesPC2vPC3)

###PCA on non-italian populations###
sampleid_nonItalian = popmap %>% filter(Population !="Italy") %>% select(IID_Col2)
samp.id.nonItalian = unlist(sampleid_nonItalian)

#LD prune
snpset_nonItalian = snpgdsLDpruning(genofile, sample.id = samp.id.nonItalian, ld.threshold=0.8, maf = 0.05, missing.rate = 0.1, slide.max.n = 50, autosome.only = F)
snpset.id.nonItalian = unlist(snpset_nonItalian)

#Run PCA
pca_nonItalian = snpgdsPCA(genofile, snp.id = snpset.id.nonItalian, sample.id = samp.id.nonItalian, autosome.only = F)

#Look at percentage of variance explained by each PC
pc.percent_nonItalian = pca_nonItalian$varprop*100
pc_nonItalian = head(round(pc.percent_nonItalian,2))
pc_nonItalian


#Make data frame with first two pcs
df_PCA_nonItalian = data.frame(sample.id = pca_nonItalian$sample.id, population = factor(population)[match(pca_nonItalian$sample.id, sample.id)], cluster = factor(cluster)[match(pca_nonItalian$sample.id, sample.id)], EV1 = pca_nonItalian$eigenvect[,1], EV2 = pca_nonItalian$eigenvect[,2], EV3 = pca_nonItalian$eigenvect[,3], EV4 = pca_nonItalian$eigenvect[,4],stringsAsFactors = FALSE)

head(df_PCA_nonItalian)
tail(df_PCA_nonItalian)

#plot the nice version
WolvesPC1vPC2_nonItalian = ggplot(df_PCA_nonItalian, aes(y=EV2, x=EV1, colour=cluster)) + geom_point(size=2) + theme_bw() + labs(y=bquote('PC2' ~'('~.(pc_nonItalian[2])~'%'~')'), x=bquote('PC1'~'('~.(pc_nonItalian[1])~'%'~')'))  + theme(axis.text.x = element_text(size  = 24), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=24), legend.text=element_text(size=24))

WolvesPC2vPC3_nonItalian = ggplot(df_PCA_nonItalian, aes(y=EV3, x=EV2, colour=population)) + geom_point(size=2) + theme_bw() + labs(y=bquote('PC3' ~'('~.(pc_nonItalian[3])~'%'~')'), x=bquote('PC2'~'('~.(pc_nonItalian[2])~'%'~')'))  + theme(axis.text.x = element_text(size  = 24), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=24), legend.text=element_text(size=24))

print(WolvesPC1vPC2_nonItalian)
print(WolvesPC2vPC3_nonItalian)


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
unrelated.ibd.robust = snpgdsIBDKING(genofile, sample.id=Unrelateds,family.id=popmap$FID_Col1[match(Unrelateds, popmap$IID_Col2)], snp.id=snpset.id, num.thread=2, autosome.only = F)
unrelatedKING = snpgdsIBDSelection(unrelated.ibd.robust)

#check that individuals areat most first cousins and write unrelateds to file
summary(unrelatedKING$kinship)
dfUnrelateds_Full = cbind.data.frame(Unrelateds, popmap$FID_Col1[match(Unrelateds, popmap$IID_Col2)],popmap$Population[match(Unrelateds, popmap$IID_Col2)], popmap$Cluster[match(Unrelateds, popmap$IID_Col2)])
names(dfUnrelateds_Full)[2] = "FID"
names(dfUnrelateds_Full)[3] = "Population"
names(dfUnrelateds_Full)[4] = "Cluster"
#write.table(dfUnrelateds_Full, "WolvesUnrelatedIndividuals_allBreeds.txt", sep = "\t", row.names = F, col.names = F, quote = F)

#Make a larger data frame with population and cluster info
UnrelatedsPerBreed_n30 = dfUnrelateds_Full %>% group_by(Cluster) %>% tally() %>% filter(n>=30) #find breeds with at least 30 unrelateds
Unrelated_sampsGrEql30 = dfUnrelateds_Full %>% filter(Cluster %in% UnrelatedsPerBreed_n30$Cluster) %>% select(FID,Unrelateds)#output in order of plink files
Unrelated_sampsGrEql30$breed = popmap$Cluster[match(Unrelated_sampsGrEql30$Unrelateds, popmap$IID_Col2)]#Save it as breed because this is colname for Rscript that parses ibd segments expects
#write.table(Unrelated_sampsGrEql30, "WolvesUnrelatedIndividuals_grEql30.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#Close file
snpgdsClose(genofile)

#Plot Unrelateds
KINGdf$Cluster1 = popmap$Cluster[match(KINGdf$ID1, popmap$IID_Col2)] 
KINGdf$Cluster2 = popmap$Cluster[match(KINGdf$ID2, popmap$IID_Col2)] 
KINGdf$SameCluster = ifelse(KINGdf$Cluster1 == KINGdf$Cluster2, 1, 0)

#Check Italy against other Clusters
ItalyVsEveryone = ggplot(subset(KINGdf, KINGdf$Cluster1 == "Italy"), aes(x=IBS0,y=kinship, colour=Cluster2)) + geom_point() + theme_bw() 

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
ItalyVsEveryBoxes = ggplot()+
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
  geom_point(data=subset(KINGdf, KINGdf$Cluster1 == "Italy"), aes(x=IBS0,y=kinship, colour=Cluster2))+
  xlab("Proportion of Zero IBS")+
  ylab("Estimated Kinship Coefficient (KING-robust)")+
  theme_bw()+
  scale_x_continuous(limits = c(0,1))+
  theme(legend.title = element_blank()) 


plot_grid(ItalyVsEveryone, ItalyVsEveryBoxes, labels = c("A", "B"), ncol = 1, align = 'h')
