#Load libraries
library(tidyverse)
library(gdsfmt)
library(SNPRelate)

#Load files
popmapDryad = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/breeds_dryad.txt")
phenotypes = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/phenotypes.txt")
unrelateds = read.table("~/Documents/DogProject_Jaz/LocalRscripts/PCA_Unrelateds/UnrelatedIndividuals_allBreeds_mergedFitakCornell.txt")
genofile = snpgdsOpen("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/MergedFile_CornellCanineFitak_UnrelatedsOnly.gds")

#Get sample ids and breed info from gds file
gdsSampIDs = read.gdsn(index.gdsn(genofile, "sample.id"))
breed.id = gsub(".*-","",gdsSampIDs)#remove everything before the dash to get breed

#function for PCA
runPCA = function(phenoColName, plotTitle){
  phenoCol = enquo(phenoColName)
  phenoOfInterest = FinalPhenotypes %>% 
    select("dogID", "breed", !!phenoCol) %>% 
    rename("status" = phenoColName) %>%
    na.omit()
  
  #grab samples for a given trait 
  #have to use gsub here because the breed is attache 
  sampsToUse = gdsSampIDs[gsub('(.*)-\\w+', '\\1',gdsSampIDs) %in% phenoOfInterest$dogID]
  samp.ids = unlist(sampsToUse)
  
  #LD prune
  snpset = snpgdsLDpruning(genofile, sample.id = samp.ids, method="corr", slide.max.n=50, ld.threshold=0.3, maf = 0.05, autosome.only = F)
  snpset.ids = unlist(snpset)
  
  #Run PCA
  pca = snpgdsPCA(genofile, snp.id = snpset.ids, sample.id = samp.ids, autosome.only = F)
  #Look at percentage of variance explained by each PC
  pc.percent = pca$varprop*100
  pc = head(round(pc.percent,2))
  pca$sample.id = gsub('(.*)-\\w+', '\\1', pca$sample.id)
  
  #Make data frame with first two pcs
  df_PCA = data.frame(sample.id = pca$sample.id, 
                      status = factor(as.character(phenoOfInterest$status))[match(pca$sample.id, phenoOfInterest$dogID)],
                      breed = factor(as.character(phenoOfInterest$breed))[match(pca$sample.id, phenoOfInterest$dogID)],
                      EV1 = pca$eigenvect[,1], 
                      EV2 = pca$eigenvect[,2], 
                      EV3 = pca$eigenvect[,3], 
                      EV4 = pca$eigenvect[,4], 
                      stringsAsFactors = FALSE) %>%
    na.omit()
  
  #Plot the pca
  pc1vs2 = ggplot(df_PCA, aes(y=EV2, x=EV1, colour=status)) +
    geom_point(size=2) +
    labs(y=bquote('PC2' ~'('~.(pc[2])~'%'~')'), x=bquote('PC1'~'('~.(pc[1])~'%'~')'), title = paste(plotTitle)) +
    theme_bw() + 
    theme(axis.text.x = element_text(size = 20), 
          axis.text.y = element_text(size = 10), 
          plot.title=element_text(size=26, face = "bold", hjust=0.5), 
          axis.title=element_text(size=20),
          legend.title=element_text(size=20), 
          legend.text=element_text(size=18))
  
  dataFrameAndPCA = list("dataFrame" = df_PCA, "plot" = pc1vs2)
  return(dataFrameAndPCA)
}

#Make phenotypes data frame
FinalPhenotypes = phenotypes %>% 
  mutate(PSVA = ifelse(is.na(PSVA), PSVA_yorkshireTerriers, PSVA), 
         MCT = ifelse(is.na(MCT), MCT_labradorRetrievers, MCT), 
         lymphoma = ifelse(is.na(lymphoma), lymphoma_goldenRetrievers, lymphoma),
         Unrelated = ifelse(dogID %in% unrelateds$V1, 1, 0), 
         breed = popmapDryad$breed[match(dogID, popmapDryad$dogID)]) %>% 
  filter(Unrelated == 1) 

#Grab phenotypes of interest
phenotypes = colnames(FinalPhenotypes)[3:13]

pdf(file = "~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/allTraits_PCA.pdf", height = 8, width = 10)

for(i in phenotypes){
  title = paste(i)
  x = runPCA(i, title)
  print(x["plot"])
}

dev.off()
