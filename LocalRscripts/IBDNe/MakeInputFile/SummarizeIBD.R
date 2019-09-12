#Load Files and libraries
library(dplyr)
library(data.table)

#####Load Files
omiaGenes = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/OMIA/omia-genes_v2.txt", fill = TRUE)
causalVars = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/OMIA/causalVars_OMIA.txt")
ibdSegs = fread("~/Documents/DogProject_Jaz/LocalRscripts/IBDNe/MakeInputFile/cornell_canine_allChroms_phasedHaplotypes_shapeIT.ibd")
Unrelated_sampsGrEql30 = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/IBDNe/MakeInputFile/UnrelatedIndividuals_grEql30.txt")
popmap = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/breeds_dryad.txt")
cladeInfo = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/BreedAndCladeInfo.txt")
censusSizes = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/CensusData_sorted_May2014.txt")
#####Modify OMIA Files

#Reformat OMIA files
causalVars$NCBI_gene_ID = omiaGenes$ncbi_gene_id[match(causalVars$Gene, omiaGenes$gene_symbol)] #add ncbi gene id
sepBreedCausalVars = causalVars %>% mutate(Breed.s.= strsplit(as.character(Breed.s.), ",")) %>% unnest(Breed.s.) #split the breeds that are comma delimited into separate rows with same info
names(sepBreedCausalVars)[19] = "Breed"
sepBreedCausalVars$Breed = tolower(sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub(" ", "_", sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("-", "_", sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("^\\_","",sepBreedCausalVars$Breed)
#Rename breeds
sepBreedCausalVars$Breed = gsub("chinese_crested_dog", "chinese_crested",sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("chinese_shar_pei", "chinese_shar-pei",sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("brittany_spaniel", "brittany", sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("french_bull_dog", "bulldog_french", sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("bull_mastiff", "bullmastiff", sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("curly_coated_retriever", "curly-coated_retriever", sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("greater_swiss_mountain", "greater_swiss_mountain_dog", sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("japanese_chin", "japanese_chin_dog", sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("standard_poodle", "poodle", sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("italian_spinone", "spinone_italiano", sepBreedCausalVars$Breed)
sepBreedCausalVars$Breed = gsub("wirehaired_fox_terrier", "wire_fox_terrier", sepBreedCausalVars$Breed)

#Count up Breeds with causal variants
BreedsWithCausalVars = sepBreedCausalVars %>% count(Breed) %>% as.data.frame()#count causal vars associated with each breed
BreedsWithCausalVars$Clade = cladeInfo$clade[match(BreedsWithCausalVars$Breed, cladeInfo$breed)]

####Modify IBD Files
#Add new cols
ibdSegs$totalLen = ibdSegs$V7 -ibdSegs$V6
ibdSegs$Breed = popmap$breed[match(ibdSegs$V1, popmap$dogID)]

#Compute IBD Segments per breed
MeanIBDperBreed = MeanIBDperBreed = ibdSegs %>% group_by(V1,Breed) %>% summarize(LenAuto = sum(as.numeric(totalLen))) %>% group_by(Breed) %>% summarize(MeanIBDperBreed = mean(LenAuto)) %>% as.data.frame()
CountIBDperBreed = ibdSegs %>% count(Breed) %>% as.data.frame()

#Add causal vars columns
MeanIBDperBreed$CountVars = BreedsWithCausalVars$n[match(MeanIBDperBreed$Breed, BreedsWithCausalVars$Breed)]
MeanIBDperBreed$CensusSize = censusSizes$CensusSize[match(MeanIBDperBreed$Breed, censusSizes$Breed)]

#Check correlation
corrIBDVars = lm(CountVars~MeanIBDperBreed, data = MeanIBDperBreed)
corrIBDNc = lm(CensusSize~MeanIBDperBreed, data = MeanIBDperBreed)

###LinearRegression Function###
ggplotRegression = function (fit) {
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + geom_point(size = 2) + stat_smooth( method = 'lm', col = "blue") +  theme_bw() + 
    labs(title = bquote(R^2== ~.(signif(summary(fit)$adj.r.squared, 5))~"&"~"p"==~.(signif(summary(fit)$coef[2,4], 5))))
  
}
#Plot
ggplotRegression(corrIBDVars) + labs(x="Average Length of Genome in IBD Segment per Breed", y="Count Causal Variants") + theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=24), legend.text=element_text(size=18), legend.position = "none") + guides(fill = guide_legend(nrow = 4))

ggplotRegression(corrIBDNc) + labs(x="Average Length of Genome in IBD Segment per Breed", y = expression(N[C])) + theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=24), legend.text=element_text(size=18), legend.position = "none") + guides(fill = guide_legend(nrow = 4))


############Redo Correlations without lab
rmLab = subset(MeanIBDperBreed, Breed != "labrador_retriever")
#Check correlation
corrIBDVars = lm(CountVars~MeanIBDperBreed, data = rmLab)
corrIBDNc = lm(CensusSize~MeanIBDperBreed, data = rmLab)

#Plot
ggplotRegression(corrIBDVars) + labs(x="Average Length of Genome in IBD Segment per Breed (No Lab)", y="Count Causal Variants") + theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=24), legend.text=element_text(size=18), legend.position = "none") + guides(fill = guide_legend(nrow = 4))

ggplotRegression(corrIBDNc) + labs(x="Average Length of Genome in IBD Segment per Breed (No Lab)", y = expression(N[C])) + theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=24), legend.text=element_text(size=18), legend.position = "none") + guides(fill = guide_legend(nrow = 4))

