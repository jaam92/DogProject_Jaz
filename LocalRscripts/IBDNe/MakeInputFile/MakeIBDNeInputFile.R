#Load Files and libraries
library(dplyr)
library(data.table)

ibdSegs = fread("~/Documents/DogProject_Jaz/IBDNe/MakeInputFile/cornell_canine_allChroms_phasedHaplotypes_shapeIT.ibd")
Unrelated_sampsGrEql30 = read.delim("~/Documents/DogProject_Jaz/IBDNe/MakeInputFile/UnrelatedIndividuals_grEql30.txt")


#Make data frame with only unrelateds
Unrelated_IBDNe_step1 = ibdSegs %>% filter(ibdSegs$V1 %in% Unrelated_sampsGrEql30$Unrelateds) 
Unrelated_IBDNe = Unrelated_IBDNe_step1 %>% filter(Unrelated_IBDNe_step1$V3 %in% Unrelated_sampsGrEql30$Unrelateds)

#Add column with breed
Unrelated_IBDNe$breed = Unrelated_sampsGrEql30$breed[match(Unrelated_IBDNe$V1, Unrelated_sampsGrEql30$Unrelateds)]
Unrelated_IBDNe$breed2 = Unrelated_sampsGrEql30$breed[match(Unrelated_IBDNe$V3, Unrelated_sampsGrEql30$Unrelateds)]
Unrelated_IBDNe$TestBreed = ifelse(Unrelated_IBDNe$breed2 == Unrelated_IBDNe$breed, 1, 0)
FinalDF = Unrelated_IBDNe[which(Unrelated_IBDNe$TestBreed != 0),]
FinalDF$breed2 = NULL
FinalDF$TestBreed = NULL

#Split on breed and use a factor so that the empty levels are dropped
splitDF = split(FinalDF , f = factor(FinalDF$breed))

for(i in seq_along(splitDF)){
  breedID = as.character(unique(splitDF[[i]][["breed"]]))
  #filename = paste(breedID, "_BeagleIBDSegs_allchroms_unrelatedsOnly.ibd",sep = "")
  filename = paste(breedID, "_IBDSeqIBDSegs_allchroms_unrelatedsOnly.ibd",sep = "")
  splitDF[[i]][["breed"]] = NULL
  write.table(splitDF[[i]], filename, sep = "\t", row.names = F, col.names = F, quote = F)
 
}
