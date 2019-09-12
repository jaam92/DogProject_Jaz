####Clade Information###
#Clades determined by elaine ostranders group and adam boyko's group
#paper titles: Genomic Analyses Reveal the Influence of Geographic Origin, Migration, and Hybridization on Modern Dog Breed Development and Genetic structure in village dogs reveals a Central Asian domestication origin

#Load Libraries
library(dplyr)
library(data.table)

#Read files in
setwd("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/")
popmap = read.delim("MergedBreeds_CornellFitak.txt")
cladeInfoWolves = read.table("popInfoWolves.txt")
cladeInfoDogs = read.delim("cladeInfo.txt")

#Make all village dogs a single breed
wolves = popmap[grep("Wolf", popmap$breed),]
dogs = popmap[!grepl("Wolf", popmap$breed),]
dogs$breed = gsub("_dog_.*","",dogs$breed) #turn village dog to just village

#Add on clades 
dogs$clade = cladeInfoDogs$Clade[match(dogs$breed, cladeInfoDogs$Breed)]

#set cocker spaniel to spaniel clade
setDT(dogs)[breed=="cocker_spaniel", clade:="Spaniel"]
#set english_toy_spaniel to Spaniel
setDT(dogs)[breed=="english_toy_spaniel", clade:="Spaniel"]

#set chesapeake_bay_retriever to retriever
setDT(dogs)[breed=="chesapeake_bay_retriever", clade:="Retriever"]

#set mini daschund to Scent_Hound clade
setDT(dogs)[breed=="dachshund_miniature", clade:="Scent_Hound"]

#set munsterlander_small to Pointer_Setter
setDT(dogs)[breed=="munsterlander_small", clade:="Pointer_Setter"]
#set dalmatian to Point_Setter
setDT(dogs)[breed=="dalmatian_croatian", clade:="Pointer_Setter"]

#set village to village
setDT(dogs)[breed=="village", clade:="Village"]

#set mastiff to European_Mastiff
setDT(dogs)[breed=="mastiff", clade:="European_Mastiff"]
#set american_pit_bull_terrier to European Mastiff
setDT(dogs)[breed=="american_pit_bull_terrier", clade:="European_Mastiff"]

#set poodles to poodle
setDT(dogs)[breed=="poodle_miniature", clade:="Poodle"]
setDT(dogs)[breed=="poodle_toy", clade:="Poodle"]

#Set mixes to mix
setDT(dogs)[breed=="mix", clade:="Mix"]

#Set mixes to mix
setDT(dogs)[breed=="large_munsterlander", clade:="Pointer_Setter"]
#Check out NAs
#View(dogs[is.na(dogs$clade),]) #24 individuals are NA

#Work with wolves to add on clades and reformat
wolves$clade = cladeInfoWolves$V3[match(wolves$dogID, cladeInfoWolves$V1)]
newWolves = cbind.data.frame(wolves$dogID,wolves$clade,wolves$breed) #I think the clusters(currently listed as breed) should be the clades 
names(newWolves)[1] = "dogID"
names(newWolves)[2] = "breed"
names(newWolves)[3] = "clade"

#Put everything together
SampleInfo = rbind.data.frame(dogs,newWolves)
#write.table(SampleInfo, file = "BreedAndCladeInfo_mergedFitakCornell.txt", quote = F, col.names = T, row.names = F, sep = "\t")
