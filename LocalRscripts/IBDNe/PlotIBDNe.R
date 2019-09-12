#setwd("~/Documents/DogProject_Jaz/IBDNe/Beagle/IBDNe_Beagle")
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(dplyr)
library(randomcoloR)
library(cowplot)

#read files in 
border_collie = read.delim("border_collie_IBDNe_usingBeagleIBDSegs_May18ver.ne")
boxer = read.delim("boxer_IBDNe_usingBeagleIBDSegs_May18ver.ne")
cairn_terrier = read.delim("cairn_terrier_IBDNe_usingBeagleIBDSegs_May18ver.ne")
cocker_spaniel = read.delim("cocker_spaniel_IBDNe_usingBeagleIBDSegs_May18ver.ne")
english_springer_spaniel  = read.delim("english_springer_spaniel_IBDNe_usingBeagleIBDSegs_May18ver.ne")
german_shepherd_dog  = read.delim("german_shepherd_dog_IBDNe_usingBeagleIBDSegs_May18ver.ne")
golden_retriever = read.delim("golden_retriever_IBDNe_usingBeagleIBDSegs_May18ver.ne")
irish_wolfhound = read.delim("irish_wolfhound_IBDNe_usingBeagleIBDSegs_May18ver.ne")
labrador_retriever = read.delim("labrador_retriever_IBDNe_usingBeagleIBDSegs_May18ver.ne")
maltese = read.delim("maltese_IBDNe_usingBeagleIBDSegs_May18ver.ne")
miniature_schnauzer  = read.delim("miniature_schnauzer_IBDNe_usingBeagleIBDSegs_May18ver.ne")
newfoundland  = read.delim("newfoundland_IBDNe_usingBeagleIBDSegs_May18ver.ne")
rottweiler  = read.delim("rottweiler_IBDNe_usingBeagleIBDSegs_May18ver.ne")
village_dog_peru  = read.delim("village_dog_peru_IBDNe_usingBeagleIBDSegs_May18ver.ne")
yorkshire_terrier = read.delim("yorkshire_terrier_IBDNe_usingBeagleIBDSegs_May18ver.ne")
wolves_italy = read.delim("Italy_IBDNE_usingBeagleIBDSegs_May18ver.ne")
wolves_northcentral_europe = read.delim("Northcentral-Europe_IBDNE_usingBeagleIBDSegs_May18ver.ne")

#Loop through and add columns we need
pops = ls()
pops
dfList = mget(pops)
for (i in seq_along(dfList)){
  dfList[[i]]$Population = pops[i]
  dfList[[i]]$Years = dfList[[i]]$GEN*4
}

#Merge to one data frame
allPopsDF = bind_rows(dfList)

#Plot
colourCount_breed = length(pops)
palette = distinctColorPalette(colourCount_breed)

#All on one
allTogether = ggplot(allPopsDF, aes(x=Years, y=NE, colour=Population)) + geom_line(size=1) + scale_colour_manual(values = c(border_collie="#77E650",boxer="#D6E1A2",cairn_terrier="#D3986D",cocker_spaniel ="#D1EB48",english_springer_spaniel="#E350D4",german_shepherd_dog="#DDC552", golden_retriever="#7EE7C2",irish_wolfhound = "#8DADD9" , labrador_retriever ="#749583", maltese = "#DB5265", miniature_schnauzer = "#D7DDD5", newfoundland = "#867BCF", rottweiler = "#867BCF", village_dog_peru = "#8447E4", wolves_italy = "#8CDC83", wolves_northcentral_europe = "#7BD9E2", yorkshire_terrier = "#DAAAC7")) + theme_bw() + scale_y_log10() + theme(axis.text.x = element_text( hjust= 0.5, vjust=1,size=20), axis.text.y = element_text(size =20), plot.title=element_text(size =24, face = "bold", hjust=0.5), axis.title=element_text(size=24)) + theme(legend.title=element_blank(), legend.text=element_text(size=14), legend.position = "bottom") + xlim(0,500)

#Split onto individual
separate = ggplot(allPopsDF, aes(x=Years, y=NE, colour=Population)) + geom_line(size=1) + geom_ribbon(aes(ymin=LWR.95.CI, ymax=UPR.95.CI), alpha=0.1) + scale_colour_manual(values = c(border_collie="#77E650",boxer="#D6E1A2",cairn_terrier="#D3986D",cocker_spaniel ="#D1EB48",english_springer_spaniel="#E350D4",german_shepherd_dog="#DDC552", golden_retriever="#7EE7C2",irish_wolfhound = "#8DADD9" , labrador_retriever ="#749583", maltese = "#DB5265", miniature_schnauzer = "#D7DDD5", newfoundland = "#867BCF", rottweiler = "#867BCF", village_dog_peru = "#8447E4", wolves_italy = "#8CDC83", wolves_northcentral_europe = "#7BD9E2", yorkshire_terrier = "#DAAAC7")) + theme_bw() + scale_y_log10() + theme(axis.text.x = element_text( hjust= 0.5, vjust=1,size=20), axis.text.y = element_text(size =20), plot.title=element_text(size =24, face = "bold", hjust=0.5), axis.title=element_text(size=24)) + theme(legend.title=element_blank(), legend.text=element_text(size=14), legend.position = "bottom") + facet_wrap(~ Population, ncol = 5) + xlim(0,500)

#Plot 
print(allTogether)
print(separate)
plot_grid(allTogether, separate)
