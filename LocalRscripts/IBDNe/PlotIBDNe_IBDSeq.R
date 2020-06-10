setwd("~/Documents/DogProject_Jaz/LocalRscripts/IBDNe/AllSitesMerge/")
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(dplyr)
library(randomcoloR)
library(cowplot)

#read files in 
boxer = read.delim("boxer_IBDNE_usingIBDSeqIBDSegs_Sept15ver.ne")
cocker_spaniel = read.delim("cocker_spaniel_IBDNe_usingIBDSeqIBDSegs_Sept15ver.ne")
german_shepherd_dog  = read.delim("german_shepherd_dog_IBDNe_usingIBDSeqIBDSegs_Sept15ver.ne")
golden_retriever = read.delim("golden_retriever_IBDNe_usingIBDSeqIBDSegs_Sept15ver.ne")
grayWolf_Europe = read.delim("grayWolf_Europe_IBDNe_usingIBDSeqIBDSegs_Sept15ver.ne")
grayWolf_NorthAmerica = read.delim("grayWolf_NorthAmerica_IBDNe_usingIBDSeqIBDSegs_Sept15ver.ne")
labrador_retriever = read.delim("labrador_retriever_IBDNe_usingIBDSeqIBDSegs_Sept15ver.ne")
maltese = read.delim("maltese_IBDNe_usingIBDSeqIBDSegs_Sept15ver.ne")
mixed  = read.delim("mix_IBDNe_usingIBDSeqIBDSegs_Sept15ver.ne")
newfoundland  = read.delim("newfoundland_IBDNe_usingIBDSeqIBDSegs_Sept15ver.ne")
poodle  = read.delim("poodle_IBDNe_usingIBDSeqIBDSegs_Sept15ver.ne")
rottweiler  = read.delim("rottweiler_IBDNe_usingIBDSeqIBDSegs_Sept15ver.ne")
village_dog_peru  = read.delim("village_dog_peru_IBDNe_usingIBDSeqIBDSegs_Sept15ver.ne")
yorkshire_terrier = read.delim("yorkshire_terrier_IBDNe_usingIBDSeqIBDSegs_Sept15ver.ne")
#wolves_italy = read.delim("Italy_IBDNE_usingIBDSeqSegs_Sept15ver.ne")
#wolves_northcentral_europe = read.delim("Northcentral-Europe_IBDNE_usingIBDSeqSegs_Sept15ver.ne")

#Loop through and add columns we need
pops = ls()
pops
dfList = mget(pops)
for (i in seq_along(dfList)){
  dfList[[i]]$Population = pops[i]
  dfList[[i]]$Years = dfList[[i]]$GEN*4
}

#Merge to one data frame and order data 
allPopsDF = bind_rows(dfList)
newOrder = c("boxer", "cocker_spaniel", "german_shepherd_dog", "golden_retriever", "labrador_retriever", "maltese", "newfoundland", "poodle", "rottweiler", "yorkshire_terrier", "mixed", "village_dog_peru", "grayWolf_Europe", "grayWolf_NorthAmerica")
allPopsDF$Population = factor(allPopsDF$Population, levels = newOrder)
#Plot
#colourCount_breed = length(pops)
#palette = distinctColorPalette(colourCount_breed)

#All on one
allTogether = ggplot(allPopsDF, aes(x=Years, y=NE, colour=Population)) + 
  geom_line(size=1) + 
  scale_y_log10() + 
  scale_colour_manual(values = c(boxer="#D6E1A2",cocker_spaniel ="#D1EB48", grayWolf_NorthAmerica="#E350D4",german_shepherd_dog="#DDC552", golden_retriever="#7EE7C2", labrador_retriever ="#749583", maltese = "#DB5265", mixed = "#D7DDD5", newfoundland = "#867BCF", rottweiler = "#867BCF", village_dog_peru = "#8447E4", grayWolf_Europe = "#8CDC83", poodle = "#7BD9E2", yorkshire_terrier = "#DAAAC7")) + 
  theme_bw() + 
  theme(axis.text.x = element_text( hjust= 0.5, vjust=1,size=20), 
        axis.text.y = element_text(size =20), 
        plot.title=element_text(size =24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=24), 
        legend.title=element_blank(), 
        legend.text=element_text(size=14), 
        legend.position = "bottom") 

#Split onto individual
separate = ggplot(allPopsDF, aes(x=Years, y=NE, colour=Population)) + 
  geom_line(size=1) + 
  geom_ribbon(aes(ymin=LWR.95.CI, ymax=UPR.95.CI), alpha=0.2) + 
  facet_wrap(~ Population, ncol = 5) + 
  scale_y_log10() + 
  scale_colour_manual(values = c(boxer="#D6E1A2",cocker_spaniel ="#D1EB48",grayWolf_NorthAmerica="#E350D4",german_shepherd_dog="#DDC552", golden_retriever="#7EE7C2", labrador_retriever ="#749583", maltese = "#DB5265", mixed = "#D7DDD5", newfoundland = "#867BCF", rottweiler = "#867BCF", village_dog_peru = "#8447E4", grayWolf_Europe = "#8CDC83", poodle = "#7BD9E2", yorkshire_terrier = "#DAAAC7")) + 
  theme_bw() + 
  theme(axis.text.x = element_text( hjust= 0.5, vjust=1,size=20), 
        axis.text.y = element_text(size =20), 
        plot.title=element_text(size =24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=24), 
        legend.title=element_blank(), 
        legend.text=element_text(size=14), 
        legend.position = "bottom") 

#Plot 
print(allTogether)
print(separate)

plotTogether = arrangeGrob(allTogether  + 
                             labs(x = NULL, y = NULL) + 
                             theme(legend.position="none"), separate + 
                             labs(x = NULL, y = NULL) +
                             theme(axis.text.x = element_text(angle = 40, vjust=0.8,size=20), 
                                   legend.position="none"), 
                           ncol = 2,
                           left=textGrob("Effective population (Ne)", gp=gpar(fontface="bold",fontsize=20), rot=90), 
                           bottom=textGrob("Time (years ago)", gp=gpar(fontface="bold",fontsize=20), 
                                           vjust=0.3), widths=c(3/4,2)) #legend on bottom

legendIBDNe = gtable_filter(ggplotGrob(allTogether), "guide-box") #pull universal legend
grid.arrange(arrangeGrob(plotTogether, legendIBDNe, 
                         heights=unit.c(unit(1, "npc") - legendIBDNe$heights, legendIBDNe$heights), 
                         ncol=1))

###Larger labels remove mixed breed dogs for figure 1
newNames <- c("boxer"="boxer","cocker_spaniel" ="cocker spaniel","grayWolf_NorthAmerica"="North American wolf","german_shepherd_dog"="german shepherd", "golden_retriever"="golden retriever", "labrador_retriever" ="labrador retriever", "maltese" = "maltese", "mixed" = "mixed", "newfoundland" = "newfoundland", "rottweiler" = "rottweiler", "village_dog_peru" = "village dog", "grayWolf_Europe" = "European wolf", "poodle" = "poodle", "yorkshire_terrier" = "yorkshire terrier")

Figure1 = ggplot(allPopsDF %>% filter(Population != "mixed"), aes(x=Years, y=NE, colour=Population)) + 
  geom_line(size=1) + 
  geom_ribbon(aes(ymin=LWR.95.CI, ymax=UPR.95.CI), alpha=0.2) + 
  facet_wrap(~ Population, labeller = as_labeller(newNames), ncol = 4) + 
  scale_y_log10() + 
  scale_colour_manual(values = c(boxer="#D6E1A2",cocker_spaniel ="#D1EB48",grayWolf_NorthAmerica="#E350D4",german_shepherd_dog="#DDC552", golden_retriever="#7EE7C2", labrador_retriever ="#749583", maltese = "#DB5265", mixed = "#D7DDD5", newfoundland = "#867BCF", rottweiler = "#867BCF", village_dog_peru = "#8447E4", grayWolf_Europe = "#8CDC83", poodle = "#7BD9E2", yorkshire_terrier = "#DAAAC7")) + 
  theme_bw() + 
  labs(x = "Years", y = bquote('Effective Population Size'~(N[E]))) + 
  theme(axis.text.x = element_text( hjust= 0.5, vjust=1,size=20), 
        axis.text.y = element_text(size =20), 
        plot.title=element_text(size =24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=24),
        strip.text = element_text(size=20),
        legend.title=element_blank(), 
        legend.text=element_text(size=14), 
        legend.position = "none") 

