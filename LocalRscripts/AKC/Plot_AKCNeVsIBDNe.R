setwd("~/Documents/DogProject_Jaz/IBDNe/AllSitesMerge/")
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

#Read in AKC data
AKC_transposed = read.delim("~/Documents/DogProject_Jaz/AKC/AKC_1926thru2005_NeTransposed.txt", check.names = F)

#pdf("~/Documents/DogProject_Jaz/AKC/AKC_vs_IBDNe.pdf")

Plot1 = ggplot() + geom_point(data = AKC_transposed, aes(x=abs(Time-2005), y= log10(boxer), colour="red")) + geom_point(data=boxer, aes(x=GEN*4, y=log10(NE), colour="blue")) + xlim(0,200)  + scale_colour_manual(name="Data Source", values=c( blue="blue",red="red"), label=c("IBDNe","AKC")) + labs(x="Time (Years Ago)", y="Log10(Population Size)") + ggtitle("boxer") + theme_bw()+ theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24),legend.title=element_text(size=24), legend.text=element_text(size=18))

Plot2 = ggplot() + geom_point(data = AKC_transposed, aes(x=abs(Time-2005), y= log10(cocker_spaniel), colour="red")) + geom_point(data=cocker_spaniel, aes(x=GEN*4, y=log10(NE), colour="blue")) + xlim(0,200)  + scale_colour_manual(name="Data Source", values=c( blue="blue",red="red"), label=c("IBDNe","AKC")) + labs(x="Time (Years Ago)", y="Log10(Population Size)") + ggtitle("cocker_spaniel") + theme_bw()+ theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24),legend.title=element_text(size=24), legend.text=element_text(size=18))

Plot3 = ggplot() + geom_point(data = AKC_transposed, aes(x=abs(Time-2005), y= log10(german_shepherd_dog), colour="red")) + geom_point(data=german_shepherd_dog, aes(x=GEN*4, y=log10(NE), colour="blue")) + xlim(0,200)  + scale_colour_manual(name="Data Source", values=c( blue="blue",red="red"), label=c("IBDNe","AKC")) + labs(x="Time (Years Ago)", y="Log10(Population Size)") + ggtitle("german_shepherd_dog") + theme_bw()+ theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24),legend.title=element_text(size=24), legend.text=element_text(size=18))

Plot4 = ggplot() + geom_point(data = AKC_transposed, aes(x=abs(Time-2005), y= log10(labrador_retriever), colour="red")) + geom_point(data=labrador_retriever, aes(x=GEN*4, y=log10(NE), colour="blue")) + xlim(0,200)  + scale_colour_manual(name="Data Source", values=c( blue="blue",red="red"), label=c("IBDNe","AKC")) + labs(x="Time (Years Ago)", y="Log10(Population Size)") + ggtitle("labrador_retriever") + theme_bw()+ theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24),legend.title=element_text(size=24), legend.text=element_text(size=18))

Plot5 = ggplot() + geom_point(data = AKC_transposed, aes(x=abs(Time-2005), y= log10(golden_retriever), colour="red")) + geom_point(data=golden_retriever, aes(x=GEN*4, y=log10(NE), colour="blue")) + xlim(0,200)  + scale_colour_manual(name="Data Source", values=c( blue="blue",red="red"), label=c("IBDNe","AKC")) + labs(x="Time (Years Ago)", y="Log10(Population Size)") + ggtitle("golden_retriever") + theme_bw()+ theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24),legend.title=element_text(size=24), legend.text=element_text(size=18))

Plot6 = ggplot() + geom_point(data = AKC_transposed, aes(x=abs(Time-2005), y= log10(maltese), colour="red")) + geom_point(data=maltese, aes(x=GEN*4, y=log10(NE), colour="blue")) + xlim(0,200)  + scale_colour_manual(name="Data Source", values=c( blue="blue",red="red"), label=c("IBDNe","AKC")) + labs(x="Time (Years Ago)", y="Log10(Population Size)") + ggtitle("maltese") + theme_bw()+ theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24),legend.title=element_text(size=24), legend.text=element_text(size=18))

Plot7 = ggplot() + geom_point(data = AKC_transposed, aes(x=abs(Time-2005), y= log10(newfoundland), colour="red")) + geom_point(data=newfoundland, aes(x=GEN*4, y=log10(NE), colour="blue")) + xlim(0,200)  + scale_colour_manual(name="Data Source", values=c( blue="blue",red="red"), label=c("IBDNe","AKC")) + labs(x="Time (Years Ago)", y="Log10(Population Size)") + ggtitle("newfoundland") + theme_bw()+ theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24),legend.title=element_text(size=24), legend.text=element_text(size=18))

Plot8 = ggplot() + geom_point(data = AKC_transposed, aes(x=abs(Time-2005), y= log10(poodle), colour="red")) + geom_point(data=poodle, aes(x=GEN*4, y=log10(NE), colour="blue")) + xlim(0,200)  + scale_colour_manual(name="Data Source", values=c( blue="blue",red="red"), label=c("IBDNe","AKC")) + labs(x="Time (Years Ago)", y="Log10(Population Size)") + ggtitle("poodle") + theme_bw()+ theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24),legend.title=element_text(size=24), legend.text=element_text(size=18))

Plot9 = ggplot() + geom_point(data = AKC_transposed, aes(x=abs(Time-2005), y= log10(rottweiler), colour="red")) + geom_point(data=rottweiler, aes(x=GEN*4, y=log10(NE), colour="blue")) + xlim(0,200)  + scale_colour_manual(name="Data Source", values=c( blue="blue",red="red"), label=c("IBDNe","AKC")) + labs(x="Time (Years Ago)", y="Log10(Population Size)") + ggtitle("rottweiler") + theme_bw()+ theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24),legend.title=element_text(size=24), legend.text=element_text(size=18))

Plot10 = ggplot() + geom_point(data = AKC_transposed, aes(x=abs(Time-2005), y= log10(yorkshire_terrier), colour="red")) + geom_point(data=yorkshire_terrier, aes(x=GEN*4, y=log10(NE), colour="blue")) + xlim(0,200)  + scale_colour_manual(name="Data Source", values=c( blue="blue",red="red"), label=c("IBDNe","AKC")) + labs(x="Time (Years Ago)", y="Log10(Population Size)") + ggtitle("yorkshire_terrier") + theme_bw()+ theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24),legend.title=element_text(size=24), legend.text=element_text(size=18))

#dev.off()


plotAll = arrangeGrob(Plot1  + labs(x = NULL, y = NULL) + theme(legend.position="none"), Plot2 + labs(x = NULL, y = NULL) + theme(legend.position="none"), Plot3 + labs(x = NULL, y = NULL) + theme(legend.position="none"), Plot4 + labs(x = NULL, y = NULL) + theme(legend.position="none"), Plot5 + labs(x = NULL, y = NULL) + theme(legend.position="none"), Plot6 + labs(x = NULL, y = NULL) + theme(legend.position="none"), Plot7 + labs(x = NULL, y = NULL) + theme(legend.position="none"), Plot8 + labs(x = NULL, y = NULL) + theme(legend.position="none"), Plot9 + labs(x = NULL, y = NULL) + theme(legend.position="none"), Plot10 + labs(x = NULL, y = NULL) + theme(legend.position="none"),ncol = 5,left=textGrob("Effective population (Ne)", gp=gpar(fontface="bold",fontsize=20), rot=90), bottom=textGrob("Time (years ago)", gp=gpar(fontface="bold",fontsize=20), vjust=0.3))

legend = gtable_filter(ggplotGrob(Plot1), "guide-box")

grid.arrange(plotAll,legend, widths=unit.c(unit(1, "npc") - legend$width, legend$width), nrow=1)
