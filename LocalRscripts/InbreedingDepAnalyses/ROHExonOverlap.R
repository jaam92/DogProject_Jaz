#Load libraries
library("dplyr")
library("ggplot2")

#set working directory and load files
setwd("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalysis")
plink = read.table("~/Documents/DogProject_Jaz/InbreedingDepression/plink/CountExonRegion_NonOverlapsROH.bed", stringsAsFactors = F, header = T)
vcftools = read.table("~/Documents/DogProject_Jaz/InbreedingDepression/vcftools/CountExonRegion_NonOverlapsROH_vcfTools.bed", stringsAsFactors = F, header = T)

#Check out overlaps
Overlaps = vcftools %>% 
  full_join(plink, by = c("GeneNames" = "GeneName"), suffix = c("_vcfTools", "_plink")) %>%
  mutate(inVCFTools = ifelse(is.na(CountExon_vcfTools),"0","1"),
         inPlink = ifelse(is.na(CountExon_plink),"0","1"),
         inBoth = ifelse(as.numeric(inPlink) + as.numeric(inVCFTools) == 2, "1", "0")) %>%
  select(GeneNames, everything()) 

#plot data
ggplot(plink, aes(x=reorder(GeneName, -CountExon), y=CountExon)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size  = 20,angle=40, vjust=1, hjust=1), axis.text.y = element_text(size  = 20), axis.title=element_text(size=24),legend.title=element_text(size=24), legend.text=element_text(size=18), plot.title = element_text(size  = 24,hjust=0.5)) +
  labs(x="Gene Name", y="Exon Count", title="Exons without ROH")

ggplot(vcftools, aes(x=reorder(GeneNames, -CountExon), y=CountExon))  + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size  = 20,angle=40, vjust=1, hjust=1), axis.text.y = element_text(size  = 20), axis.title=element_text(size=24),legend.title=element_text(size=24), legend.text=element_text(size=18), plot.title = element_text(size  = 24,hjust=0.5)) +
  labs(x="Gene Name", y="Exon Count", title="Exons without ROH")

#Merged bar plot
plotDF = Overlaps %>%
  select(GeneNames, CountExon_vcfTools, CountExon_plink) %>%
  melt() %>%
  mutate(variable = gsub("CountExon_", "", variable),
         value = ifelse(is.na(value), 0, value)) 

ggplot(data = plotDF, aes(x=reorder(GeneNames, -value), y= value, fill=variable)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c(plink = "blue", vcfTools = "orange"), name = "Method") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 16, angle=40, vjust=1, hjust=1), axis.text.y = element_text(size  = 20), axis.title=element_text(size=24),legend.title=element_text(size=24), legend.text=element_text(size=18), plot.title = element_text(size  = 24,hjust=0.5)) +
  labs(x="Gene Name", y="Exon Count", title="Exons without ROH")
