#Load libraries
library(tidyverse)

#####Only ran this commented section one time to generate the bedtools file where coordinates were mapped back to the original from the plasmind like autosome coordinates, saved new file to output file #####
#library(data.table)
#Empty list and dataframe to hold hash table and updated roh coords
#LookupList = list()
#rohList = list()
#Load files
#autosome = read.delim(file = "~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/chromosomeLengths.txt", check.names = F, stringsAsFactors = F, sep = " ")  %>%
#  mutate(end = cumsum(as.numeric(LENGTH)),
#         start = lag(end, default = 0),
#         check = end - start) 
#perms = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/SummaryFile_Merge_CountPermutedOverlaps_100Kb_AutosomalSplits_bedtoolsShuffle.txt") 
#make the hashtable
#for(i in 1:38){
#  hash = cbind.data.frame(c(0:autosome[i,2]), c(autosome[i,4]:autosome[i,3])) %>%
#    mutate(chrom = paste0("chr",i)) 
#  LookupList[[i]] <- hash
#}
#map the new coordinates back to their original positions on a non-plasmid like genome
#for (i in seq_along(LookupList)) {
#start mapping
#  perms$mappedBackStart = LookupList[[i]][,1][match(perms$start,LookupList[[i]][,2])] #og chrom start
#  perms$mappedBackChromStart = LookupList[[i]][,3][match(perms$start,LookupList[[i]][,2])] #og chrom
#  x = perms[rowSums(is.na(perms)) < 1, ] #remove rows with 1 more NAs, want to allow 2 NAs for those ROH that start on one chrom and end on another
#  rohList[[i]] <- x #append dataframe to list 
#}
#merge everything back together
#FinalDF = as.data.frame(rbindlist(rohList)) 
#write out
#write.table(FinalDF, file = "~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/SummaryFile_Merge_CountPermutedOverlaps_100Kb_AutosomalSplits_bedtoolsShuffle_mappedBack.txt", quote = F, row.names = F, col.names = T, sep = "\t")

#Read empirical data in
df = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/CountEmpiricalOverlaps_100Kb_AutosomalSplits.bed") %>%
  mutate(chrom = as.numeric(gsub("chr", "", V1)),
         proportion = V4/4342) %>%
  arrange(chrom) 

perms = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/SummaryFile_Merge_CountPermutedOverlaps_100Kb_AutosomalSplits_bedtoolsShuffle.txt") 

avgPerChrom = df %>%
  group_by(chrom) %>%
  summarise(mean = mean(proportion)) 

#Look at the distribution of empirical data overlaps
ggplot(df, aes(x=V4)) +
  geom_density() +
  labs(x="Number of ROH overlaps\nper 100Kb segment") +
  theme_bw() + 
  theme(axis.text.x = element_text(size=18), 
        axis.text.y = element_text(size =18), 
        plot.title=element_text(size = 22, face = "bold", hjust=0.5), 
        axis.title=element_text(size= 20)) 

#Everything on same scale
ggplot(data = df, aes(x=V3,y=proportion)) +
  geom_point() +
  geom_hline(aes(yintercept=mean(proportion), linetype = "Genome-wide average"), colour= 'red') +
  geom_hline(data = avgPerChrom, aes(yintercept=mean, linetype="Chromosome average"), colour= 'blue') +
  facet_wrap(~chrom) +
  scale_linetype_manual(name = "Region", values = c(1, 2), 
                        guide = guide_legend(override.aes = list(color = c("blue", "red")))) +
  scale_y_continuous(labels=scales::percent) +
  
  labs(x="Position along chromosome (base pairs)", y="Number of overlapping ROHs") +
  theme_bw()  + 
  theme(axis.text.x = element_text(hjust= 0.8, vjust=0.8, angle = 45, size=14), 
        axis.text.y = element_text(size = 16), 
        plot.title=element_text(size = 20, face = "bold", hjust=0.5), 
        axis.title=element_text(size= 18),
        strip.text = element_text(size= 14),
        legend.title=element_text(size= 18), 
        legend.text=element_text(size= 16)) 

#Everything on different scale
ggplot(data = df, aes(x=V3,y=proportion)) +
  geom_point() +
  geom_hline(aes(yintercept=mean(proportion), linetype = "Genome-wide average"), colour= 'red') +
  geom_hline(data = avgPerChrom, aes(yintercept=mean, linetype="Chromosome average"), colour= 'blue') +
  facet_wrap(~chrom, scales = "free") +
  scale_linetype_manual(name = "Region", values = c(1, 2), 
                        guide = guide_legend(override.aes = list(color = c("blue", "red")))) +
  scale_y_continuous(labels=scales::percent) +
  
  labs(x="Position along chromosome (base pairs)", y="Number of overlapping ROHs") +
  theme_bw()  + 
  theme(axis.text.x = element_text(hjust= 0.8, vjust=0.8, angle = 45, size=14), 
        axis.text.y = element_text(size = 16), 
        plot.title=element_text(size = 20, face = "bold", hjust=0.5), 
        axis.title=element_text(size= 18),
        strip.text = element_text(size= 14),
        legend.title=element_text(size= 18), 
        legend.text=element_text(size= 16)) 

#plot the permutations
ggplot(data = perms, aes(x=start,y=avg)) +
  geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.1) +
  geom_point(data = df, aes(x=V3,y=proportion), colour="red") +
  facet_wrap(~chrom, scales = "free") +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Position along chromosome (base pairs)", y="Number of overlapping ROHs") +
  theme_bw()  + 
  theme(axis.text.x = element_text(hjust= 0.8, vjust=0.8, angle = 45, size=14), 
        axis.text.y = element_text(size = 16), 
        plot.title=element_text(size = 20, face = "bold", hjust=0.5), 
        axis.title=element_text(size= 18),
        strip.text = element_text(size= 14),
        legend.title=element_text(size= 18), 
        legend.text=element_text(size= 16)) 
