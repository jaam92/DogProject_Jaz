#Load libraries
library(tidyverse)

#Read empirical data in
df = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/CountEmpiricalOverlaps_100Kb_AutosomalSplits.bed") %>%
  mutate(chrom = as.numeric(gsub("chr", "", V1)),
         proportion = V4/4342) %>%
  arrange(chrom) 

perms = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/SummaryFile_Merge_CountPermutedOverlaps_100Kb_AutosomalSplits_chromFirst.txt") 

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
  geom_hline(aes(yintercept=mean(df$proportion), linetype = "Genome-wide average"), colour= 'red') +
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
  geom_hline(aes(yintercept=mean(df$proportion), linetype = "Genome-wide average"), colour= 'red') +
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
