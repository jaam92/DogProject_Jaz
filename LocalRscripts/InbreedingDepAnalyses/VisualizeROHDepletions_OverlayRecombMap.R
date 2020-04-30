#Libraries
library(data.table)
library(tidyverse)


#Empty list for recom map
recomMapList = list() 

#Process recombination map files that have a header
filenames = list.files(path="~/Documents/DogProject_Jaz/LocalRscripts/dog_genetic_maps/rmChromCol/", pattern = "^chr*")
for(i in seq_along(filenames)){
  #Compute the average recombination rate per bin using the current bin and previous bin info
  df = read.delim(file = paste0("~/Documents/DogProject_Jaz/LocalRscripts/dog_genetic_maps/rmChromCol/",filenames[i]), stringsAsFactors = F) %>%
    mutate(chrom = parse_number(filenames[i]),
           avgPerBin = (cM - lag(cM, default = NA))/(pos - lag(pos, default = NA))) %>%
    select(chrom, pos, avgPerBin)
  df[is.na(df)] <- as.numeric(0) #replace intial NA with 0
  #Append dataframe to list
  recomMapList[[i]] <- df 
}
recomMapAllChroms = as.data.frame(rbindlist(recomMapList)) %>%
  arrange(chrom)

#Now process the empirical data
#Read empirical data in
df = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/CountEmpiricalOverlaps_100Kb_AutosomalSplits.bed", stringsAsFactors = F) %>%
  mutate(chrom = as.numeric(gsub("chr", "", V1)),
         proportion = V4/4342) %>%
  arrange(chrom) 

avgPerChrom = df %>%
  group_by(chrom) %>%
  summarise(mean = mean(proportion)) 


exons = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/ExonRegion_NonOverlapsROH_cornellData.bed", stringsAsFactors = F) %>% 
  mutate(chrom = as.numeric(gsub("chr", "", V1)), 
         proportion = as.numeric(0.22), 
         indicator = V4,
         dummy = gsub("ANKH|FYTTD1|PRMT2", "1", indicator),
         data = ifelse(dummy == "1", as.character(V4), "no ROH")) %>%
  select(-c(indicator, dummy))

foldChange = df %>%
  group_by(chrom) %>%
  mutate(adjV4 = V4 + 0.05,
         proportion = adjV4/mean(adjV4)) %>%
  select(-c(adjV4))

foldChange_exons = exons %>%
  mutate(proportion = as.numeric(0.9))

bindDF_foldChange = foldChange %>%
  mutate(data="empirical") %>%
  rbind.data.frame(foldChange_exons) %>%
  na.omit()

#Everything on different scale
ggplot(data = recomMapAllChroms, aes(x=pos, y=avgPerBin*1e5)) +
  geom_line(colour = "orange") +
  geom_point(data = bindDF_foldChange, aes(x=V3, y=proportion, colour=data)) +
  facet_wrap(~chrom, scales = "free") +
  scale_color_manual(values = c("empirical" = "black", "ANKH" = "blue",  "exon" = "steelblue", "FYTTD1" = "gray50", "no ROH" = "purple", "PRMT2" = "red")) +  
  labs(x="Position along chromosome (base pairs)", y="Fold-change overlapping ROHs (relative to chromosomal mean)") +
  theme_bw()  + 
  theme(axis.text.x = element_text(hjust= 0.8, vjust=0.8, angle = 45, size=14), 
        axis.text.y = element_text(size = 16), 
        plot.title=element_text(size = 20, face = "bold", hjust=0.5), 
        axis.title=element_text(size= 18),
        strip.text = element_text(size= 14),
        legend.title=element_text(size= 18), 
        legend.text=element_text(size= 16)) + ylim(0,2)
