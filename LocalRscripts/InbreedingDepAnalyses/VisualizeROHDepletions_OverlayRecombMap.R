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
           avgPerBin = (cM - lag(cM, default = NA))/(pos - lag(pos, default = NA))*1e6) %>% #convert bp to Mb
    select(chrom, pos, avgPerBin)
  df[is.na(df)] <- as.numeric(0) #replace intial NA with 0
  #Append dataframe to list
  recomMapList[[i]] <- df 
}
recomMapAllChroms = as.data.frame(rbindlist(recomMapList)) %>%
  arrange(chrom)

#Now process the empirical data
#Read empirical data in
vcftools = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/CountEmpiricalOverlaps_100Kb_AutosomalSplits.bed", stringsAsFactors = F) %>%
  mutate(chrom = as.numeric(gsub("chr", "", V1)),
         proportion = V4/4342,
         data="VCFTools") %>%
  arrange(chrom) 

plink = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/CountEmpiricalOverlaps_100Kb_AutosomalSplits_plink.bed", stringsAsFactors = F) %>%
  mutate(chrom = as.numeric(gsub("chr", "", V1)),
         proportion = V4/4342,
         data="PLINK") %>%
  arrange(chrom) 

df = rbind.data.frame(vcftools, plink)

avgPerChrom = df %>%
  group_by(chrom, data) %>%
  summarise(mean = mean(proportion)) 

exons = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/VettingResults/ExonRegion_NonOverlapsROH_cornellData.bed", stringsAsFactors = F) %>% 
  mutate(chrom = as.numeric(gsub("chr", "", V1)), 
         proportion = as.numeric(0.22), 
         indicator = V4,
         dummy = gsub("ANKH|FYTTD1|PRMT2", "1", indicator),
         data = ifelse(dummy == "1", as.character(V4), "no ROH")) %>%
  select(-c(indicator, dummy))

foldChange = df %>%
  group_by(chrom, data) %>%
  mutate(adjV4 = V4 + 0.05,
         proportion = adjV4/mean(adjV4)) %>%
  select(-c(adjV4))

foldChange_exons = exons %>%
  mutate(proportion = as.numeric(0.9))

bindDF_foldChange = foldChange %>%
  rbind.data.frame(foldChange_exons) %>%
  na.omit() %>%
  rename("chromo" = V1 , "window_start"= V2 , "window_end" = V3, "overlaps"= V4 )

#Everything on different scale and split by chromosome
ggplot(data = bindDF_foldChange, aes(x=window_start, y=proportion, colour=data)) +
  geom_point() +
  facet_wrap(~chrom, scales = "free") +
  scale_color_manual(values = c("VCFTools" = "gray50", "PLINK" = "darkorange", "ANKH" = "blue",  "exon" = "gray50", "FYTTD1" = "steelblue", "no ROH" = "purple", "PRMT2" = "red")) +  
  labs(x="Position along chromosome (base pairs)", y="Fold-change overlapping ROHs (relative to chromosomal mean)") +
  theme_bw()  + 
  theme(axis.text.x = element_text(hjust= 0.8, vjust=0.8, angle = 45, size=14), 
        axis.text.y = element_text(size = 16), 
        plot.title=element_text(size = 20, face = "bold", hjust=0.5), 
        axis.title=element_text(size= 18),
        strip.text = element_text(size= 14),
        legend.title=element_text(size= 18), 
        legend.text=element_text(size= 16)) #+ 
  #ylim(0,2) 

#Plot separated by chromosome with recombination map
ggplot(data = recomMapAllChroms, aes(x=pos, y=avgPerBin)) +
  geom_line(colour = "orange") +
  geom_point(data = bindDF_foldChange, aes(x=window_start, y=proportion, colour=data)) +
  facet_wrap(~chrom, scales = "free") +
  scale_color_manual(values = c("VCFTools" = "gray50", "PLINK" = "darkorange", "ANKH" = "blue",  "exon" = "steelblue", "FYTTD1" = "gray50", "no ROH" = "purple", "PRMT2" = "red")) +  
  labs(x="Position along chromosome (base pairs)", y="Fold-change overlapping ROHs (relative to chromosomal mean)") +
  theme_bw()  + 
  theme(axis.text.x = element_text(hjust= 0.8, vjust=0.8, angle = 45, size=14), 
        axis.text.y = element_text(size = 16), 
        plot.title=element_text(size = 20, face = "bold", hjust=0.5), 
        axis.title=element_text(size= 18),
        strip.text = element_text(size= 14),
        legend.title=element_text(size= 18), 
        legend.text=element_text(size= 16)) #+ 
  #ylim(0,2)


#plot all the chromosomes together like the sliding window plots

#reformat data and rescale windows
rescaledWindows = bindDF_foldChange %>% 
  
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(window_start)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(bindDF_foldChange, ., by=c("chrom"="chrom")) %>%
  
  # Add a cumulative position of each window
  arrange(chrom, window_start) %>%
  mutate( newWinStart=window_start+tot)

#Function to plot data for each individual
plotFunction = function(dataFrame,dataWanted, color1, color2) {
  #filter data wanted
  newDataFrame = dataFrame %>%
    filter(data == dataWanted)
  #Generate x axis with any one of the data frames 
  axisdf = newDataFrame %>% 
    group_by(chrom) %>% 
    summarize(center=(max(newWinStart) + min(newWinStart) ) / 2 )
  #Now plot with the axis
  rohAlongChrom = ggplot() + 
    geom_point(data = newDataFrame, aes(x=newWinStart, y=proportion, color=as.factor(chrom))) +
    scale_color_manual(values = rep(c(color1, color2), 38 )) +
    #custom X axis:
    scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center ) +
    labs(x = "Chromosome", y = "Fold-change overlapping ROHs (relative to chromosomal mean)") +
    theme_bw() +
    theme(axis.text.x = element_text(size=16), 
          axis.text.y = element_text(size = 16), 
          axis.title=element_text(size= 18),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(), 
          legend.position = "none")
  return(rohAlongChrom)
}

rohAcrossGenome_VCFTools = plotFunction(rescaledWindows,"VCFTools", "gray50","darkmagenta")

rohAcrossGenome_PLINK = plotFunction(rescaledWindows,"PLINK", "darkorange","darkmagenta")
