#Load library
library(tidyverse)
library(randomcoloR)
library(RColorBrewer)

#Load files
popmapMerge = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/BreedCladeInfo/BreedAndCladeInfo_mergedFitakCornell.txt")
rohs = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/ROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt") %>%
  group_by(INDV) %>%
  summarise(ROHBurden = sum(AUTO_LEN))

#evaluate case-control data
fnames = paste0("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/splitPhenotypeFile/IncludeMixedBreeds/", grep(list.files(path="~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/splitPhenotypeFile/IncludeMixedBreeds/"), pattern='_', inv=T, value=T)) #paste the path in front of the filename and only include files with all breeds rather than those separated by breed

pdf(file ="~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/ROHBurden_CaseControl_includeMixedBreeds.pdf", height = 14, width = 28) #open pdf

for (i in seq_along(fnames)) {
  trait = gsub(".*[/]([^.]+)[.txt].*", "\\1", fnames[i]) #keep track of trait and breed tested
  
  traitDF = read.delim(fnames[i]) %>%
    mutate(ROH = rohs$ROHBurden[match(dogID, rohs$INDV)]) #make a data frame with trait of interest and ROH per indiv
  names(traitDF)[3] = "status" #rename the trait column
  
  #Get ready for plotting
  colourCount = length(unique(traitDF$breed))
  palette = distinctColorPalette(colourCount)
  
  plotROHBurden = ggplot(traitDF, aes(x=status, y=ROH/10^6, group=status, colour=status)) +
    geom_boxplot() +
    geom_point() +
    theme_bw() +
    facet_wrap(~breed, scale="free") +
    labs(x ="Case-Control Status", y="ROH Burden (Mb)", title = paste(trait)) +
    scale_x_continuous(breaks = c(1,2), labels = c("1"="Controls", "2"="Cases")) +
    scale_color_gradient(low = "blue", high = "red") +
    theme(axis.text.x = element_text(size = 18), 
          axis.text.y = element_text(size = 18), 
          plot.title= element_text(size=24, face = "bold", hjust=0.5), 
          axis.title= element_text(size=20),
          strip.text = element_text(size = 14),
          legend.position = "none") +
    stat_summary(fun.y = "median", geom = "point", size=2, color="black") #plot ROH burden separated by breed and case-control status
  
  #Density plots 
  plotROHBurdenDensities = ggplot() +
    geom_density(data = traitDF, aes(x=ROH/10^6, fill=breed),alpha=0.5, adjust=2) +
    theme_bw() +
    scale_fill_manual(values = palette, na.value = "black") +
    labs(x ="ROH Burden (Mb)", title = paste(trait)) +
    theme(axis.text.x = element_text(size = 18), 
          axis.text.y = element_text(size = 18), 
          plot.title= element_text(size=24, face = "bold", hjust=0.5), 
          axis.title= element_text(size=20),
          legend.position = "none")  
    
  #put in pdf
  print(plotROHBurdenDensities)
  print(plotROHBurden)
}


dev.off()
