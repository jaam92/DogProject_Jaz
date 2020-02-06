#Load library
library("tidyverse")

#Load files
popmapMerge = read.delim("~/DogProject_Jaz/LocalRscripts/BreedCladeInfo/BreedAndCladeInfo_mergedFitakCornell.txt")
rohs = read.delim("~/DogProject_Jaz/LocalRscripts/ROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt") %>%
  group_by(INDV) %>%
  summarise(ROHBurden = sum(AUTO_LEN))

#Grab case-control data
psva = read.delim("~/DogProject_Jaz/LocalRscripts/CaseControlROH/splitPhenotypeFile/PSVA_cairn_terrier.txt") %>%
  mutate(ROH = rohs$ROHBurden[match(dogID, rohs$INDV)])

ggplot(psva, aes(x=PSVA, y=ROH, group=PSVA, colour=PSVA)) +
  geom_violin() +
  geom_jitter() +
  theme_bw() +
  labs(x ="Case-Control Status", y="ROH Burden") +
  scale_x_continuous(breaks = c(1,2), labels = c("1"="Controls", "2"="Cases")) +
  scale_color_gradient(low = "blue",high = "red") +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 18), 
        plot.title=element_text(size=24, face = "bold", hjust=0.5), 
        axis.title=element_text(size=20),
        legend.position = "none") +
  stat_summary(fun.y = "median", geom = "point", size=2, color="black")

