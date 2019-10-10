#Load Libraries
library(tidyverse)
library(pracma)

#Read file in and reformat
setwd("~/Documents/DogProject_Jaz/LocalRscripts/AKC")

AKC = read.table("~/Documents/DogProject_Jaz/LocalRscripts/AKC/AKC.txt", check.names = F, header = T, fill = T, stringsAsFactors = F) %>%
  pivot_longer(-breed, names_to = "date", values_to = "count") %>% #reformat data frame
 filter(count != "N/A") %>% #remove NA
  group_by(breed) %>%
  slice(-1) %>% #remove the first occurance because there is a big registration influx when a group is recognized 
  mutate(count = as.numeric(count),
         date = as.numeric(date))

#Compute the percentile rank 
percRankPerYear = AKC %>%
  group_by(date) %>%
  mutate(perRank = percent_rank(count)) %>% #get the percentile rank of each breed present for a given year
  ungroup() 

#Popularity of breed
breedPop = percRankPerYear %>%
  group_by(breed) %>%
  summarise(popularity = trapz(x=date, y=perRank)) #for the trapz fxn make sure that the x variable is ordered ascending or the value will be negative

#Plot some quick comparisons
ggplot(percRankPerYear %>% 
         filter(breed=="golden_retriever" | breed=="beagle" | breed == "affenpinscher" | breed== "yorkshire_terrier"), 
       aes(x=date,y=perRank, colour=breed, group=breed)) +
  geom_line(size=1) +
  theme_bw()


#write out
write.table(breedPop, "AKC_breedPopularity_1926thru2005.txt", sep = "\t", row.names = F, col.names = T, quote = F)



