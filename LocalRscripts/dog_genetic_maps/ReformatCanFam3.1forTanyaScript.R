setwd("~/Documents/DogProject_Jaz/dog_genetic_maps/rmChromCol/")

filenames = list.files(pattern ="chr*") 

for (i in 1:length(filenames))
{
  #Import file and give col names 
  dogGenMap = read.delim(filenames[i])
  
  #Find difference between current and previous site
  dogGenMap$diff = ave(dogGenMap$cM, FUN=function(x) c(0, diff(x)))
  #Remove the rate column 
  dogGenMap$rate = NULL
  #Rename position column
  names(dogGenMap)[1] = "Pos"
  #Generate file for Tanya's script 
  write.table(dogGenMap, file= paste("ReformattedTanyaScript_", filenames[i], sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
