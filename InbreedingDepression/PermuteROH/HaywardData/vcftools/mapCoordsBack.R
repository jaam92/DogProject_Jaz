#load libraries and turn off scientific notation
library(tidyverse)
library(data.table)
options(scipen=999)
library(optparse)

#Input parameters from command line if not found set to default
option_list <- list( 
  make_option(c("--SGETaskID"), type="integer", 
              help="SGE Task ID for job array"),
  make_option(c("--genomeFile"), type="character", default="/u/scratch/j/jmooney3/PermuteROH/chromosomeLengths.txt", 
              help="file with chrom number in one column and length of chrom in second column", metavar="character"),
  make_option(c("--rohInfile"), type="character", default=NULL, 
              help="file with new ROHs in it from my pipeline", metavar="character"),
  make_option(c("--outFilePath"), type="character", default="/u/scratch/j/jmooney3/",
              help="specify output file path", metavar="character"))

#Parse the options
opt = parse_args(OptionParser(option_list=option_list))
SGETaskID = opt$SGETaskID
outFilePath = opt$outFilePath
autosome = read.delim(file = opt$genomeFile, check.names = F, stringsAsFactors = F, sep = " ")  %>%
  mutate(end = cumsum(as.numeric(LENGTH)),
         start = lag(end, default = 0),
         check = end - start) 
df = read.table(file = opt$rohInfile)

#Empty list and dataframe to hold hash table and updated roh coords
LookupList = list()
rohList = list()

#Fxn to merge rows where start is on one chrom and end is on another
coalesce_by_column <- function(df) {
  return(dplyr::coalesce(!!! as.list(df)))
}

#Read autosome length and roh files in 
#autosome = read.delim(file = "../chromosomeLengths.txt", check.names = F, stringsAsFactors = F, sep = " ")  %>%
  #mutate(end = cumsum(as.numeric(LENGTH)),
  #       start = lag(end, default = 0),
  #       check = end - start) 

#df = read.table("newCoords1.bed") #read in coordinates

#make the hashtable
for(i in 1:38){
  hash = cbind.data.frame(c(0:autosome[i,2]), c(autosome[i,4]:autosome[i,3])) %>%
    mutate(chrom = paste0("chr",i)) 
  LookupList[[i]] <- hash
}

#map the new coordinates back to their original positions on a non-plasmid like genome
for (i in seq_along(LookupList)) {
  #start mapping
  df$mapStart = LookupList[[i]][,1][match(df$V2,LookupList[[i]][,2])] #og chrom start
  df$mapEnd = LookupList[[i]][,1][match(df$V3,LookupList[[i]][,2])] #og chrom end
  df$mapStartChrom = LookupList[[i]][,3][match(df$V2,LookupList[[i]][,2])] #og chrom
  df$mapEndChrom = LookupList[[i]][,3][match(df$V3,LookupList[[i]][,2])] #og chrom
  x = df[rowSums(is.na(df)) < 3, ] #remove rows with 3 more NAs, want to allow 2 NAs for those ROH that start on one chrom and end on another
  rohList[[i]] <- x #append dataframe to list 
}

FinalDF = as.data.frame(rbindlist(rohList)) %>% 
  group_by(V2,V3) %>% 
  summarise_all(coalesce_by_column) %>%
  ungroup()

#check whether start and end og chroms match
FinalDF$match = ifelse(FinalDF$mapEndChrom == FinalDF$mapStartChrom, FinalDF$mapStartChrom, paste(FinalDF$mapStartChrom, FinalDF$mapEndChrom,sep=","))
FinalDF$updateStartEnd = ifelse(FinalDF$mapEndChrom == FinalDF$mapStartChrom, "N", paste("E","S", sep=",")) #indicator for what will be replaced once rows are separated needs to follow order start and end chroms are given in paste for match variable

#split the dataframe into separate rows based on match and update cols
FinalDF2 = FinalDF %>%
  ungroup() %>%
  separate_rows(match, updateStartEnd, convert = T) 

#update starts and ends
FinalDF2$mapEnd = ifelse(FinalDF2$updateStartEnd == "E", autosome$LENGTH[match(FinalDF2$mapStartChrom, autosome$CHROM)], FinalDF2$mapEnd)
FinalDF2$mapStart = ifelse(FinalDF2$updateStartEnd == "S", as.numeric(0), FinalDF2$mapStart)

#subset down and output
writeOutDF = FinalDF2 %>%
  select(match, mapStart, mapEnd)

#write out
write.table(writeOutDF, paste0(outFilePath,"newCoords",SGETaskID,"_mappedBack.bed"), quote = F, row.names = F, col.names = F, sep = "\t")
