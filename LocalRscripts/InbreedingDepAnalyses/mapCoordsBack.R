#load libraries and turn off scientific notation
library(tidyverse)
library(data.table)
options(scipen=999)

#read in hashtables
files = list.files(path = "/u/scratch/j/jmooney3/PermuteROH/test",
                   pattern = glob2rx("hash*.txt$"))
hashTable = rbindlist(sapply(files, fread, simplify = FALSE))

#read in autosome
autosome = read.delim(file = "~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/chromosomeLengths.txt", check.names = F, stringsAsFactors = F, sep = " ")

#read in coordinates
df = read.table("newCoords1.bed")

#start mapping
df$mapStart = hashTable$V1[match(df$V2,hashTable$V2)] #og chrom start
df$mapEnd = hashTable$V1[match(df$V3,hashTable$V2)] #og chrom end
df$mapStartChrom = hashTable$V3[match(df$V2,hashTable$V2)] #og chrom
df$mapEndChrom = hashTable$V3[match(df$V3,hashTable$V2)] #og chrom

#check whether start and end og chroms match
df$match = ifelse(df$mapEndChrom == df$mapStartChrom, df$mapStartChrom, paste(df$mapStartChrom, df$mapEndChrom,sep=","))

#indicator for what will be replaced once rows are separated needs to follow order start and end chroms are given in paste for match variable
df$updateStartEnd = ifelse(df$mapEndChrom == df$mapStartChrom, "N", paste("E","S", sep=","))

#split the dataframe into separate rows based on match and update cols
df2 = df %>%
  separate_rows(match, updateStartEnd, convert = T)

#update starts and ends
df2$mapEnd = ifelse(df2$updateStartEnd == "E", autosome$LENGTH[match(df2$mapStartChrom, autosome$CHROM)], df2$mapEnd)
df2$mapStart = ifelse(df2$updateStartEnd == "S", as.numeric(0), df2$mapStart)

#subset down and output
finalDF = df2 %>%
  select(match, mapStart, mapEnd)

#write out
write.table(finalDF, "mapBackNewCoords.bed", quote = F, row.names = F, col.names = F, sep = "\t")
