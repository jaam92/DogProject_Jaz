#Load Libraries 
suppressPackageStartupMessages(library(tidyverse))
library(optparse)
library(tictoc)

#Input parameters from command line if not found set to default
option_list <- list( 
  make_option(c("--SGETaskID"), type="integer", 
              help="SGE Task ID for job array"),
  make_option(c("--genomeFile"), type="character", default="/u/scratch/j/jmooney3/PermuteROH/test/chromosomeLengths.txt", 
              help="file with chrom number in one column and length of chrom in second column", metavar="character"),
  make_option(c("--rohInfile"), type="character", default="/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/InbreedingDepression/HaywardData/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFile_Cornell_allChroms_vcfToolsROH_rmROHlessThan50snps.txt", 
              help="file with true ROHs in it from my pipeline", metavar="character"),
  make_option(c("--outFilePath"), type="character", default="/u/scratch/j/jmooney3/",
              help="specify output file path", metavar="character"))

#Parse the options
opt = parse_args(OptionParser(option_list=option_list))
SGETaskID = opt$SGETaskID
outFilePath = opt$outFilePath
#autosome = read.delim(file = opt$genomeFile, check.names = F, stringsAsFactors = F, sep = " ") 
#rohs = read.delim(file = opt$rohInfile, stringsAsFactors = F)

#Read autosome length and roh files in 
autosome = read.delim(file = "~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/chromosomeLengths.txt", check.names = F, stringsAsFactors = F, sep = " ")  #%>%
  #mutate(end = cumsum(as.numeric(LENGTH)),
  #       start = lag(end, default = 0),
  #       check = end - start) 
rohs = read.delim(file = "~/Documents/DogProject_Jaz/LocalRscripts/ROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt", stringsAsFactors = F)

#Fxn for rounding 
round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}

#Fxn to turn each chromosome into a vector of 0 or number aka nonROH/ROH and then permute the location of ROH with samp
permuteROH = function(dataFrame){
  totalNonROH = sum(autosome$LENGTH) - sum(dataFrame$AUTO_LEN)
  #calculate the total number of 100 bp nonROH blocks by rounding totalNonROH to the nearest 100 then dividing by 100
  total_nonROHBlocks = floor(round_any(x = totalNonROH, accuracy = 5) / 100)
  #make a vector of the chromosome where 0 is each 100bp nonROH block and ROHs are represented as their lengths
  recodeChrom = c(rep(0, total_nonROHBlocks), dataFrame$AUTO_LEN)
  #permute rohs with sample  
  #recode with rle format and make a data frame
  #convert data frame to chromosome coords
  #multiply number of 0 by 100 since we spaced every 100 bp else leave alone
  #start is where previous element ends and the end is cummulative sum minus 1 because our default on lag is 0-based
  #output roh only
  permuted = data.frame(unclass(rle(sample(recodeChrom)))) %>% 
    mutate(lengths = ifelse(values == "0", lengths*100, values),
           end = cumsum(as.numeric(lengths)),
           start = lag(end, default = 0),
           check = end - start) %>% #the check is to make sure the start and end lengths are correct
    filter(values != "0") %>%
    mutate(chrom = "chr1") %>%
    select(chrom, start, end)
}

#Permute ROH for all individuals and on all chroms
tic("total time all indivs and chroms")

mylist = rohs %>%
  group_by(INDV) %>%
  group_map(~ permuteROH(dataFrame = .x))

lapply(mylist, function(x) write.table(x,"~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/test.bed", quote = F, sep = "\t", row.names = F, col.names = F, append = T))


toc()
