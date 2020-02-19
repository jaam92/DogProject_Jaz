#Load libraries
library(tidyverse)
library(GenomicRanges)

#Load breed info, ROHs, and kinship matrix
rohs = read.delim(file = "~/Documents/DogProject_Jaz/LocalRscripts/ROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt", stringsAsFactors = F) %>% 
  filter(INDV == "30") 

autosome = read.delim(file = "~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/chromosomeLengths.txt", check.names = F, stringsAsFactors = F, sep = " ") %>%
  mutate(CHROM = as.numeric(gsub("chr", "", CHROM)))

autoLen = sum(as.numeric(autosome$LENGTH))

#Move runs on same chrom by sampling uniformly
newCoords = data.frame() #keep first set of rohs
updatedROHs = data.frame() #keep rohs when overlaps are detected
finalROH = data.frame() #permuted rohs

for(i in rohs$AUTO_LEN){
  rohLength = i
  start = runif(1, min = 0, max = autoLen)
  end = start + rohLength
  #Next check whether rohs fall within a chromosome, if they do not make new roh until it does
    if (end <= autoLen) {
      newROH = cbind(start, end)
      newCoords = rbind(newCoords, newROH)
    }else{
      while (end > autoLen) {
        start = runif(1, min = 0, max = autoLen)
        end = start + rohLength
                            }
      newROH = cbind(start, end)
      newCoords = rbind(newCoords, newROH)
    }
  
  #Check for overlapping rohs if there are overlaps redo everything from the top until there are none
  finalROH = newCoords %>%
    mutate(start = start/10^3,
             end = end/10^3) #scale everything to kb so that integers in IRanges don't overflow
  finalROH$chrom = "chr1"
  rohRanges = with(finalROH, GRanges(chrom, IRanges(start=start, end =end)))
  noOverlaps = rohRanges[unique(findOverlaps(rohRanges, type = "any", select = "first"))] %>%
    as.data.frame() %>%
    select(start, end)
  
  if (nrow(finalROH) == nrow(newCoords)) {
    cat(sprintf("new set rohs created"), sep = "\n")

  }else{
    while (nrow(finalROH) < nrow(newCoords)) {
      cat(sprintf("overlaps detected generating new set of rohs"), sep = "\n")
      
      #Re-run everything from the top
      start = runif(1, min = 0, max = autoLen)
      end = start + rohLength
      if (end <= autoLen) {
        retryROH = cbind(start, end)
        updatedROHs = rbind(finalROH, retryROH)
      }else{
        while (end > autoLen) {
          start = runif(1, min = 0, max = autoLen)
          end = start + rohLength
        }
        retryROH = cbind(start, end)
        updatedROHs = rbind(finalROH, retryROH)
        
      }
      #Check for overlaps again
      makeNewROHs = updatedROHs %>%
        mutate(chrom = "chr1",
               start = start/10^3,
               end = end/10^3) #scale everything to kb so that integers in IRanges don't overflow
      rohRanges = with(makeNewROHs, GRanges(chrom, IRanges(start=start, end =end)))
      finalROH = rohRanges[unique(findOverlaps(rohRanges, type = "any", select = "first"))] %>%
        as.data.frame() %>%
        select(start, end)
    }
    
  }
#print(finalROH)
}

autosome$START = lag(autosome$LENGTH, default = 0)+1