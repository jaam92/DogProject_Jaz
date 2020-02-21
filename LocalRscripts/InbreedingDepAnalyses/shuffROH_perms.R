#Load libraries
library(tidyverse)
library(GenomicRanges)

#Empty data frames for shuffling roh
newCoords = data.frame() #keep first set of rohs
updatedROHs = data.frame() #keep rohs when overlaps are detected
finalROH = data.frame() #permuted rohs
rohsMappedBack = data.frame() #permuted rohs that are mapped back to original chroms
everyonesROHs = list() #keep everyones new rohs

#Load dog info and ROHs
dogInfo = read.delim("~/DogProject_Jaz/LocalRscripts/BreedCladeInfo/breeds_dryad.txt")
rohs = read.delim(file = "~/DogProject_Jaz/LocalRscripts/ROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt", stringsAsFactors = F) %>%
  filter(INDV %in% dogInfo$dogID) #subset to Hayward data only
indivs = unique(rohs$INDV) #grab all individuals

#use the cumulative sum as the end since we are shuffling over the entire length of the genome 
autosome = read.delim(file = "~/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/chromosomeLengths.txt", check.names = F, stringsAsFactors = F, sep = " ") %>%
  mutate(CHROM = as.numeric(gsub("chr", "", CHROM)),
         endUnif = cumsum(as.numeric(LENGTH)),
         START = as.numeric(endUnif) - as.numeric(LENGTH))

#data frame to map everything back to the autosome (scale into kb for Genomic Ranges) after shuffling
mapChrom = autosome %>%
  mutate(start = START/10^3,
         end = endUnif/10^3,
         chrom = "chr1") %>%
  select(start, end, chrom,CHROM)

autoLen = sum(as.numeric(autosome$LENGTH)) #compute the length of the autsome since we are making the genome circular to eliminate edge effects then we will map everything back at the end

#Shuffle the ROH
#Move runs on same chrom by sampling uniformly
for(j in indivs){
  
  #subset rohs to contain only individual of interest
  df = rohs %>%
    filter(INDV == j)
  
  for(i in df$AUTO_LEN){
    rohLength = i
    start = runif(1, min = 1, max = autoLen)
    end = start + rohLength
    #Next check whether rohs fall within a chromosome, if they do not make new roh until it does
    if (end <= autoLen) {
      newROH = cbind(start, end)
      newCoords = rbind(newCoords, newROH)
    }else{
      while (end > autoLen) {
        start = runif(1, min = 1, max = autoLen)
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
        start = runif(1, min = 1, max = autoLen)
        end = start + rohLength
        if (end <= autoLen) {
          retryROH = cbind(start, end)
          updatedROHs = rbind(finalROH, retryROH)
        }else{
          while (end > autoLen) {
            start = runif(1, min = 1, max = autoLen)
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
    #map rohs back to chromosomes
    rangeA = with(finalROH, GRanges(chrom, IRanges(start=start, end =end)))
    rangeB = with(mapChrom, GRanges(chrom, IRanges(start=start, end =end, names= CHROM)))
    type2 = findOverlaps(query = rangeA, subject = rangeB, type = 'any')
    rohsMappedBack = data.frame(finalROH[queryHits(type2),], mapChrom[subjectHits(type2),]) %>%
      select(CHROM, start, end) %>%
      mutate(start = start*10^3,
             end = end*10^3)
  }
  everyonesROHs[[j]] = list(rohsMappedBack)
}
