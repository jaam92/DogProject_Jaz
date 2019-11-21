#Load libraries
library(GenometriCorr)

#Load files
rohs = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFile_Cornell_allChroms_vcfToolsROH_rmROHlessThan50snps.bed", stringsAsFactors = F)
exons = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/ForAbi_EnsemblGenes_CanFam3.1_SingleTranscript.bed", stringsAsFactors = F)
chromInfo = read.table("~/Documents/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/chromosomeLengths.txt",  header = T, stringsAsFactors = F)

#Make genomic ranges objects
rohs_query = with(rohs, GRanges(V1, IRanges(start=V2, end = V3)))
exons_ref = with(exons, GRanges(V1, IRanges(start=V2, end = V3)))
chromLengths = as.numeric(chromInfo$LENGTH)
names(chromLengths) = as.character(chromInfo$CHROM)   

#Check overlaps
pn.area = 10
pn.dist = 10
pn.jacc = 10

overlaps_exons = GenometriCorrelation(rohs_query, 
                                     exons_ref,
                                     ecdf.area.permut.number=pn.area,
                                     mean.distance.permut.number=pn.dist,
                                     jaccard.measure.permut.number=pn.jacc,
                                     keep.distributions=FALSE,
                                     showProgressBar=TRUE)