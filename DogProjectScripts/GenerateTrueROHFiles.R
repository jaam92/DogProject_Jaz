library(dplyr)
library(data.table)

#Read ROH Files in 
df = read.delim("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/SubsetStronen2013Sites/ROH/StronenSitesOnly_mergedFile_FitakCornell_allChroms_vcfToolsROH_rmROHlessThan50snps.txt")
#df = read.delim("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/vcfToolsROH/mergedFile_FitakCornell_allChroms_vcfToolsROH_rmROHlessThan50snps.txt")
#df = read.delim("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/Stronen2013_Wolves/vcfToolsROH/Stronen2013_Wolves_Stronen2013_allChroms_vcfToolsROH_rmROHlessThan50snps.txt")

##Length of autozygous segment
df$AUTO_END = as.numeric(as.character(df$AUTO_END))
df$AUTO_START = as.numeric(as.character(df$AUTO_START))
df$AUTO_LEN = df$AUTO_END - df$AUTO_START
df$PropCovered = df$N_VARIANTS_BETWEEN_MAX_BOUNDARIES/df$AUTO_LEN

##Identify ROH greater than 100Kb and remove ROH where the avg. prop covered by SNPs is within 1 SD of the mean 
ROHgr100kb = df[which(df$AUTO_LEN >= 100000),][c("MIN_START", "MAX_END","AUTO_START", "AUTO_END", "AUTO_LEN", "INDV", "PropCovered", "CHROM")]
z = data.table(ROHgr100kb)
z[,ToKeep := abs(ROHgr100kb$PropCovered - mean(ROHgr100kb$PropCovered)) < sd(ROHgr100kb$PropCovered)][ToKeep  == TRUE] #ID ROH to keep
FinalDF_allROH = subset(z, z$ToKeep == "TRUE") #Subset out true ROH

#Write final ROH to outfile
write.table(x = FinalDF_allROH,file = "Stronen2013SitesOnlyTrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(x = FinalDF_allROH,file = "TrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(x = FinalDF_allROH,file = "Stronen2013_WolvesTrueROH_propCoveredwithin1SDMean_allChroms.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
