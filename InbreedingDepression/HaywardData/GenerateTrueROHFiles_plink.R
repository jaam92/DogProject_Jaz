library(dplyr)
library(data.table)

#Read ROH Files in 
df = read.delim("/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/InbreedingDepression/HaywardData/plink/ROH_CornellCanine_plink_rmLt50Rm.txt")

#Look at proportion covered
df$PropCovered = df$NSNP/df$KB

##Identify ROH greater than 100Kb and remove ROH where the avg. prop covered by SNPs is within 1 SD of the mean 
ROHgr100kb = df[which(df$KB >= 100),]
z = data.table(ROHgr100kb)
z[,ToKeep := abs(ROHgr100kb$PropCovered - mean(ROHgr100kb$PropCovered)) < sd(ROHgr100kb$PropCovered)][ToKeep  == TRUE] #ID ROH to keep
FinalDF_allROH = subset(z, z$ToKeep == "TRUE") #Subset out true ROH

#Write final ROH to outfile
write.table(x = FinalDF_allROH,file = "/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/InbreedingDepression/HaywardData/plink/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFile_Cornell_allChroms_plinkROH_rmROHlessThan50snps.txt" , sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

