#split files (need to remove the slash in N/A)
sed -e 's/N\/A/NA/g' TrueROH_propCoveredwithin1SDMean_allChroms_mergedFile_Cornell_allChroms_vcfToolsROH_rmROHlessThan50snps.bed | awk '{print >> ($4".bed")}'

#list of file names with path
ls -1 splitFiles/*.bed > fnames.txt
