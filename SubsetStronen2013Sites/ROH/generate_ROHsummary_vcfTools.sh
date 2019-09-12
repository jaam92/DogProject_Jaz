
#Filter ROH
#awk '$6 > 49 {print $0}' mergedFile_*.LROH | grep -v "CHROM" | sed 1i'CHROM\tAUTO_START\tAUTO_END\tMIN_START\tMAX_END\tN_VARIANTS_BETWEEN_MAX_BOUNDARIES\tN_MISMATCHES\tINDV' > StronenSitesOnly_mergedFile_FitakCornell_allChroms_vcfToolsROH_rmROHlessThan50snps.txt

#Count of ROH per indiv
#awk '{print $8}' StronenSitesOnly_mergedFile_FitakCornell_allChroms_vcfToolsROH_rmROHlessThan50snps.txt | sort | uniq -c | sed 1i'countROH\tID' > countROHperIndiv.txt

. /u/local/Modules/default/init/modules.sh
module load R/3.4.2
Rscript GenerateTrueROHFiles.R 
