#load bedtools
. /u/local/Modules/default/init/modules.sh
module load bedtools

#Intersect ROH and GERP scores and keep only non-NA gerp scores
bedtools intersect -a ../masterFileCanFam3_withGerpScore.bed -b SortedROH_withIndivCol.bed -wa -wb | awk '$4 != "NA" {print $0}' > allNonZeroGERPScore_ROHOverlaps.bed 
