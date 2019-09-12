#Load vcftools
. /u/local/Modules/default/init/modules.sh
module load vcftools

for f in {1..38} 
do 

#generate ROH for merged data
vcftools --vcf splitVCF/MergedFile_CornellCanineFitak_allIndivs_chr"$f".vcf --LROH --chr "$f" --out vcfToolsROH/mergedFile_chr"$f"

done
