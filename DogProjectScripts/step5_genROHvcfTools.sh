#Load vcftools
. /u/local/Modules/default/init/modules.sh
module load vcftools

for f in {1..38} 
do 

#generate ROH for merged data
#all individuals
#vcftools --gzvcf splitVCF/MergedFile_CornellCanineFitak_allIndivs_chr"$f".vcf.gz --LROH --chr "$f" --out vcfToolsROH/allIndivs/mergedFile_chr"$f"

#at least 6 individuals per group
vcftools --gzvcf splitVCF/MergedFile_CornellCanineFitak_allIndivs_chr"$f".vcf.gz --keep /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/IndividualFiles/AllIndivsByBreed/sampSize_grEqlto6_indivs.txt --LROH --chr "$f" --out vcfToolsROH/mergedFile_chr"$f"_sampSizeGrEqlto6

done
