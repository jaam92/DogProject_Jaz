#Load vcftools
. /u/local/Modules/default/init/modules.sh
module load vcftools

while read p
do
for f in {1..38} 
do 

#remove the txt extension
outFileName="${p/.txt/}"

#generate ROH for merged data
vcftools --gzvcf MergeFitkakAndCornell/MergedFile_CornellCanineFitak_allIndivs.vcf.gz --keep /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/IndividualFiles/AllIndivsByBreed/$p --LROH --chr "$f" --out vcfToolsROH/splitByBreed/"$outFileName"_chr"$f"

done
done < /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/IndividualFiles/AllIndivsByBreed/sampSize_grEqlto6_fileNames.txt
