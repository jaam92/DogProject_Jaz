#. /u/local/Modules/default/init/modules.sh
#module load plink

#Subset merged files to stronen sites
#plink --bfile MergeFitkakAndCornell/MergedFile_CornellCanineFitak_UnrelatedsOnly --dog --extract SubsetStronen2013Sites/Stronen2013Sites.txt --make-bed --out SubsetStronen2013Sites/MergedFile_CornellCanineFitak_UnrelatedIndivs_Stronen2013Sites
#plink --bfile MergeFitkakAndCornell/MergedFile_CornellCanineFitak_allIndivs --dog --extract SubsetStronen2013Sites/Stronen2013Sites.txt --make-bed --out SubsetStronen2013Sites/MergedFile_CornellCanineFitak_allIndivs_Stronen2013Sites

#Make vcfs
#for f in {1..38} 
#do 

#do not gzip vcf otherwise vcftools LROH won't work
#plink --dog --bfile SubsetStronen2013Sites/MergedFile_CornellCanineFitak_allIndivs_Stronen2013Sites --recode vcf-fid --chr "$f" --out SubsetStronen2013Sites/VCFs/MergedFile_CornellCanineFitak_allIndivs_chr"$f"
#plink --dog --bfile SubsetStronen2013Sites/MergedFile_CornellCanineFitak_UnrelatedIndivs_Stronen2013Sites --recode vcf-fid --chr "$f" --out SubsetStronen2013Sites/VCFs/MergedFile_CornellCanineFitak_UnrelatedIndivs_Stronen2013Sites_chr"$f"

#done

#Call ROHs
. /u/local/Modules/default/init/modules.sh
module load vcftools 
for f in {1..38} 
do
vcftools --vcf SubsetStronen2013Sites/VCFs/MergedFile_CornellCanineFitak_allIndivs_chr"$f".vcf --chr "$f" --LROH --out SubsetStronen2013Sites/ROH/mergedFile_CornellCanineFitak_StronenSites_chr"$f"
done 
