#. /u/local/Modules/default/init/modules.sh
#module load plink

for f in {1..38} 
do 

#do not gzip vcf otherwise vcftools LROH won't work
plink --dog --bfile MergeFitkakAndCornell/MergedFile_CornellCanineFitak_allIndivs --recode vcf-fid --chr "$f" --out splitVCF/MergedFile_CornellCanineFitak_allIndivs_chr"$f"

plink --dog --bfile MergeFitkakAndCornell/MergedFile_CornellCanineFitak_UnrelatedsOnly --recode vcf-fid --chr "$f" --out splitVCF/MergedFile_CornellCanineFitak_UnrelatedsOnly_chr"$f"

done
