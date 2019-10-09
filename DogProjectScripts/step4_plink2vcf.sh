#. /u/local/Modules/default/init/modules.sh
#module load plink

for f in {1..38} 
do 

#make vcf the gzip
plink --dog --bfile MergeFitkakAndCornell/MergedFile_CornellCanineFitak_allIndivs --recode vcf-fid --chr "$f" --out splitVCF/MergedFile_CornellCanineFitak_allIndivs_chr"$f"

bgzip -c splitVCF/MergedFile_CornellCanineFitak_allIndivs_chr"$f".vcf > splitVCF/MergedFile_CornellCanineFitak_allIndivs_chr"$f".vcf.gz
tabix -p vcf splitVCF/MergedFile_CornellCanineFitak_allIndivs_chr"$f".vcf.gz

plink --dog --bfile MergeFitkakAndCornell/MergedFile_CornellCanineFitak_UnrelatedsOnly --recode vcf-fid --chr "$f" --out splitVCF/MergedFile_CornellCanineFitak_UnrelatedsOnly_chr"$f"

bgzip -c splitVCF/MergedFile_CornellCanineFitak_UnrelatedsOnly_chr"$f".vcf > splitVCF/MergedFile_CornellCanineFitak_UnrelatedsOnly_chr"$f".vcf.gz 
tabix -p vcf splitVCF/MergedFile_CornellCanineFitak_UnrelatedsOnly_chr"$f".vcf.gz

#delete unzipped vcf 
rm -f splitVCF/MergedFile_CornellCanineFitak_allIndivs_chr"$f".vcf
rm -f splitVCF/MergedFile_CornellCanineFitak_UnrelatedsOnly_chr"$f".vcf

done
