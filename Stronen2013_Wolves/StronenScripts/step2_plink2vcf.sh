#Load plink
. /u/local/Modules/default/init/modules.sh
module load plink
#Convert to vcf
for f in {1..38}
do 
plink --dog --bfile Stronen2013_CanFam3_MergedFitakCornellSharedSites --recode vcf-iid --chr "$f" --out splitVCF/Stronen2013_wolves_chr"$f"
done
