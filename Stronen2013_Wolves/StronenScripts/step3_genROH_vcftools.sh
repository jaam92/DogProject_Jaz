#Load vcftools
. /u/local/Modules/default/init/modules.sh
module load vcftools

#Run command locally (ran in about 1 min)
for f in {1..38} 
do 
vcftools --vcf splitVCF/Stronen2013_wolves_chr"$f".vcf --LROH --chr "$f" --out vcfToolsROH/Stronen2013_wolves_chr"$f"
done
