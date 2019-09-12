#Create a vcf with all chromosomes and all indivudals and annotate with snpEff

#Generate vcf on scratch
#Load plink
#. /u/local/Modules/default/init/modules.sh
#module load plink
#make vcf
#plink --dog --bfile MergeFitkakAndCornell/MergedFile_CornellCanineFitak_allIndivs --recode vcf-fid --out $SCRATCH/MergedFile_CornellCanineFitak_allIndivs_allChroms



#Annotate vcf
#. /u/local/Modules/default/init/modules.sh
#module load java
#java -Xmx4g -jar /u/home/j/jmooney3/klohmueldata/jazlyn_data/software/snpEff/snpEff.jar CanFam3.1.86 /u/flashscratch/j/jmooney3/MergedFile_CornellCanineFitak_allIndivs_allChroms.vcf > /u/flashscratch/j/jmooney3/MergedFile_CornellCanineFitak_allIndivs_allChroms_ANNOTATEDsnpEFF.vcf
