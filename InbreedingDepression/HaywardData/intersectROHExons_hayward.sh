#!/bin/bash
#$ -cwd
#$ -V
#$ -N callROH
#$ -l h_data=1G,time=02:00:00
#$ -M eplau
#$ -m bea

#Load applications
. /u/local/Modules/default/init/modules.sh
module load plink
module load vcftools
module load R/3.4.2
module load bedtools

#make vcf
#plink --dog --bfile /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/ogFiles_Hayward2016/cornell_canine --recode vcf-fid --out cornell_canine

#for f in {1..38} 
#do 
#call roh
#vcftools --vcf cornell_canine.vcf --LROH --chr "$f" --out cornell_canine_chr"$f"
#done

#Filter ROH
#cat cornell_canine_chr*.LROH | awk '$6 > 49 {print $0}'| grep -v "CHROM" | sed 1i'CHROM\tAUTO_START\tAUTO_END\tMIN_START\tMAX_END\tN_VARIANTS_BETWEEN_MAX_BOUNDARIES\tN_MISMATCHES\tINDV' > mergedFile_Cornell_allChroms_vcfToolsROH_rmROHlessThan50snps.txt

#generate true ROH
#Rscript GenerateTrueROHFiles.R

#convert to bed format
awk '{print "chr"$8"\t"$3-1"\t"$4"\t"$6}' TrueROH_propCoveredwithin1SDMean_allChroms_mergedFile_Cornell_allChroms_vcfToolsROH_rmROHlessThan50snps.txt | grep -v "CHROM" > TrueROH_propCoveredwithin1SDMean_allChroms_mergedFile_Cornell_allChroms_vcfToolsROH_rmROHlessThan50snps.bed

#sorted all the ROHs by the chromosome
bedtools sort -i TrueROH_propCoveredwithin1SDMean_allChroms_mergedFile_Cornell_allChroms_vcfToolsROH_rmROHlessThan50snps.bed | bedtools intersect -a ForAbi_EnsemblGenes_CanFam3.1_SingleTranscript.bed -b stdin -v > ExonRegion_NonOverlapsROH_cornellData.bed

#count how many exons per gene have no ROH
awk '{print $4}' ExonRegion_NonOverlapsROH_cornellData.bed | uniq -c | sed -e 1i'CountExon\tGeneNames' > CountExonRegion_NonOverlapsROH_cornellData.bed

#overlaps between all three
cat CountExonRegion_NonOverlapsROH_cornellData.bed ../plink/CountExonRegion_NonOverlapsROH.bed ../vcftools/CountExonRegion_NonOverlapsROH_vcfTools.bed | grep -v "Count" | awk '{print $2}' | sort | uniq -c | awk '$1 ==3 {print}'
