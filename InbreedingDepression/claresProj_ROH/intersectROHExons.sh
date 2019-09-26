#grab data without IR
#grep -v "IR" /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/FinalROHQCd_min10Kb_allIndivs_Dec2018_DogProjClare.LROH | awk '{print $1"\t"$2-1"\t"$3"\t"$11}' | grep -v "CHROM" > FinalROHQCd_min10Kb_allIndivs_Dec2018_DogProjClare_rmIR.bed

#grab data without IR
#grep -v "IR" /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Clare/FinalROHQCd_min1Mb_allIndivs_Dec2018_DogProjClare.LROH | awk '{print $1"\t"$2-1"\t"$3"\t"$11}' | grep -v "CHROM" > FinalROHQCd_min1Mb_allIndivs_Dec2018_DogProjClare_rmIR.bed

#load bedtools
. /u/local/Modules/default/init/modules.sh
module load bedtools

#sorted all the ROHs by the chromosome
bedtools sort -i FinalROHQCd_min10Kb_allIndivs_Dec2018_DogProjClare_rmIR.bed | bedtools intersect -a ForAbi_EnsemblGenes_CanFam3.1_SingleTranscript.bed -b stdin -v > ExonRegion_NonOverlapsROH_seqData.bed

#count how many exons per gene have no ROH
awk '{print $4}' ExonRegion_NonOverlapsROH_seqData.bed | uniq -c | sed -e 1i'CountExon\tGeneNames' > CountExonRegion_NonOverlapsROH_seqData.bed

#overlaps between all three
cat CountExonRegion_NonOverlapsROH.bed ../plink/CountExonRegion_NonOverlapsROH.bed ../vcftools/CountExonRegion_NonOverlapsROH_vcfTools.bed | grep -v "Count" | awk '{print $2}' | sort | uniq -c | awk '$1 ==3 {print}'


##### 1MB min for ROH #####
#sorted all the ROHs by the chromosome
bedtools sort -i FinalROHQCd_min1Mb_allIndivs_Dec2018_DogProjClare_rmIR.bed | bedtools intersect -a ForAbi_EnsemblGenes_CanFam3.1_SingleTranscript.bed -b stdin -v > ExonRegion_NonOverlapsROH_seqData_1Mb.bed

#count how many exons per gene have no ROH
awk '{print $4}' ExonRegion_NonOverlapsROH_seqData_1Mb.bed | uniq -c | sed -e 1i'CountExon\tGeneNames' > CountExonRegion_NonOverlapsROH_seqData_1Mb.bed

#overlaps between all three
cat CountExonRegion_NonOverlapsROH_1Mb.bed ../plink/CountExonRegion_NonOverlapsROH.bed ../vcftools/CountExonRegion_NonOverlapsROH_vcfTools.bed | grep -v "Count" | awk '{print $2}' | sort | uniq -c | awk '$1 ==3 {print}'
