#load bedtools
. /u/local/Modules/default/init/modules.sh
module load bedtools

#make file
# grep -v "INDV" /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/vcfToolsROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFitakCornell_sampSizeGrEqlto6.txt | awk '{print "chr"$8"\t"$3-1"\t"$4"\t"$6}' > allROH_withIndivCol_vcfTools.bed

#sorted all the ROHs by the chromosome
bedtools sort -i allROH_withIndivCol_vcfTools.bed > SortedROH_withIndivCol_vcfTools.bed

#intersects exon coordinate with ROHs and outputs the exonic region where there are no ROH overlap
bedtools intersect -a ForAbi_EnsemblGenes_CanFam3.1_SingleTranscript.bed -b SortedROH_withIndivCol_vcfTools.bed -v > ExonRegion_NonOverlapsROH_vcfTools.bed

#count how many exons per gene have no ROH
awk '{print $4}' ExonRegion_NonOverlapsROH_vcfTools.bed | uniq -c | sed -e 1i'CountExon\tGeneNames' > CountExonRegion_NonOverlapsROH_vcfTools.bed

