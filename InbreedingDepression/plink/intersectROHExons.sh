#load bedtools
. /u/local/Modules/default/init/modules.sh
module load bedtools

#sorted all the ROHs by the chromosome
bedtools sort -i allROH_withIndivCol.bed > SortedROH_withIndicCol.bed

#intersects exon coordinate with ROHs and outputs the exonic region where there are no ROH overlap
bedtools intersect -a ../ForAbi_EnsemblGenes_CanFam3.1_SingleTranscript.bed -b SortedROH_withIndivCol.bed -v > ExonRegion_NonOverlapsROH.bed

#count how many exons per gene have no ROH
awk '{print $4}' ExonRegion_NonOverlapsROH.bed | uniq -c | sed -e 1i'CountExon\tGeneNames' > CountExonRegion_NonOverlapsROH.bed

