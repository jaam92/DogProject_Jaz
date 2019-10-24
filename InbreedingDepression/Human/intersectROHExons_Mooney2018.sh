#Make ROH bed file
#awk '{print "chr"$7"\t"$1-1"\t"$2}' TrueROH_Unrelated_perIndivFINAL_Aug17_centromereOverlapRm.txt | grep -v "AUTO_END" > TrueROH_Unrelated_perIndivFINAL_Aug17_centromereOverlapRm.bed

#Add gene names to exons
#bedtools intersect -wa -wb  -a protein_coding_genes_hg19_HGNC.bed -b exome.canonical.strict.bed | awk '{print $1"\t"$7"\t"$8"\t"$4}' > exome.canonical.strict.HGNCGeneNames.bed 

#Intersect ROH with exons
#bedtools sort -i TrueROH_Unrelated_perIndivFINAL_Aug17_centromereOverlapRm.bed | bedtools intersect -a exome.canonical.strict.HGNCGeneNames.bed -b stdin -v > ExonRegion_NonOverlapsROH_Mooney.bed

##count how many exons per gene have no ROH
#awk '{print $4}' ExonRegion_NonOverlapsROH_Mooney.bed | uniq -c | sed -e 1i'CountExon\tGeneNames' > CountExonRegion_NonOverlapsROH_MooneyData.bed
