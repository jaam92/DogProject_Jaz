
#Load applications
. /u/local/Modules/default/init/modules.sh
module load xz

#split autosome into 100kb windows and then interesect with exons

/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/bedtools2/bin/bedtools makewindows -g AutosomeCoords.bed -w 100000 | /u/home/j/jmooney3/klohmueldata/jazlyn_data/software/bedtools2/bin/bedtools intersect -a ForAbi_EnsemblGenes_CanFam3.1_SingleTranscript.bed -b stdin -wa -wb |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7}'> OverlapExons_and_Autosome_100KbWindows.bed

#map plasmid genome back to previous genome coordinates
/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/bedtools2/bin/bedtools makewindows -g AutosomeCoord_plasmid.bed -w 100000 | /u/home/j/jmooney3/klohmueldata/jazlyn_data/software/bedtools2/bin/bedtools intersect -a stdin -b mapChromsBackFromPlasmid.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$7}' > Autosome_100KbWindows_mapBack.bed 

