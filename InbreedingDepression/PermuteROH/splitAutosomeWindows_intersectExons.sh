
#Load applications
. /u/local/Modules/default/init/modules.sh
module load xz

#split autosome into 100kb windows and then interesect with exons
/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/bedtools2/bin/bedtools makewindows -g AutosomeCoords.bed -w 100000 | /u/home/j/jmooney3/klohmueldata/jazlyn_data/software/bedtools2/bin/bedtools intersect -a ForAbi_EnsemblGenes_CanFam3.1_SingleTranscript.bed -b stdin -wa -wb |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7}'> OverlapExons_and_Autosome_100KbWindows.bed

#map plasmid genome back to previous genome coordinates
/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/bedtools2/bin/bedtools makewindows -g AutosomeCoord_plasmid.bed -w 100000 | /u/home/j/jmooney3/klohmueldata/jazlyn_data/software/bedtools2/bin/bedtools intersect -a stdin -b mapChromsBackFromPlasmid.bed -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$7}' > Autosome_100KbWindows_mapBack.bed 


#map the exon coords to the plasmid like genome
#first paste together the chrom lengths with the plasmid-like genome locations of each chrom
#use awk to add a 0 to file so that chromosome lengths are in bed format and the plamid-like genome location will be treated as an info col. in bedtools 
#intersect the exons with the new genome file we just created
#use awk to map the exon locations to their locations in the plasmid-like genome by adding the plamid-start to the actual start and end of the exons, since the location of exons is based on distance from 0  we treat the plasmid-start as our new 0
paste AutosomeCoords.bed mapChromsBackFromPlasmid.bed | awk '{print $1"\t""0""\t"$2"\t"$6"\t"$4"\t"$5}' | /u/home/j/jmooney3/klohmueldata/jazlyn_data/software/bedtools2/bin/bedtools intersect -a ForAbi_EnsemblGenes_CanFam3.1_SingleTranscript.bed -b stdin -wa -wb | awk '{print "chr1""\t"$2+$9"\t"$3+$9"\t"$4}' > ForAbi_EnsemblGenes_CanFam3.1_SingleTranscript_mapToPlasmid.bed

