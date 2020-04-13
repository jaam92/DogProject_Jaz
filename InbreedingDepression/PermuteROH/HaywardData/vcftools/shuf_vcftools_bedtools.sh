#!/bin/bash
#$ -l h_data=4G,h_rt=00:20:00
#$ -cwd
#$ -A jmooney3
#$ -N exonPerms
#$ -e /u/flashscratch/j/jmooney3
#$ -o /u/flashscratch/j/jmooney3
#$ -m ea

#Directories 
ref_dir='/u/scratch/j/jmooney3/PermuteROH'
in_dir='/u/scratch/j/jmooney3/PermuteROH/HaywardData/vcftools'
out_dir='/u/scratch/j/jmooney3/PermuteROH/HaywardData/vcftools/CountsPerms'
software_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/bedtools2Beta/bin'

#Load applications
. /u/local/Modules/default/init/modules.sh
module load xz

#Shuffle ROH in each individual and then concatenate in a single file
while read p 
do 

${software_dir}/bedtools shuffle -i $p -g ${ref_dir}/'AutosomeCoord_plasmid.bed' -allowBeyondChromEnd -noOverlapping >> ${out_dir}/'newCoords'$SGE_TASK_ID'.bed'
 
done < ${in_dir}/'fnames.txt' 

##sorted all the ROHs by the chromosome
${software_dir}/bedtools sort -i ${out_dir}/'newCoords'$SGE_TASK_ID'.bed' | ${software_dir}/bedtools intersect -a ${ref_dir}/'ForAbi_EnsemblGenes_CanFam3.1_SingleTranscript_mapToPlasmid.bed' -b stdin -v > ${out_dir}/'ExonRegion_NonOverlapsROH_cornellData'$SGE_TASK_ID'.bed'

##count how many exons and genes have no ROH
wc -l ${out_dir}/'ExonRegion_NonOverlapsROH_cornellData'$SGE_TASK_ID'.bed' > ${out_dir}/'count_nonROH_exons'$SGE_TASK_ID'.txt'
awk '{print $4}' ${out_dir}/'ExonRegion_NonOverlapsROH_cornellData'$SGE_TASK_ID'.bed' | uniq | wc -l > ${out_dir}/'count_nonROH_genes'$SGE_TASK_ID'.txt'

#sorted all the ROHs by the chromosome
#then count how many windows they overlap
${software_dir}/bedtools sort -i ${out_dir}/'newCoords'$SGE_TASK_ID'.bed' | ${software_dir}/bedtools intersect -a ${ref_dir}/'Autosome_100KbWindows_mapBack.bed' -b stdin -c > ${out_dir}/'CountPermutedOverlaps'$SGE_TASK_ID'_100Kb_AutosomalSplits.bed'

rm ${out_dir}/'newCoords'$SGE_TASK_ID'.bed'

sleep 3m
