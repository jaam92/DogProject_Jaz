#!/bin/bash
#$ -l highp,h_data=28G,h_rt=03:20:00,h_vmem=56G 
#$ -pe shared 2
#$ -cwd
#$ -A jmooney3
#$ -N set2VROHperms
#$ -e /u/flashscratch/j/jmooney3
#$ -o /u/flashscratch/j/jmooney3
#$ -m a

#Directories 
ref_dir='/u/scratch/j/jmooney3/PermuteROH'
in_dir='/u/scratch/j/jmooney3/PermuteROH/HaywardData/vcftools'
out_dir='/u/scratch/j/jmooney3/PermuteROH/HaywardData/vcftools/CountsPerms'
software_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/bedtools2/bin'

#Load applications
. /u/local/Modules/default/init/modules.sh
module load xz
module load R/3.4.2

#use R to shuffle roh
Rscript --vanilla shuffROH_final.R --SGETaskID $SGE_TASK_ID --genomeFile ${ref_dir}/'chromosomeLengths.txt' --rohInfile /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/InbreedingDepression/HaywardData/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFile_Cornell_allChroms_vcfToolsROH_rmROHlessThan50snps.txt --outFilePath /u/scratch/j/jmooney3/PermuteROH/HaywardData/vcftools/CountsPerms/

#use R to map shuffled ROH back to locations on genome
Rscript --vanilla mapCoordsBack.R --SGETaskID $SGE_TASK_ID --genomeFile ${ref_dir}/'chromosomeLengths.txt' --rohInfile ${out_dir}/'newCoords'$SGE_TASK_ID'.bed' --outFilePath /u/scratch/j/jmooney3/PermuteROH/HaywardData/vcftools/CountsPerms/


#sorted all the ROHs by the chromosome
${software_dir}/bedtools sort -i ${out_dir}/'newCoords'$SGE_TASK_ID'_mappedBack.bed' | ${software_dir}/bedtools intersect -a ${ref_dir}/'ForAbi_EnsemblGenes_CanFam3.1_SingleTranscript.bed' -b stdin -v > ${out_dir}/'ExonRegion_NonOverlapsROH_cornellData'$SGE_TASK_ID'.bed'

#count how many exons and genes have no ROH
wc -l ${out_dir}/'ExonRegion_NonOverlapsROH_cornellData'$SGE_TASK_ID'.bed' > ${out_dir}/'count_nonROH_exons'$SGE_TASK_ID'.txt'
awk '{print $4}' ${out_dir}/'ExonRegion_NonOverlapsROH_cornellData'$SGE_TASK_ID'.bed' | uniq | wc -l > ${out_dir}/'count_nonROH_genes'$SGE_TASK_ID'.txt'

#delete permuted rohs
##rm ${out_dir}/'newCoords'$SGE_TASK_ID'.bed'

#sleep
sleep 2m
