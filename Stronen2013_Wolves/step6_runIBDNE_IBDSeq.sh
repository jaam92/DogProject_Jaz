#!/bin/bash
#$ -cwd 
#$ -V 
#$ -N Wolves_IBDNe_allpops
#$ -l highp,time=20:00:00,h_data=10G 
#$ -pe shared 2
#$ -M eplau 
#$ -m bea


#Load Java
. /u/local/Modules/default/init/modules.sh
module load java

for i in {Italy,Northcentral-Europe}
do


#Run IBDNE on IBDSeq segments (add mincm=4) this comes from IBDNe paper where they call segments on array data using IBDSeq
cat /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/Stronen2013_Wolves/IBDNE/"$i"_IBDSeqIBDSegs_allchroms_unrelatedsOnly.ibd | java -jar /u/home/j/jmooney3/klohmueldata/jazlyn_data/software/Phasing_IBD/IBDNe/ibdne.04Sep15.e78.jar map=/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/dog_genetic_maps/canFam3.1_average_allChroms_IBDNe.map minibd=4 out=/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/Stronen2013_Wolves/IBDNE/"$i"_IBDNE_usingIBDSeqSegs_Sep15ver 

cat /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/Stronen2013_Wolves/IBDNE/"$i"_IBDSeqIBDSegs_allchroms_unrelatedsOnly.ibd | java -jar /u/home/j/jmooney3/klohmueldata/jazlyn_data/software/Phasing_IBD/IBDNe/ibdne.07May18.6a4.jar map=/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/dog_genetic_maps/canFam3.1_average_allChroms_IBDNe.map mincm=4 out=/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/Stronen2013_Wolves/IBDNE/"$i"_IBDNE_usingIBDSeqSegs_May18ver

done
