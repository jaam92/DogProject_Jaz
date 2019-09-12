#!/bin/bash
#$ -cwd 
#$ -V 
#$ -N S15_IBDSeqIBDNe_allpops
#$ -l highp,time=50:00:00,h_data=10G 
#$ -pe shared 2
#$ -M eplau 
#$ -m bea


#Load Java
. /u/local/Modules/default/init/modules.sh
module load java

for i in {boxer,cocker_spaniel,german_shepherd_dog,golden_retriever,grayWolf,labrador_retriever,maltese,mix,newfoundland,poodle,rottweiler,village_dog_peru,yorkshire_terrier}

do

#Run IBDNe of IBDSeq segments (set the minibd=4cm) This is suggested in the IBDNE paper when using IBDSeq on array based data
cat IBDSeqIBDNe/"$i"_IBDSeqIBDSegs_allchroms_unrelatedsOnly.ibd | java -jar /u/home/j/jmooney3/klohmueldata/jazlyn_data/software/Phasing_IBD/IBDNe/ibdne.04Sep15.e78.jar map=/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/dog_genetic_maps/canFam3.1_average_allChroms_IBDNe.map minibd=4 out=IBDSeqIBDNe/"$i"_IBDNE_usingIBDSeqIBDSegs_Sept15ver

#cat IBDSeqIBDNe/"$i"_IBDSeqIBDSegs_allchroms_unrelatedsOnly.ibd | java -jar /u/home/j/jmooney3/klohmueldata/jazlyn_data/software/Phasing_IBD/IBDNe/ibdne.07May18.6a4.jar map=/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/dog_genetic_maps/canFam3.1_average_allChroms_IBDNe.map mincm=4 out=IBDSeqIBDNe/"$i"_IBDNE_usingIBDSeqIBDSegs_May18ver

done
