#!/bin/bash
#$ -l h_data=2G,h_rt=22:00:00
#$ -cwd
#$ -A jmooney3
#$ -N GRM
#$ -e /u/flashscratch/j/jmooney3
#$ -o /u/flashscratch/j/jmooney3
#$ -m a

#Load applications
. /u/local/Modules/default/init/modules.sh
module load xz

##SGE_TASK_ID=942

#Use bedtools to intersect ROH from individuals 
#Count the number of overlaps
#Use sed to insert the information of two individuals being compared without the file paths

while read -r a b 
do 

/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/bedtools2/bin/bedtools coverage -a splitFiles/"$a" -b splitFiles/"$b" | awk '{sum +=6} END {print "S=" sum}' | sed -e "s|^|$a\t$b\t|g" >> Output/'overlaps'$SGE_TASK_ID'.bed' 

#Let us know where we are
echo "$a" "$b" 

done < inputs/'fnames_haywardComps'$SGE_TASK_ID
