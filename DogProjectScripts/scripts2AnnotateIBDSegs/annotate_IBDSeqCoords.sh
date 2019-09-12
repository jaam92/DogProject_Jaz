


##qsub -cwd -V -N IBDAnnot -l highp,time=100:00:00,h_data=12G -t 1-38:1 -M eplau -m a annotate_BeagleIBDCoords.sh 


source /u/local/Modules/default/init/modules.sh
module load python 

CHROM=${SGE_TASK_ID}

#Directories
sites_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/splitVCF'
in_dir='/u/flashscratch/j/jmooney3/PhasingDogData/IBDSeq/IBDSegs'
out_dir='/u/flashscratch/j/jmooney3/PhasingDogData/IBDSeq/IBDSegs/Annotate_Segs'




	
	IBDrange_in=${in_dir}/'cornell_canine_chr'$CHROM'_IBDSeqcoords_split09.txt'

	sites_in=${sites_dir}/'cornell_canine_chr'$CHROM'_SitesInVCF.txt'

	count_out=${out_dir}/'Annotated_cornell_canine_chr'$CHROM'_IBDSeqcoords_split09.txt'

	python IBDSeqQualCounts.py ${IBDrange_in} ${sites_in} 1 ${count_out}
	




