


##qsub -cwd -V -N IBDAnnot -l highp,time=00:45:00,h_data=12G -t 1-38:1 -M eplau -m ea annotate_BeagleIBDCoords.sh 


source /u/local/Modules/default/init/modules.sh
module load python 

CHROM=${SGE_TASK_ID}

#Directories
sites_dir='/u/flashscratch/j/jmooney3/PhasingDogData/SHAPEIT_PhasedVCF'
in_dir='/u/flashscratch/j/jmooney3/PhasingDogData/BeagleIBD'
out_dir='/u/flashscratch/j/jmooney3/PhasingDogData/BeagleIBD/Annotate_Segs'




	
	IBDrange_in=${in_dir}/'cornell_canine_chr'$CHROM'_BeagleCoords_sharedIBD.txt'

	sites_in=${sites_dir}/'cornell_canine_chr'$CHROM'_SitesInVCF.txt'

	count_out=${out_dir}/'Annotated_cornell_canine_chr'$CHROM'_BeagleCoords_sharedIBD_rd2.txt'

	python BeagleQualCounts.py ${IBDrange_in} ${sites_in} 1 ${count_out}
	


sleep 400

