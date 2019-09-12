##qsub -cwd -V -N IBDtracts -l highp,time=02:00:00,h_data=15G -t 1-38:1 -M eplau -m bea submit_BeagleIBD.sh


#Directories
in_dir='/u/flashscratch/j/jmooney3/PhasingDogData/SHAPEIT_PhasedVCF'
out_dir='/u/flashscratch/j/jmooney3/PhasingDogData/BeagleIBD'
beagle_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/Phasing_IBD/Beagle_v4.1'
recomMap_dir='/u/flashscratch/j/jmooney3/PhasingDogData/dog_genetic_maps'
CHROM=${SGE_TASK_ID}

#Files
recomMap_in=${recomMap_dir}/'chr'${CHROM}'_average_canFam3.1_noHeader_beagleFormat.map'

vcfin=${in_dir}/'cornell_canine_chr'${CHROM}'_phasedHaplotypes_shapeIT.vcf'
ibdout=${out_dir}/'cornell_canine_chr'${CHROM}'_phasedHaplotypes_shapeIT'

#Load Java Module
. /u/local/Modules/default/init/modules.sh
module load java/1.8.0_77

#Command
java -Xss5m -Xmx4G -jar ${beagle_dir}/beagle.08Jun17.d8b.jar gt=${vcfin} impute=false ibd=true gprobs=false lowmem=true overlap=25 ibdtrim=10 out=${ibdout} map=${recomMap_in}

sleep 200
