##qsub -cwd -V -N dogPhase -l highp,time=20:00:00,h_data=15G -t 1-38:1 -M eplau -m a submit_SHAPEIT_doggos.sh
##qsub -cwd -V -N dogHaplo -l highp,time=00:10:00,h_data=5G -t 1-38:1 -M eplau -m a submit_SHAPEIT_doggos.sh

##directories
SHAPEIT_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/Phasing_IBD/SHAPEIT/bin'
RefFiles_dir='/u/flashscratch/j/jmooney3/PhasingDogData/dog_genetic_maps'
out_dir='/u/flashscratch/j/jmooney3/PhasingDogData/SHAPEIT_Haplotypes'
vcf_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/splitVCF'
CHROM=${SGE_TASK_ID}

##input files
vcfInFile=${vcf_dir}/'cornell_canine_chr'${CHROM}'.vcf'
geneticMapFile=${RefFiles_dir}/'chr'${CHROM}'_average_canFam3.1_rmChromCol.txt'
haplotypeOutFile=${out_dir}/'cornell_canine_chr'${CHROM}'_phasedHaplotypes_shapeIT'
vcfOutFile=${out_dir}/'cornell_canine_chr'${CHROM}'_phasedShapeIT.vcf'

#logfiles
logFile=${out_dir}/'haplotype_Chr'${CHROM}'_dog'

##generate phased data without phasing against a reference
${SHAPEIT_dir}/shapeit -V ${vcfInFile} \
	-M ${geneticMapFile} \
	--output-log ${logFile} \
        -O ${haplotypeOutFile}

###***ENDED UP UP RUNNING THIS PART IN A FOR LOOP INSTEAD OF SUBMITTING A JOB (finished all chroms in about 5mins)######
##convert haplotype to VCF
##${SHAPEIT_dir}/shapeit -convert \
 #       --input-haps ${haplotypeOutFile} \
#	--output-vcf ${vcfOutFile}

sleep 200

