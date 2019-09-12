##qsub -cwd -V -N IBDtracts -l highp,time=02:30:00,h_data=10G -t 1-38:1 -M eplau -m a step6_runIBDSeq.sh



###chromosome 1 is going to need 5 hours###
#Directories
in_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/SubsetStronen2013Sites/VCFs'
out_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/SubsetStronen2013Sites/IBDSeq'
ibdSeq_dir='/u/home/j/jmooney3/klohmueldata/jazlyn_data/software/Phasing_IBD/IBDSeq'

CHROM=${SGE_TASK_ID}

#Files
vcfin=${in_dir}/'MergedFile_CornellCanineFitak_UnrelatedIndivs_Stronen2013Sites_chr'${CHROM}'.vcf'
ibdout=${out_dir}/'MergedFile_CornellCanineFitak_UnrelatedIndivs_Stronen2013Sites_chr'${CHROM}'_phasedHaplotypes_IBDSeq'

#Load Java Module
. /u/local/Modules/default/init/modules.sh
module load java/1.8.0_77

#IBDSeq command
java -Xmx4G -jar ${ibdSeq_dir}/ibdseq.r1206.jar gt=${vcfin} out=${ibdout}

sleep 200

