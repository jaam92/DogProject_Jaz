#Ran locally finished in about 5 minutes

#Load Java Module
. /u/local/Modules/default/init/modules.sh
module load java/1.8.0_77

#IBDSeq command
for f in {1..38}; do java -Xmx4G -jar /u/home/j/jmooney3/klohmueldata/jazlyn_data/software/Phasing_IBD/IBDSeq/ibdseq.r1206.jar gt=/u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/Stronen2013_Wolves/splitVCF/Stronen2013_wolves_chr"$f".vcf out=IBDSeq/Stronen2013_wolves_chr"$f"_Haplotypes_IBDSeq;done
