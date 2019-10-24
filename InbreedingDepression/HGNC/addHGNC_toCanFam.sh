#Malika Kumar, updated Jan 24 2018. This script takes in a raw genelist (one column, list of gene names) and (1) converts all names to HGNC approved symbols, (2) identifies coordinates for each gene, (3) saves coordinates and new gene list into two outfiles. Two intermediate tempfiles are created, and removed at the end of the script.
#usage: ./format_genelist.sh <input_file>



#Load bedtools
. /u/local/Modules/default/init/modules.sh
module load bedtools

INFILE=$1 #you could also change this to a path, and avoid having to pass the script a parameter

#output directory & files
OUT_DIR=/u/scratch/j/jmooney3/ForAbi #change these to your own!
COORD_FILE=${OUT_DIR}/EnsemblGenes_CanFam3.1_SingleTranscript_hg19_HGNCrename.bed
GENES_FILE=${OUT_DIR}/EnsemblGenes_CanFam3.1_SingleTranscript_hg19_HGNCrename.txt

#reference directory & files
JAZLYN_DIR=/u/scratch/j/jmooney3/ForAbi
SYM_REF=${JAZLYN_DIR}/hgnc_symbol_pairs.txt
COORD_REF=${JAZLYN_DIR}/EnsemblGenes_CanFam3.1_SingleTranscript_v2.bed

#Join list of input genes with HGNC approved, previous, and synonym symbols to identify HGNC approved symbols for each gene. Then identify chromosomal location from protein_coding_genes.bed (COORD_REF).
join -1 1 -2 1 <(sort -k1,1 $INFILE) <(sort -k1,1 $SYM_REF) | awk '{print $2}' | sort -u | join -1 1 -2 4 - <(sort -k4,4 $COORD_REF) | awk '{print $2, $3, $4}' OFS="\t" > ${COORD_FILE}.temp

#Pull all entries from protein_coding_genes.bed that match chromosomal location with these genes to catch any genes that share chromosomal location but whose symbols hadn't been linked.
while read p; do grep -w "${p}" $COORD_REF >> ${COORD_FILE}.temp2 ; done < ${COORD_FILE}.temp

#Sort to remove duplicates
sort -u ${COORD_FILE}.temp2 | bedtools sort -i > $COORD_FILE

#Extract gene names
#awk '{print $4}' $COORD_FILE | sort -u > $GENES_FILE #comment out this line if you don't need a file with the list of gene names (one column)

#Remove intermediate files
rm ${COORD_FILE}.*    
