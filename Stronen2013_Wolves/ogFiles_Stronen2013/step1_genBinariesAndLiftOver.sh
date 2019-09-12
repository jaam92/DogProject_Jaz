###Download Files###
#wget https://datadryad.org/bitstream/handle/10255/dryad.174676/WolvesAll_22feb2012.ped?sequence=2 
#wget https://datadryad.org/bitstream/handle/10255/dryad.174677/WolvesAll_22feb2012.map?sequence=2
#wget https://datadryad.org/bitstream/handle/10255/dryad.174682/Ukrainian_Steppe.txt?sequence=2
#wget https://datadryad.org/bitstream/handle/10255/dryad.174681/Northcentral_Europe.txt?sequence=2
#wget https://datadryad.org/bitstream/handle/10255/dryad.174680/Carpathians.txt?sequence=2
#wget https://datadryad.org/bitstream/handle/10255/dryad.174679/Dinaric-Balkan.txt?sequence=2
#wget https://datadryad.org/bitstream/handle/10255/dryad.174678/Italy.txt?sequence=2

###Generate Binaries######
#Get the list of filtered SNPs used in the Stronen et al. PLoS One 2013 paper
#all the .hwe files have the same amout of SNPs so pick one to make the snp list  

#Command to check snps
#cat *.hwe | awk '{print $2}' | sort | uniq -c | head

#Command to output snplist
#awk '{print $2}' HWE_files/Carpathians.hwe | grep -v "SNP" > SNPList_Filtered_Stronen2013.txt

#Turn ped/map files into binary and only keep autosomes
#. /u/local/Modules/default/init/modules.sh
#module load plink
#plink --dog --ped WolvesAll_22feb2012.ped --map WolvesAll_22feb2012.map --chr 1-38 --extract SNPList_Filtered_Stronen2013.txt --remove RmIndivsGrThan10perMissing.txt --make-bed --out WolvesAll_22feb2012

######LIFT OVER#########
#Create bedfile for ucsc to switch from canfam 2 to canfam 3
#awk '{print $1, $4-1, $4, $2}' WolvesAll_22feb2012.bim | sed 's/^/chr/' > LiftOver/Wolves_LiftOver_ConvertCanFam2toCanFam3.txt

#Create snp list for recoding and excluding
#awk '{print $4"\t"$3}' LiftOver/hglft_genome_ConvertCanFam2toCanFam3.bed > LiftOver/recodeSNPs_CanFam3.txt
#grep -v "#" LiftOver/SitesThatFailed_n114.txt | awk '{print $4}' > LiftOver/removeSites_CanFam3.txt

#Create new plink file
#. /u/local/Modules/default/init/modules.sh
#module load plink
#plink --bfile WolvesAll_22feb2012 --dog --update-map LiftOver/recodeSNPs_CanFam3.txt --exclude LiftOver/removeSites_CanFam3.txt --make-bed  --out LiftOver/IntermidiateFile/WolvesAll_22feb2012_LiftOvercanFam3


#####Subset down to sites shared by Merged Fitak and Hayward data#####
#awk '{print $2}' ../../MergeFitkakAndCornell/MergedFile_CornellCanineFitak_UnrelatedsOnly.bim > LiftOver/IntermidiateFile/snpList_MergedFitakCornell.txt
#plink --bfile LiftOver/IntermidiateFile/WolvesAll_22feb2012_LiftOvercanFam3 --dog --extract LiftOver/IntermidiateFile/snpList_MergedFitakCornell.txt  --make-bed  --out ../Stronen2013_CanFam3_MergedFitakCornellSharedSites
