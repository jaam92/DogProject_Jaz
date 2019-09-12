#####Download data####
#wget https://datadryad.org/bitstream/handle/10255/dryad.171257/MERGED.clean.pruned.ped?sequence=1 
#wget https://datadryad.org/bitstream/handle/10255/dryad.171258/MERGED.clean.pruned.map?sequence=1

#Rename
#mv MERGED.clean.pruned.map?sequence=1 MERGED.clean.Fitak2018.map
#mv MERGED.clean.pruned.ped?sequence=1 MERGED.clean.Fitak2018.ped

#####THESE COMMANDS MAKE THE FINAL DATA FILE FOR THE WOLVES FROM FITAK 2018######
#. /u/local/Modules/default/init/modules.sh
#module load plink

#Convert to binary format
#plink  --dog --ped MERGED.clean.Fitak2018.ped --map MERGED.clean.Fitak2018.map --chr 1-38 --make-bed --out MERGED.clean.Fitak2018

#Lift over from canfam 2 to canfam 3
#awk '{print $1, $4-1, $4, $2}' MERGED.clean.Fitak2018.bim | sed 's/^/chr/' > Fitak2018_LiftOver_ConvertCanFam2toCanFam3.txt

#Sites to keep
#awk '{print $4"\t"$3}' LiftOver/hglft_genome_FitakLiftOverToCanFam3.bed   > LiftOver/recodeSNPs_Fitak_CanFam3.txt

#Sites to remove
#grep -v "#" LiftOver/SitesThatFailed_n210.txt | awk '{print $4}' > LiftOver/removeSites_Fitak_CanFam3.txt

#Remake plink file 
#plink --bfile MERGED.clean.Fitak2018  --dog --update-map LiftOver/recodeSNPs_Fitak_CanFam3.txt  --exclude LiftOver/removeSites_Fitak_CanFam3.txt --make-bed  --out MERGED.clean.Fitak2018.LiftOverToCanFam3

######Update the rsIDs with the ones from cornell canine
#Get a list of those RSIDs that are shared
#awk 'FNR==NR{a[$1,$4]=$2;next}{print $0,a[$1,$4]?a[$1,$4]:"NA"}' /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/cornell_canine_updatedFID.bim MERGED.clean.Fitak2018.LiftOverToCanFam3.bim | awk '{print $2"\t"$7}' | grep -v "NA" > updateRSids.txt

#update into file plink file and delete old plink file
#plink --bfile MERGED.clean.Fitak2018.LiftOverToCanFam3 --dog --update-name LiftOver/updateRSids.txt --make-bed --out MERGED.clean.Fitak2018.LiftOverToCanFam3.updatedRSids


###Update the individual ID and FID to include breed and species
#plink --bfile MERGED.clean.Fitak2018.LiftOverToCanFam3.updatedRSids --dog --update-ids update_IIDs_Fatik_final.txt --make-bed --out MERGED.clean.Fitak2018.LiftOverToCanFam3.updatedRSids.updateIID

