#Make final file with updated FID and remove duplicate snps (kept a single occurance) and remove a mislabelled individual (says rhodesian ridgeback in dryad but it is actually a golden retriever)
#. /u/local/Modules/default/init/modules.sh
#module load plink
#plink --bfile cornell_canine --dog --exclude excludeSNPs.txt --remove removeIndividual.txt --update-ids recode.txt --make-bed --out ../CornellData/cornell_canine_updatedFID
