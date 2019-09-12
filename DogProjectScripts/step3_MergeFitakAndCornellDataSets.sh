#The first time we attempt to merge there is an error because of strand flips on the chip
#We ID the errors and then we will correct them in step 2
#In step 3 we will create the final merged vcf file

#Load plink
#. /u/local/Modules/default/init/modules.sh
#module load plink

#Step 1 (id allele flip errors)
#plink --dog --bfile FitakData/MERGED.clean.Fitak2018.LiftOverToCanFam3.updatedRSids.updateIID --bmerge CornellCanine/cornell_canine_updatedFID.bed CornellCanine/cornell_canine_updatedFID.bim CornellCanine/cornell_canine_updatedFID.fam --make-bed --out FindAlleleFlips

#Step 2 fix errors
#plink --bfile ../FitakData/IntermidiateFile/MERGED.clean.Fitak2018.LiftOverToCanFam3.updatedRSids.updateIID --dog --flip FindAlleleFlips-merge.missnp --make-bed --out ../FitakData/FitakData_FinalCanFam3

#Step 3 final merged file, only keep sites where at least 90% of individuals have GT
# plink --dog --bfile ../FitakData/FitakData_FinalCanFam3 --bmerge ../CornellData/cornell_canine_updatedFID.bed ../CornellData/cornell_canine_updatedFID.bim ../CornellData/cornell_canine_updatedFID.fam --geno 0.1 --make-bed --out IntermidiateFile/MergedFile_CornellCanineFitak

#Step 4 identify sets of unrelated indiviudals
#module load R/3.4.2
#Rscript ID_duplicateIndivs.R 

#Step 5 subset down to a file without duplicates and a file with  unrelateds only
#plink --bfile IntermidiateFile/MergedFile_CornellCanineFitak --dog --keep Individuals_allBreeds_mergedFitakCornell.txt --make-bed --out MergedFile_CornellCanineFitak_allIndivs
#plink --bfile IntermidiateFile/MergedFile_CornellCanineFitak --dog --keep UnrelatedIndividuals_allBreeds_mergedFitakCornell.txt --make-bed --out MergedFile_CornellCanineFitak_UnrelatedsOnly

