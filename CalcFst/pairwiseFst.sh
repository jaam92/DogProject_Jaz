
#load plink

#while read p
#do

#fname=${p// /_vs_}
#echo $fname

#plink --bfile ancestry/MergedFile_CornellCanineFitak_Unrelateds_grEql30 --dog --within make.clst --keep-cluster-names $p --fst --out $fname 

#done < pairwise.clst

#Collected the pairwise fst measures per comparison
# grep "Weighted" *.log | sed 's/.log:Weighted Fst estimate: /\t/g' | sed -e 's/_vs_/\t/g' | sed 1i'breed1\tbreed2\tWeightedFst' > pairwiseWeighted.fst
