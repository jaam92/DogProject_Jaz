#split file by breed
#awk '{print>$2".txt"}' ../Individuals_allBreeds_mergedFitakCornell.txt

#grab file names by minimum sample size
#for i in {6..10}; do wc -l *.txt | grep -v "total" | awk '$1>='$i' {print $2}' > sampSize_grEqlto"$i"_fileNames.txt;done
