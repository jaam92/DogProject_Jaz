This folder contains scripts for making a GRM using overlapping ROH within the Hayward data
This script will count up the total number of base-pairs of overlapping ROHs between two target individuals that are being compared.
The output contains the two individuals being compared and the total number of bp of shared ROH between the two samples.
Final output ROHGRM.ped is the concatenated data across all 942 comparison files with a header.


Now to the scripts and files:

The split files are the exact same split files that we are using to Permute ROH so I simply copied this directory: /u/home/j/jmooney3/klohmueldata/jazlyn_data/DogProject_Jaz/InbreedingDepression/PermuteROH/HaywardData/vcftools/splitFiles
The file fnames_haywardComps.txt contains all pairwise comparisons of individuals from the Hayward data, minus comparing to same sample. I got this file from my Rscript visualizeROHBurdenandKinshipMatrices.R, which is stored in LocalRscripts.

Step 1:
run script step1_splitFnames_comparison.sh
This script will split fnames_haywardComps.txt file into piecies so there are only 20000 comparisons in each file
Next, rename the files as fnames_haywardComps, then change extension 000 to 942, and remove the 0s from the begining of digits in the filename 
Lastly, move these files into a directory called inputs 

Step 2:
make directory called Output 
submit script to do comparisons as follows qsub -t 1:942 step2_pairwiseComps_array.sh
