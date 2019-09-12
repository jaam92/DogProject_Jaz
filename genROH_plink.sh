

#Call ROHs using parameters from biorxiv pre-print: Fine-scale resolution and analysis of runs of homozygosity in domestic dogs by Aaron Sams and Adam Boyko

#Should run in about 1 min

#Load plink
#module load plink

#Original file no LD pruning
plink --dog --bfile cornell_canine_updatedFID --not-chr 39,41 --homozyg-window-het 0 --homozyg-snp 41 --homozyg-window-snp 41 --homozyg-window-missing 0 --homozyg-window-threshold 0.05 --homozyg-kb 500 --homozyg-density 5000 --homozyg-gap 1000 --out ROH_ogFile_autosomesOnly

#File that is LD pruned
plink --dog --bfile smartPCA/cornell_canine_smartPCA_LDpruned --not-chr 39,41 --homozyg-window-het 0 --homozyg-snp 41 --homozyg-window-snp 41 --homozyg-window-missing 0 --homozyg-window-threshold 0.05 --homozyg-kb 500 --homozyg-density 5000 --homozyg-gap 1000 --out ROH_ldPrunedFile_autosomesOnly
