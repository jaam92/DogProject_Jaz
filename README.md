# The impact of identity by descent on fitness and disease in dogs

Code and summary files accompanying:

> Mooney, J. A., Yohannes, A., & Lohmueller, K. E. (2021). The impact of identity by descent on fitness and disease in dogs. *Proceedings of the National Academy of Sciences*, 118(16), e2019116118. [Full text on PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC8072400/)

This repository contains the scripts used to merge and process SNP array data from dogs and wolves, call runs of homozygosity (ROH) and identity-by-descent (IBD) segments, infer recent effective population size, and test how ROH and IBD relate to breed popularity, disease, and inbreeding depression.

## Data availability

| Resource | Location |
|---|---|
| SNP genotype data from the original projects | Dryad: https://datadryad.org/stash/dataset/doi:10.5061/dryad.g68k008 and https://datadryad.org/stash/dataset/doi:10.5061/dryad.266k4 |
| Merged dog and wolf dataset used in the paper | This repository: `MergeFitkakAndCornell/FinalFiles.tar.bz2` (see below) |
| Code to process the original data and generate intermediate files | This repository |

In the code, the original datasets are referred to as the Cornell/Hayward data (Hayward et al. 2016) and the Fitak data (Fitak et al. 2018). A separate wolf dataset (Stronen et al. 2013, lifted over from CanFam2 to CanFam3.1) is processed in `Stronen2013_Wolves/`.

## Repository structure

### Processing pipeline (cluster scripts)

The numbered scripts in `DogProjectScripts/` run the main pipeline in order. Copies of several of them also sit in the top level of the repository and in analysis-specific folders.

| Step | Scripts | What it does |
|---|---|---|
| 1–3 | `step1_updatedCornellCanine.sh`, `step2_processFitak2018.sh`, `step3_MergeFitakAndCornellDataSets.sh` | Clean the Hayward and Fitak datasets (update IDs, remove duplicate SNPs and a mislabeled individual) and merge them. |
| 4–5 | `step4_plink2vcf.sh`, `step5_genROHvcfTools.sh` | Convert PLINK files to VCF and call ROH with VCFtools. |
| 6–9 | `step6_runIBDSeq.sh` through `step9B_runIBDNE_IBDSeq_Sept15ver.sh` | Call IBD segments with IBDSeq, concatenate them, and estimate recent effective population size with IBDNe (two IBDNe versions). |
| 10–13 | `step10_subsetToStronen2013Analyses.sh` through `step13B_...` | Repeat ROH, IBD, and IBDNe on dogs restricted to the sites shared with the Stronen et al. 2013 wolf data, for dog–wolf comparisons. |
| 14–15 | `step14_runIBDSeq_allIndivs4CaseControl.sh`, `step15_consensusROH_allIndivs4CaseControl.sh` | IBD and consensus ROH in all individuals for case–control analyses. |

Other cluster-side folders:

| Directory | Contents |
|---|---|
| `ogFiles_Hayward2016/`, `ogFiles_Fitak2018/` | Processing of each original dataset. |
| `MergeFitkakAndCornell/` | Merging scripts and the final merged files. |
| `CornellData/`, `IndividualFiles/` | Identifying unrelated individuals and splitting samples by breed. |
| `vcfToolsROH/`, `IBDSeq/`, `IBDNE/` | ROH, IBD, and IBDNe runs on the merged data. |
| `SubsetStronen2013Sites/`, `Stronen2013_Wolves/` | Analyses restricted to Stronen 2013 sites, and processing of the wolf data. |
| `DogProjectScripts/PhasingSHAPEIT/` | Phasing with SHAPEIT and IBD calling with Beagle, used as a comparison with IBDSeq. |
| `CalcFst/` | Pairwise F<sub>ST</sub> between breeds. |
| `ComputeG12/` | G12 haplotype homozygosity scans. |
| `InbreedingDepression/` | Intersecting ROH with exons and GERP scores, adding HGNC gene names, and permuting ROH across the genome to find regions depleted of ROH. |
| `ROHGRM/` | A relationship matrix based on shared ROH (see `README_ROHGRM.txt`). |
| `DidNotUse_SamsData/` | Exploratory scripts not used in the final paper. |

### Analysis and plotting (`LocalRscripts/`)

All R scripts used for the analyses and figures are in `LocalRscripts/`, with subdirectories labeled by analysis: `AKC` (breed popularity and effective population size), `BreedCladeInfo`, `CaseControlROH` (ROH burden and association tests for disease phenotypes), `Fst`, `IBDNe`, `IBDSegs` (IBD scores), `InbreedingDepAnalyses` (ROH depletion near exons and validation of candidate genes), `OMIA` (causal disease variants and disease prevalence), `PCA_Unrelateds`, `ROH` (ROH scores and comparisons with IBD), and `dog_genetic_maps`.

## Summary files

| File | Description |
|---|---|
| `MergeFitkakAndCornell/FinalFiles.tar.bz2` | Final merged PLINK files (`MergedFile_CornellCanineFitak_allIndivs.bed/.bim/.fam`) used in the paper, plus breed and clade assignments per dog (clades from Parker et al. 2004), the list of unrelated individuals (no closer than third-degree relatives by SNPRelate KING kinship), and the geographic origin of each wolf. A `README.info` inside describes each file. |
| `LocalRscripts/ROH/TrueROH_propCoveredwithin1SDMean_allChroms_mergedFile_Cornell_allChroms_vcfToolsROH_rmROHlessThan50snps_HaywardDataOnly.txt` | Filtered ROH calls (VCFtools) for the Hayward data, with ROH containing fewer than 50 SNPs removed: start and end positions, length, individual, chromosome, proportion of the ROH covered by SNPs, and whether the ROH was retained. |
| `LocalRscripts/OMIA/processedCausalVarsOMIA.txt` | Causal disease variants by breed from Online Mendelian Inheritance in Animals (OMIA): OMIA ID, breed, gene, phenotype, and NCBI gene ID. |
| `LocalRscripts/AKC/AKC_breedPopularity_1926thru2005.txt` | American Kennel Club breed popularity from 1926 to 2005, per breed. |
| `LocalRscripts/BreedCladeInfo/OrderCluster.txt` | Ordering of breed clades used for plotting. |
| `LocalRscripts/InbreedingDepAnalyses/VettingResults/SupplementaryFileFormat_ExonRegion_NonOverlapsROH_cornellData.bed` | Genes with at least one exon not overlapped by any ROH: gene name, chromosome, and exon start and end coordinates. |

## Software

Command line: PLINK, VCFtools, IBDSeq (r1206), IBDNe (04Sep15 and 07May18 releases), SHAPEIT, Beagle, bedtools, liftOver, and SnpEff.

R: tidyverse, data.table, SNPRelate, GENESIS, GWASTools, GMMAT, GenometriCorr, GenomicRanges, biomaRt, ggpubr, cowplot, and others loaded at the top of each script. See the Methods and SI Appendix of the paper for versions and parameters.

## Usage notes

The cluster scripts were written for UCLA's Hoffman2 cluster and contain hard-coded paths and `module load` commands; update these before running. R scripts in `LocalRscripts/` were run locally and may also need their working directories changed.

## Citation

If you use this code or data, please cite the paper:

```bibtex
@article{mooney2021dogibd,
  title   = {The impact of identity by descent on fitness and disease in dogs},
  author  = {Mooney, Jazlyn A. and Yohannes, Abigail and Lohmueller, Kirk E.},
  journal = {Proceedings of the National Academy of Sciences},
  volume  = {118},
  number  = {16},
  pages   = {e2019116118},
  year    = {2021},
  doi     = {10.1073/pnas.2019116118}
}
```

## Contact

For questions, please open an issue on this repository or contact the corresponding authors listed in the paper.
