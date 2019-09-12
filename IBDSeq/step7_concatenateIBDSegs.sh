#concatenate IBD segments for IBDNe
#for f in {1..38}; do awk '$1 != $3 {print $0}' MergedFile_CornellCanineFitak_UnrelatedsOnly_chr"$f"_Haplotypes_IBDSeq.ibd >> MergedFitakCornell_allChroms_Haplotypes_IBDSeq.ibd;done
