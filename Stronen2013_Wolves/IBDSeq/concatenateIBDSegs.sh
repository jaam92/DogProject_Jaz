for f in {1..38}; do awk '$1 != $3 {print $0}' Stronen2013_wolves_chr"$f"_Haplotypes_IBDSeq.ibd >> Stronen2013Wolves_allChroms_Haplotypes_IBDSeq.ibd;done
