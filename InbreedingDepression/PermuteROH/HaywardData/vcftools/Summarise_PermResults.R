library(tidyverse)

perms = read_delim(gzfile("Merge_CountPermutedOverlaps_100Kb_AutosomalSplits_bedtoolsShuffle.bed.gz"), delim = "\t", col_names = c("start", "stop", "chromosome", "overlaps")) %>%
  select(chromosome,start,overlaps) %>%
  group_by(chromosome, start) %>%
  mutate(proportion = overlaps/4342) %>%
  summarize(avg=mean(proportion),
            n=n(),
            sd=sd(proportion),
            se=sd/sqrt(n)) %>%
  mutate(chrom = as.numeric(gsub("chr", "", chromosome)),
         upper = avg + 1.96*sd,
         lower = avg - 1.96*sd) %>%
  arrange(chrom)

write.table(perms, file="SummaryFile_Merge_CountPermutedOverlaps_100Kb_AutosomalSplits_bedtoolsShuffle.txt",quote=F, sep="\t", row.names=F)


