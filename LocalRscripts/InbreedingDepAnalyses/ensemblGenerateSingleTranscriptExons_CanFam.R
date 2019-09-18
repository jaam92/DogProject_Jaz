#Load libraries
library(dplyr)
library(tidyr)
####Read files in
genes = read.delim("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/EnsemblGenes_CanFam3.1.bed") #These genes come from downloading Ensembl genes from UCSC in bed format
gene_names = read.table("~/Documents/DogProject_Jaz/LocalRscripts/CaseControlROH/EnsemblGenes_CanFam3.1_geneNames.txt")

####Make file with genes of interest
####Filter to only chromosomes 1-38 
####keep only longest transcript (need to use distinct too bc some transcripts have equal length and we only want to keep one entry)

#Annotate gene set for Ensembl
GeneSet = genes %>% 
mutate(chrom = gsub("chr", "", chrom), 
       AbbrevName = gene_names$V2[match(name, gene_names$V1)], 
       transcript_length = transcriptionEnd - transcriptionStart) %>% 
  group_by(AbbrevName) %>% 
  filter(transcript_length == max(transcript_length) & as.numeric(chrom) <= 38) %>%
  distinct(AbbrevName,.keep_all= TRUE) %>% 
  as.data.frame()

#Split by exons
SplitExons = GeneSet %>%
  separate_rows(exonStarts, exonEnds) %>% #make each comma separation a separate row
  select(chrom, exonStarts, exonEnds, AbbrevName) %>%
  mutate(chrom = paste0("chr", chrom))

#Remove blank rows
FinalSplitExons = SplitExons[!(SplitExons$exonStarts == ""), ]

#write exons to file
write.table(FinalSplitExons, "~/Documents/DogProject_Jaz/InbreedingDepression/ForAbi_EnsemblGenes_CanFam3.1_SingleTranscript.bed", quote = F, row.names = F, col.names = F, sep = "\t")
