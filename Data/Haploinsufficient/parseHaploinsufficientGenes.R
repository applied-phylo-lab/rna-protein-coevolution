library(tidyverse)

haplo.genes <- read_tsv("haploinsufficient.tsv")
ensembl.id <- haplo.genes %>% 
  dplyr::select(Ensembl_Gene_ID	) %>%
  deframe()

mrna.prot <- read_tsv("../mean_se_rna_protein.tsv")

mrna.prot.haplo <- mrna.prot %>% 
  filter(Gene_ID %in% ensembl.id)

mrna.prot.non.haplo <- mrna.prot %>% 
  filter(!Gene_ID %in% ensembl.id)

write_tsv(mrna.prot.haplo,"haploinsufficient_mean_se_rna_protein.tsv")
write_tsv(mrna.prot.non.haplo,"haplosufficient_mean_se_rna_protein.tsv")
