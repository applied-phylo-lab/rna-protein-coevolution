library(tidyverse)

essential.genes <- read_csv("essential.csv")
ensembl.id <- essential.genes %>% 
  dplyr::select(ensg) %>%
  deframe()

mrna.prot <- read_tsv("../mean_se_rna_protein.tsv")

mrna.prot.ess <- mrna.prot %>% 
  filter(Gene_ID %in% ensembl.id)

mrna.prot.non.ess <- mrna.prot %>% 
  filter(!Gene_ID %in% ensembl.id)

write_tsv(mrna.prot.ess,"essential_mean_se_rna_protein.tsv")
write_tsv(mrna.prot.non.ess,"non_essential_mean_se_rna_protein.tsv")
