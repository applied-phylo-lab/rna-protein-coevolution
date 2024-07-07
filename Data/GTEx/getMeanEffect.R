library(tidyverse)



gtex.files <- list.files("GTEx_Analysis_v8_eQTL/",pattern = "signif_variant_gene_pairs.txt",full.names = T)

gtex.df <- purrr::map(gtex.files,read_tsv) %>%
  bind_rows()

gtex.df <- gtex.df %>%
  mutate(gene_id = str_remove(gene_id,"\\.[0-9]+"))

gtex.summary <- gtex.df %>%
  group_by(gene_id) %>%
  summarize(Num.eQTL = n(),
            Mean.effect = mean(slope),
            Mean.abs.effect = mean(abs(slope)))
write_tsv(gtex.summary,"eQTL_summary.tsv")