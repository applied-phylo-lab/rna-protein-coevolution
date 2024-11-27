library(tidyverse)

only.essential <- F
only.low.se.genes <- T

hsapien.prot <- read_tsv("paxdb_human_whole_organism_prot_abundance_integrated.tsv",comment="#") %>%
  mutate(string_external_id = str_remove(string_external_id,"9606."))

essential <- read_csv("../Essential/essential.csv")
low.se.genes <- readLines("../low_se_genes.txt")

ensembl.prot.to.gene <- read_tsv("hsapien_gene_mapping.tsv")

hsapien.prot <- hsapien.prot %>%
  left_join(ensembl.prot.to.gene,by=c("string_external_id"="ensembl_peptide_id"))

ba.data <- read_tsv("../mean_se_rna_protein.tsv")
ba.genes <- ba.data %>%
  dplyr::select(Gene_ID) %>%
  deframe() %>%
  unique()

hsapien.prot.filt <- hsapien.prot %>%
  filter(ensembl_gene_id %in% ba.genes) %>% 
  group_by(ensembl_gene_id) %>% 
  summarize(abundance = sum(abundance))

missing.ids <- ba.genes[which(!ba.genes %in% hsapien.prot.filt$ensembl_gene_id)]

test <- ba.data %>% 
  filter(Species == "Homo_sapiens") %>%
  left_join(hsapien.prot.filt,by=c("Gene_ID" = "ensembl_gene_id"))

hsapien.prot.filt <- hsapien.prot.filt %>%
  mutate(log10.abundance = log10(abundance),
         Bin = Hmisc::cut2(log10.abundance,g = 10))



data.split <- ba.data %>%
  left_join(hsapien.prot.filt %>% dplyr::select(ensembl_gene_id,Bin),by=c("Gene_ID"="ensembl_gene_id")) %>%
  filter(!is.na(Bin)) 

if (only.low.se.genes)
{
  data.split <- data.split %>%
    filter(Gene_ID %in% low.se.genes)
}

if (only.essential)
{
  data.split <- data.split %>%
    filter(Gene_ID %in% essential$ensg) %>% 
    group_by(Bin) %>%
    group_split()
} else{
  data.split <- data.split %>%
    group_by(Bin) %>%
    group_split()
}
  


bin.names <- c("0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80-90","90-100")


lapply(1:length(bin.names),function(x)
{
  write_tsv(data.split[[x]] %>% dplyr::select(-Bin),paste0("../Bin_based_on_paxdb_human_prot_abund/10_groups_low_se/mean_se_rna_protein_",bin.names[x],".tsv"))
})




