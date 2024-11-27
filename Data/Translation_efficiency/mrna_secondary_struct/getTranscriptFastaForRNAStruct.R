library(tidyverse)
library(Biostrings)

species <- "human"
if (species == "mouse")
{
  target.fasta <- "gencode.vM35.pc_transcripts.fa.gz"
  gene.regex <- "ENSMUSG[0-9]+"
  transcript.regex <- "ENSMUST[0-9]+"
  output.file <- "mouse_transcripts_for_genes_in_ba_etal.fasta"
} else 
{
  target.fasta <- "gencode.v46.pc_transcripts.fa.gz"
  gene.regex <- "ENSG[0-9]+"
  transcript.regex <- "ENST[0-9]+"
  output.file <- "human_transcripts_for_genes_in_ba_etal.fasta"
}
  
genes.in.data <- read_tsv("../../mean_se_rna_protein.tsv")
genes.in.data <- unique(genes.in.data$Gene_ID)

if (species == "mouse")
{
  library(biomaRt)

  mart1 <- useMart("ensembl", dataset="hsapiens_gene_ensembl",host = "dec2021.archive.ensembl.org")
  mart2 <- useMart("ensembl", dataset="mmusculus_gene_ensembl",host = "dec2021.archive.ensembl.org") 
  
  # human / mouse
  mouse.to.human <- getLDS(attributes=c("ensembl_gene_id"),
                           filters="ensembl_gene_id", values=genes.in.data, mart=mart1,
                           attributesL=c("ensembl_gene_id"), martL=mart2)
  genes.in.data <- mouse.to.human %>%
    dplyr::rename(Human = Gene.stable.ID,Mouse= Gene.stable.ID.1) %>%
    filter(!is.na(Mouse)) %>%
    dplyr::select(Mouse) %>%
    unlist() %>%
    unique()

}

transcripts <- readDNAStringSet(target.fasta)
gene.ids <- str_extract(names(transcripts),gene.regex)
to.keep <- which(gene.ids %in% genes.in.data)
to.keep <- intersect(to.keep,which(str_detect(names(transcripts),"UTR5")))
transcripts.to.keep <- transcripts[to.keep]


cds.start <- as.numeric(str_extract(names(transcripts.to.keep),"(?<=CDS:)[0-9]+(?=-[0-9]+)"))
transcripts.to.keep <- transcripts.to.keep[which(cds.start >= 51)]
transcript.id <- str_extract(names(transcripts.to.keep),transcript.regex)
cds.loc <- str_extract(names(transcripts.to.keep),"CDS:[0-9]+-[0-9]+")
names(transcripts.to.keep) <- paste(transcript.id,cds.loc)
writeXStringSet(transcripts.to.keep,output.file)
