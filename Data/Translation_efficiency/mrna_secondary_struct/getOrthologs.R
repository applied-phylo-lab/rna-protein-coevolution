library(tidyverse)
library(biomaRt,lib.loc = "~/R_dev/")

data <- read_tsv("../../../../mean_se_rna_protein.tsv")
human.appris <- read_tsv("../human_appris.tsv") %>%
  dplyr::select(Gene_ID,Transcript_ID,APPRIS_Annotation)
mouse.appris <- read_tsv("../mouse_appris.tsv") %>%
  dplyr::select(Gene_ID,Transcript_ID,APPRIS_Annotation)



mart.1 <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org")

res <- getBM(attributes = c('ensembl_transcript_id', 
                            'ensembl_gene_id'),
             filters = 'ensembl_gene_id', 
             values = unique(data$Gene_ID),
             mart = mart.1)


mart.2 <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org")


mouse.to.human <- getLDS(attributes=c("ensembl_transcript_id","ensembl_gene_id"),
                         filters="ensembl_transcript_id", values=unique(res$ensembl_transcript_id), 
                         mart=mart.1,
                         attributesL=c("ensembl_transcript_id","ensembl_gene_id"), 
                         martL=mart.2)
mouse.to.human <- mouse.to.human %>%
  dplyr::rename(Human.Transcript.ID = Transcript.stable.ID,
                Human.Gene.ID = Gene.stable.ID,
                Mouse.Transcript.ID = Transcript.stable.ID.1,
                Mouse.Gene.ID = Gene.stable.ID.1) 

mouse.to.human.w.appris <- mouse.to.human %>%
  left_join(human.appris,by=c("Human.Transcript.ID" = "Transcript_ID","Human.Gene.ID" = "Gene_ID")) %>%
  left_join(mouse.appris,by=c("Mouse.Transcript.ID" = "Transcript_ID","Mouse.Gene.ID" = "Gene_ID"),suffix=c("_Human","_Mouse"))


isoform.per.gene.human <- mouse.to.human.w.appris %>%
  group_by(Human.Gene.ID,Human.Transcript.ID) %>%
  summarize(Isoform.per.gene = n())

isoform.per.gene.mouse <- mouse.to.human.w.appris %>%
  group_by(Mouse.Gene.ID,Mouse.Transcript.ID) %>%
  summarize(Isoform.per.gene = n())


test <- mouse.to.human.w.appris %>% 
  filter(!is.na(APPRIS_Annotation_Human) & !is.na(APPRIS_Annotation_Mouse)) %>%
  left_join(isoform.per.gene.human,by = c("Human.Transcript.ID", "Human.Gene.ID")) %>% 
  left_join(isoform.per.gene.mouse,by = c("Mouse.Transcript.ID", "Mouse.Gene.ID"),suffix=c(".Human",".Mouse"))
