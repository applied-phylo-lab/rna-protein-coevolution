library(tidyverse)
library(biomaRt)

efe.around.start.codon <- T
filter.by.princ.iso <- T
use.chew <- F
getMRNASecStructPerSpecies <- function(file.path.ss,dataset = "hsapiens_gene_ensembl",pivot=T,tsv=F)
{
  if (tsv)
  {
    mrna.ss <- read_tsv(file.path.ss)
  } else {
    mrna.ss <- read_csv(file.path.ss)
  }
  mart <- useMart(biomart = "ensembl", dataset = dataset,host = "https://dec2021.archive.ensembl.org")
  
  res <- getBM(attributes = c('ensembl_transcript_id', 
                              'ensembl_gene_id'),
               filters = 'ensembl_transcript_id', 
               values = mrna.ss$TranscriptID,
               mart = mart)
  
  mrna.ss <- mrna.ss %>%
    right_join(res,by=c("TranscriptID" = "ensembl_transcript_id")) %>%
    relocate(ensembl_gene_id)
  
  if (pivot)
  {
    mrna.ss.long <- mrna.ss %>%
      pivot_longer(-c(ensembl_gene_id,TranscriptID),names_to = "Distance",values_to="CDS_SS_EFE")
  } else {
    mrna.ss.long <- mrna.ss %>%
      dplyr::rename(CDS_SS_EFE = ss_efe_profile)
  }
  return(mrna.ss.long)
}

# for Chew et al.
#human.ss <- getMRNASecStructPerSpecies("hs_35_RNA_fold_CDS_start_profiles_transpose.csv",dataset = "hsapiens_gene_ensembl")
#mouse.ss <- getMRNASecStructPerSpecies("mm_35_RNA_fold_CDS_start_profiles_transpose.csv",dataset = "mmusculus_gene_ensembl")

# for updated start codon analysis
human.ss <- getMRNASecStructPerSpecies("human_ss_around_start_site.tsv",dataset = "hsapiens_gene_ensembl",tsv=T,pivot=T)
mouse.ss <- getMRNASecStructPerSpecies("mouse_ss_around_start_site.tsv",dataset = "mmusculus_gene_ensembl",tsv=T,pivot=T)

if (use.chew)
{
  human <- readLines("human_ids_chew_etal.txt") 
  mouse <- readLines("mouse_ids_chew_etal.txt") 
  
  human.ss <- human.ss %>%
    filter(TranscriptID %in% human) %>%
    dplyr::select(-TranscriptID)
  mouse.ss <- mouse.ss %>%
    filter(TranscriptID %in% mouse) %>%
    dplyr::select(-TranscriptID)
  
  mouse.id <- unique(mouse.ss$ensembl_gene_id)
  
  mart1 <- useMart("ensembl", dataset="hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org")
  mart2 <- useMart("ensembl", dataset="mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org") 
  
  # human / mouse
  mouse.to.human <- getLDS(attributes=c("ensembl_gene_id"),
                           filters="ensembl_gene_id", values=mouse.id, mart=mart2,
                           attributesL=c("ensembl_gene_id"), martL=mart1)
  
  mouse.to.human.one.to.one <- mouse.to.human %>%
    dplyr::rename(Mouse = Gene.stable.ID,Human = Gene.stable.ID.1)
  
  mouse.ss <- mouse.ss %>%
    inner_join(mouse.to.human.one.to.one,by=c("ensembl_gene_id" = "Mouse"))
  
  mouse.ss <- mouse.ss %>%
    dplyr::select(-ensembl_gene_id)
  
  if (efe.around.start.codon){
    human.mouse.ss <- human.ss %>%
      left_join(mouse.ss,by=c("ensembl_gene_id" = "Human","Distance"),suffix = c("_Human","_Mouse")) %>%
      arrange(ensembl_gene_id)
  } else {
    human.mouse.ss <- human.ss %>%
      left_join(mouse.ss,by=c("ensembl_gene_id" = "Human"),suffix = c("_Human","_Mouse")) %>%
      arrange(ensembl_gene_id)
  }
} else {
  mouse.id <- unique(mouse.ss$TranscriptID)
  
  mart1 <- useMart("ensembl", dataset="hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org")
  mart2 <- useMart("ensembl", dataset="mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org") 
  
  # human / mouse
  mouse.to.human <- getLDS(attributes=c("ensembl_transcript_id"),
                           filters="ensembl_transcript_id", values=mouse.id, mart=mart2,
                           attributesL=c("ensembl_transcript_id"), martL=mart1)
  
  mouse.to.human.one.to.one <- mouse.to.human %>%
    dplyr::rename(Mouse = Transcript.stable.ID,Human = Transcript.stable.ID.1)
  
  mouse.ss <- mouse.ss %>%
    inner_join(mouse.to.human.one.to.one,by=c("TranscriptID" = "Mouse"))
  
  mouse.ss <- mouse.ss %>%
     dplyr::select(-ensembl_gene_id)
  
  if (efe.around.start.codon){
    human.mouse.ss <- human.ss %>%
      left_join(mouse.ss,by=c("TranscriptID" = "Human","Distance"),suffix = c("_Human","_Mouse")) %>%
      arrange(ensembl_gene_id)
  }
  
}

if (filter.by.princ.iso)
{
  human.iso <- read_tsv("../../../human_appris.tsv") %>%
    dplyr::select(Transcript_ID,APPRIS_Annotation,MANE)
  mouse.iso <- read_tsv("../../../mouse_appris.tsv") %>%
    dplyr::select(Transcript_ID,APPRIS_Annotation)
  human.mouse.ss <- human.mouse.ss %>%
    left_join(human.iso,by=c("TranscriptID" = "Transcript_ID")) %>%
    left_join(mouse.iso,by=c("TranscriptID_Mouse" = "Transcript_ID"),suffix=c("_Human","_Mouse")) 
  human.mouse.ss <- human.mouse.ss %>%
    #filter(APPRIS_Annotation_Human %in% c("PRINCIPAL:1","PRINCIPAL:2","PRINCIPAL:3") & APPRIS_Annotation_Mouse %in% c("PRINCIPAL:1","PRINCIPAL:2","PRINCIPAL:3") & !is.na(MANE))
    filter(APPRIS_Annotation_Human == "PRINCIPAL:1" & APPRIS_Annotation_Mouse == "PRINCIPAL:1" & !is.na(MANE))
  
} 

write_tsv(human.mouse.ss,"human_vs_mouse_start_codon_ss_transcript_id_appris_1_mane_select_2.tsv")
