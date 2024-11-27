library(tidyverse)
library(rtracklayer)
library(biomaRt)

getUTRs <- function(input="human/gencode.v46.basic.annotation.gff3.gz",transcripts.to.keep=NULL)
{
  gff3 <- readGFFAsGRanges(input)
  verified.loci <- gff3[gff3$type=="five_prime_UTR"] %>% # & (gff3$level == "1" | gff3$level == "2")] %>%
    as.data.frame()
  length.utr5 <- verified.loci$width
  utr.df <- data.frame(Transcript = verified.loci$transcript_id, Gene=verified.loci$gene_id,UTR5.Length=length.utr5)
  utr.df <- utr.df %>%
    mutate(Gene = str_remove(Gene,"\\.[0-9]+"))
  utr.df <- utr.df %>%
    group_by(Transcript) %>%
    summarize(UTR5.Length = sum(UTR5.Length),
              Gene=unique(Gene)) %>%
    mutate(Transcript = str_remove(Transcript,"\\.[0-9]+"))
  if (!is.null(transcripts.to.keep))
  {
    utr.df <- utr.df %>%
      filter(Transcript %in% transcripts.to.keep)
  }
  utr.df.summ <- utr.df %>%
    group_by(Gene) %>%
    summarize(Mean.UTR5.Length = mean(UTR5.Length),
              Median.UTR5.Length = median(UTR5.Length),
              Max.UTR5.Length = max(UTR5.Length))
  return(utr.df.summ)
}

human.iso <- read_tsv("../human_appris.tsv") %>%
  filter(MANE == "MANE_Select" & APPRIS_Annotation == "PRINCIPAL:1") %>%
  dplyr::select(Transcript_ID) %>%
  deframe()
mouse.iso <- read_tsv("../mouse_appris.tsv") %>%
  filter(APPRIS_Annotation == "PRINCIPAL:1") %>%
  dplyr::select(Transcript_ID) %>%
  deframe()


human.utr <- getUTRs(input="human/gencode.v46.basic.annotation.gff3.gz",human.iso)
mouse.utr <- getUTRs(input="mouse/gencode.vM35.basic.annotation.gff3.gz",mouse.iso)

mouse.id <- unique(mouse.utr$Gene)

mart1 <- useMart("ensembl", dataset="hsapiens_gene_ensembl",host = "dec2021.archive.ensembl.org")
mart2 <- useMart("ensembl", dataset="mmusculus_gene_ensembl",host = "dec2021.archive.ensembl.org") 

# human / mouse
mouse.to.human <- getLDS(attributes=c("ensembl_gene_id"),
                         filters="ensembl_gene_id", values=mouse.id, mart=mart2,
                         attributesL=c("ensembl_gene_id"), martL=mart1)%>%
  dplyr::rename(Mouse = Gene.stable.ID,Human = Gene.stable.ID.1) 


final.utr.df <- mouse.utr %>%
  full_join(mouse.to.human,by=c("Gene"="Mouse")) %>%
  full_join(human.utr,by=c("Human"="Gene"),suffix=c("_Mouse","_Human"))
  
write_tsv(final.utr.df,"human_vs_mouse_utr5_length_level_all_by_all_princ_isoforms.tsv")
