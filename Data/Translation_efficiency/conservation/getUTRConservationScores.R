library(gridExtra)
library(grid)
library(ape)
library(stringr)
library(ggplot2)
library(stats)
library(splitstackshape)
library(dplyr)
library(GenomicScores)
library(seqinr)

df1 <-
  read.gff(
    "MANE.GRCh38.v1.4.ensembl_genomic.gff",
    na.strings = c(".", "?"),
    GFF3 = TRUE
  )[c(1, 3:5, 9)]

#Calculate mean 5' and 3' length
fives <- subset(df1, type == "five_prime_UTR")

library(GenomicScores)

phast <- getGScores("phyloP100way.UCSC.hg38")

#3. Obtain and average scores for each UTR in the 5' (Positive scores — Measure conservation, Negative scores — Measure acceleration.)
fives$PGS <- "NA" #Introduce new column with NA values
for (i in 1:(nrow(fives))) {
  #for all elements of dataset
  gs1 <-
    gscores(phast, GRanges(seqnames = (as.character(fives$seqid[i])), IRanges(
      start = (as.numeric(fives$start[i])):(as.numeric(fives$end[i])),
      width = 1
    ))) #Obtain a set of scores from within UTR coordinates
  gs1 <- na.omit(gs1) #do not count missing data (avoid skewing)
  scores1 <-
    gs1@elementMetadata@listData[["default"]] #Extract scores from GRanges object
  pgs1 <-
    (sum(scores1) / length(scores1)) #Divide the sum of scores by number of scores in the calculation
  fives$PGS[i] <- list(pgs1) #add to the output data the mean
}
fives$PGS <- as.numeric(fives$PGS)


#split attributes out
five_split <- cSplit(fives, "attributes", ";")
#remove titles
names(five_split)[8] <- "gene_id"
names(five_split)[9] <- "transcript_id"
names(five_split)[11] <- "gene"
names(five_split)[14] <- "exon_number"
names(five_split)[16] <- "mane_type"

titles <-
  c("gene_id", "transcript_id", "gene", "exon_number", "mane_type")

five_split[, titles] <-
  lapply(five_split[, titles, with = FALSE], function(x)
    gsub("\\w+=", "", x))

#want mane select only
five_split_mane <-
  five_split[grepl("MANE_Selec", five_split$mane_type),]

#now need to get PGS scores for whole UTR. now they are per 5'UTR exon
phylop_final <-
  five_split_mane[, .(agg_PGS = sum(PGS) / .N), by = transcript_id]


phylop_final <- phylop_final %>%
  as.data.frame() %>%
  mutate(transcript_id = str_remove(transcript_id,"\\.[0-9]+"))
mart1 <- useMart("ensembl", dataset="hsapiens_gene_ensembl",host = "dec2021.archive.ensembl.org")
id.map <- getBM(attributes=c('ensembl_transcript_id',"ensembl_gene_id"),
                values=phylop_final$transcript_id,
                filters="ensembl_transcript_id", 
                mart = mart1)

phylop_final_w_gene_id <- phylop_final %>% 
  left_join(id.map,by=c("transcript_id"="ensembl_transcript_id")) %>%
  filter(!is.na(ensembl_gene_id) & !is.na(agg_PGS)) %>%
  dplyr::select(ensembl_gene_id,agg_PGS)

write_tsv(phylop_final_w_gene_id,"2024-11-11_human_utr5_phylop.tsv")

