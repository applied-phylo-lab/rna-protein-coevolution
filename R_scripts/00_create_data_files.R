library(tidyverse)

cleanData <- function(df)
{
  df <- df %>%
    dplyr::select(GeneStableID,where(is.numeric))
  return(df)
}

calculateSummary <- function(df)
{
  df <- df %>% 
    pivot_longer(-GeneStableID,names_to="Measurement",values_to="Abundance") %>%
    separate(col=Measurement,into=c("Species","Replicate"),sep="_") %>% 
    mutate(Species = str_remove(Species,"SF"),
           Replicate = str_remove(Species,"SF")) %>%
    group_by(GeneStableID,Species) %>%
    summarize(Mean = mean(Abundance,na.rm=T),
              SD = sd(Abundance,na.rm=T),
              Missing = sum(is.na(Abundance)),
              SE = SD/sqrt(n() - Missing)) %>%
    mutate(Species = case_when(
      Species == "Cat" ~ "Felis_catus",
      Species == "Cow" ~ "Bos_taurus",
      Species == "Dog" ~ "Canis_lupus",
      Species == "Horse" ~ "Equus_caballus",
      Species == "Human" ~ "Homo_sapiens",
      Species == "Monkey" ~ "Macaca_mulatta",
      Species == "Sheep" ~ "Ovis_aries",
      Species == "Pig" ~ "Sus_scrofa",
      Species == "Rat" ~ "Rattus_norvegicus",
      Species == "Opossum" | Species == "Opposum" ~ "Monodelphis_domestica",
      Species == "Rabbit" ~ "Oryctolagus_cuniculus"
    ))
  return(df)
}



tree <- read.nexus("../Data/tree_11sp_noGpig.nex")
tree <- force.ultrametric(tree,method = "extend")

rna <- read_tsv("../Data/rna_abundances.tsv")
prot <- read_tsv("../Data/protein_abundances.tsv")


rna <- cleanData(rna)
prot <- cleanData(prot)
prot <- prot %>%
  dplyr::select(-contains("_mix"))
rna <- rna %>% 
  filter(GeneStableID %in% prot$GeneStableID)


total.tpm.per.species <- colSums(rna[,2:ncol(rna)])

rna <- rna %>% 
  mutate(across(where(is.numeric),~ifelse(. == 0,NA,.)),
         across(where(is.numeric),~log(.)),
         across(where(is.numeric),~(. - mean(.,na.rm=T))/sd(.,na.rm=T))
  )



prot <- prot %>% 
  mutate(across(where(is.numeric),~ifelse(. == 0,NA,.)),
         across(where(is.numeric),~log(.)),
         across(where(is.numeric),~(. - mean(.,na.rm=T))/sd(.,na.rm=T))
  )


rna.summary <- calculateSummary(rna)
prot.summary <- calculateSummary(prot)

rna.missing <- rna.summary %>% 
  ungroup() %>%
  group_by(GeneStableID) %>%
  summarize(Total.Missing = sum(Missing)) %>%
  filter(Total.Missing == 0)

prot.missing <- prot.summary %>% 
  ungroup() %>%
  group_by(GeneStableID) %>%
  summarize(Total.Missing = sum(Missing)) %>%
  filter(Total.Missing == 0)

rna.summary <- rna.summary %>%
  filter(GeneStableID %in% rna.missing$GeneStableID) %>%
  dplyr::rename(Gene_ID = GeneStableID) 
prot.summary <- prot.summary %>%
  filter(GeneStableID %in% prot.missing$GeneStableID) %>%
  dplyr::rename(Gene_ID = GeneStableID) 

expr.mean.df <- rna.summary %>% 
  dplyr::select(Gene_ID,Species,Mean) %>%
  inner_join(prot.summary %>% dplyr::select(Gene_ID,Species,Mean),
            by=c("Gene_ID","Species"),
            suffix = c("_RNA","_Protein"))

expr.se.df <- rna.summary %>% 
  dplyr::select(Gene_ID,Species,SE) %>%
  inner_join(prot.summary %>% dplyr::select(Gene_ID,Species,SE),
            by=c("Gene_ID","Species"),
            suffix = c("_RNA","_Protein"))


expr.mean.se <- expr.mean.df %>%
  left_join(expr.se.df,by=c("Gene_ID","Species"))

write_tsv(expr.mean.se,"../Data/mean_se_rna_protein.tsv")
