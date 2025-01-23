library(tidyverse)

cleanData <- function(df)
{
  df <- df %>%
    dplyr::select(GeneStableID,where(is.numeric))
  return(df)
}

convertZeroToNA <- function(df)
{
  df <- df %>% 
    mutate(across(where(is.numeric),~ifelse(. == 0,NA,.)))
  return(df)
}

transformLog <- function(df)
{
  df <- df %>% 
    mutate(across(where(is.numeric),~log(.)))
  return(df)
}

transformLog10 <- function(df)
{
  df <- df %>% 
    mutate(across(where(is.numeric),~log10(.)))
  return(df)
}

transformCenterMean <- function(df)
{
  df <- df %>% 
    mutate(across(where(is.numeric),~(. - mean(.,na.rm=T))))
  return(df)
}

transformCenterMedian <- function(df)
{
  df <- df %>% 
    mutate(across(where(is.numeric),~(. - median(.,na.rm=T))))
  return(df)
}


transformZScore <- function(df)
{
  df <- df %>% 
    mutate(across(where(is.numeric),~(. - mean(.,na.rm=T))/sd(.,na.rm=T))
    )
  return(df)
}

transformSqrt <- function(df)
{
  df <- df %>% 
          mutate(across(where(is.numeric),~sqrt(.)))
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


rna <- read_tsv("../Data/rna_abundances.tsv")
prot <- read_tsv("../Data/protein_abundances.tsv")

rna <- cleanData(rna)
prot <- cleanData(prot)

prot <- prot %>%
  dplyr::select(-contains("_mix"))

## Filtering before calculating Z-scores
rna <- rna %>%
  filter(GeneStableID %in% prot$GeneStableID)

prot <- prot %>%
  filter(GeneStableID %in% rna$GeneStableID)


rna <- convertZeroToNA(rna)
prot <- convertZeroToNA(prot)

rna <- transformLog(rna)
prot <- transformLog(prot)

# 
# rna <- transformZScore(rna)
# prot <- transformZScore(prot)


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

expr.mean.df <- expr.mean.df %>% 
  mutate(Mean_RNA = exp(Mean_RNA),
         Mean_Protein = exp(Mean_Protein))

write_tsv(expr.mean.se,"../Data/mean_rna_protein.tsv")
