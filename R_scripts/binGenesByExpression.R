library(tidyverse)
library(Hmisc)

data <- read_tsv("../Data/mean_se_rna_protein.tsv")
data.summ <- data %>%
  filter(Species != "Homo_sapiens") %>%
  group_by(Gene_ID) %>%
  dplyr::summarize(Mean_Protein = mean(Mean_Protein),
            Mean_RNA = mean(Mean_RNA))

data.summ <- data.summ %>% 
  mutate(Bin = cut2(Mean_Protein,g = 5))

data.split <- data %>%
  left_join(data.summ %>% dplyr::select(Gene_ID,Bin)) %>%
  group_by(Bin) %>%
  group_split()

bin.names <- c("0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80-90","90-100")
#bin.names <- c("0-20","20-40","40-60","60-80","80-100")


lapply(1:length(bin.names),function(x)
  {
   write_tsv(data.split[[x]] %>% dplyr::select(-Bin),paste0("../Data/Bin_based_on_mrna/5_groups/mean_se_rna_protein_",bin.names[x],".tsv"))
})

  