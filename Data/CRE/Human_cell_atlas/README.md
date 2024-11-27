CREs were downloaded from the [Human Cell Atlas](https://doi.org/10.1101/2023.11.13.566791), specifically [Supplemental Table 2](https://doi.org/10.6084/m9.figshare.c.6926944). This file was converted to a CSV file using Microsoft Excel and then processed using the following command in R

```
human_atlas <- df %>% 
	dplyr::select(tCRE_ID,tCRE_type,link_gene_ID) %>%
	separate(link_gene_ID,sep=";",into=paste("Test",1:93),fill="right") %>% 
	pivot_longer(starts_with("Test"),names_to="Case",values_to="Gene_ID") %>% 
	filter(!is.na(Gene_ID) & Gene_ID != "unlinked") %>% 
	dplyr::select(-Case) %>%
	group_by(Gene_ID) %>% 
	summarize(Num.Enhancers=n()) %>% 
	mutate(Gene_ID = str_remove(Gene_ID,"\\.[0-9]+"))

```

where `df` was the resulting `data.frame` read in using the function `tidyverse::read_csv`. Note that using `Num.Enhancers` is just to make these files consistent with the files analyzed from EnhancerAtlas 2.0, but in reality this includes both promoters and enhancers. 