# `Results/`

All results presented in the manuscript are based on the results from dropping human from our analysis. This directory contains the result of our analysis, included figures showing traces from the MCMC. Currently, objects from `LaplacesDemon` are saved as .Rda files, which can be quite large. We will soon add FLAT files containing the summary parameters (some runs may already have these). We are happy to provide these for specific runs on requests. Each run should also contain a README.md file providing the command used to generate the analyses. 

## `Results/Drop_human/`

1. `Pop_average_model_all` - analyses fitting all 3 models (protein-driven, mRNA-driven, and independent) to all genes at the same time
2. `Essential_vs_non_essential` - analyses fitting the protein-driven model to essential vs. non-essential genes
3. `Haplosufficient_vs_insufficient` - analyses fitting the protein-driven model to examine the impact of haploinsufficiency on gene expression evolution
4. `Bin_based_on_paxdb_human_prot_abund` - analyses fitting the protein-driven model across gene expression decile bins based on human protein abundances obtained from PaxDB, shows how parameters vary with gene expression
5. `Bin_based_on_paxdb_human_prot_abund_essential` - same as `Bin_based_on_paxdb_human_prot_abund`, but using only genes marked as essential
6. `Bin_based_on_species_mean_protein` - same as `Bin_based_on_paxdb_human_prot_abund`, but using across species mean protein abundances for determining the binds
