# `R_scripts/`

This directory contains some of the key R scripts used for this project. 

## `runMCMC_OUOU_parallel`

The primary script for running the MCMC analysis. Was primarily used along with the bash script `../runMCMC_OUOU.sh`. Can use arguments from the command line as input using the R package `argparse`. Also uses functions found in the `local_functions.R` script.

## `createDataFiles.R`

This script takes in the `Data/rna_abundances.tsv` and `Data/protein_abundances.tsv` (which were obtained from the supplemental material from [Ba et al. Science Advances 2022](https://doi.org/10.1126/sciadv.abn0756)) and processes it. Here, we did Z-scores for mRNA and protein abundances, but some other options are available in the script.

## `binGenesByExpression.R`

This script bins genes in deciles based on the species-wise mean of protein abundances. This was included as a supplemental result of the final version of the manuscript to complement analyses based on the PaxDB human protein abundances. Can be modified to use mRNA abundances. 