# rna-protein-coevolution

This repository contains all data, code, and results necessary to recreate the results of Macroevolutionary divergence of gene expression driven by selection on protein abundance by Cope, Schraiber, and Pennell. A link to an early version of the manuscript posted as a preprint can be found [here](https://doi.org/10.1101/2024.07.08.602411).


## Running the code

All code can be run in the R programming environment. Most figures were generated using RStudio (version 2024.04.0+735) with `R 4.4.0.`. Below is a table of key software and packages necessary to run the MCMCs.

| Software      	| version     	|
|---------------	|-------------	|
| PCMBase	     	| 1.2.13        |
| PCMBaseCpp       	| 0.1.9       	|
| LaplacesDemon     | 16.1.6       	|
| phytools          | 2.3-0      	|
| ape               | 5.8       	|
| argparse          | 2.2.3         |
| tidyverse         | 2.0.0         |
| R             	| 4.4.0       	|

We note we used some local versions of code from `LaplacesDemon` to correct errors in the code and expand the functionality to help speed up our MCMC runs. These local version can be found in `R_scripts/local_functions.R`. The bash script `runMCMC_OUOU.sh` can be modified to specify the setup for the MCMC analyses.

## `Data/`

The `Data/` directory contains the necessary processed data and phylogenetic trees necessary to recreate the analyses described in Cope et al.

## `Results/`

Contains all results applied to real data presented in the manuscript. Each result will contain a `README.md` file that contains the command used to run the MCMC.

## `Model_testing/`

The `Model_testing/` directory contains simulated data and model fits to simulated data to assess that the MCMC implementation is working correctly. Figures for the model testing can be created using the `compareEstimatesToTruth.Rmd`.  