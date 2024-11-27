Model was run with bash command
```
Rscript --vanilla R_scripts/runMCMC_OUOU_parallel.R -i Data/mean_se_rna_protein.tsv -o Results/Drop_human/Pop_average_model_all//Independent//2024-10-16_independent_prior_based_on_fossil_tree_random_start/ --burnin_samples 0 --adapt_samples 0 --samples 50000 --prev_adapt_run Results/Drop_human/Pop_average_model_all//Independent//2024-10-16_independent_prior_based_on_fossil_tree_random_start//adaptive_mcmc.Rda --thin 5 --independent_evo --alpha_star 0.44 --num_cores 8 --drop_species Homo_sapiens --randomize_psi_start_sd 0 --newick --tree Data/vertlife_mcc_dna_only_fossil_bd.tre
```
