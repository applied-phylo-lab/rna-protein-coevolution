Model was run with bash command
```
Rscript --vanilla R_scripts/runMCMC_OUOU_parallel.R -i Data/Essential/non_essential_mean_se_rna_protein.tsv -o Results/Drop_human/Essential_vs_non_essential/Non_essential//Independent//2024-06-13_independent_prior_based_on_tree_random_start/ --burnin_samples 1000 --adapt_samples 25000 --samples 50000 --thin 5 --independent_evo --alpha_star 0.44 --num_cores 4 --drop_species Homo_sapiens --randomize_psi_start_sd 0 --tree Data/tree_11sp_noGpig.nex
```