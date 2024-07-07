Model was run with bash command
```
Rscript --vanilla R_scripts/runMCMC_OUOU_parallel.R -i Data/Essential/essential_mean_se_rna_protein.tsv -o Results/Drop_human/Essential_vs_non_essential/Essential//Protein_selection_strong//2024-06-13_prot_determines_mrna_prior_based_on_tree_random_start/ --burnin_samples 0 --adapt_samples 0 --samples 100000 --thin 5 --alpha_star 0.44 --num_cores 4 --prev_rwm_run Results/Drop_human/Essential_vs_non_essential/Essential//Protein_selection_strong//2024-06-13_prot_determines_mrna_prior_based_on_tree_random_start//rwm_mcmc_1.Rda --rwm_run_number 2 --drop_species Homo_sapiens --randomize_psi_start_sd 0 --tree Data/tree_11sp_noGpig.nex
```
