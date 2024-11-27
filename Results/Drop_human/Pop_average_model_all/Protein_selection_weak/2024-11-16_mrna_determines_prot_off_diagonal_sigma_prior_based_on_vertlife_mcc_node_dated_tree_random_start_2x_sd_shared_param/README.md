Model was run with bash command
```
Rscript --vanilla R_scripts/runMCMC_OUOU_parallel.R -i Data/mean_se_rna_protein.tsv -o Results/Drop_human/Pop_average_model_all//Protein_selection_weak//2024-11-16_mrna_determines_prot_off_diagonal_sigma_prior_based_on_vertlife_mcc_node_dated_tree_random_start_2x_sd_shared_param/ --burnin_samples 0 --adapt_samples 30000 --samples 50000 --thin 5 --alpha_star 0.44 --sigma_off_diagonal --reverse_trait_order --num_cores 8 --h_sd 3 --sigma_sd 3 --q_sd 2 --c_sd 2 --drop_species Homo_sapiens --randomize_psi_start_sd 0 --newick --tree Data/vertlife_mcc_dna_only_node_dated.tre
```
