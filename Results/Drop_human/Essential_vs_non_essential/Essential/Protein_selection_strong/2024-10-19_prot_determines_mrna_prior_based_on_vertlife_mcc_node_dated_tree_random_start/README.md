Model was run with bash command
```
Rscript --vanilla R_scripts/runMCMC_OUOU_parallel.R -i Data/Essential/non_essential_mean_se_rna_protein.tsv -o Results/Drop_human/Essential_vs_non_essential/Non_essential//Protein_selection_strong//2024-10-19_prot_determines_mrna_prior_based_on_vertlife_mcc_node_dated_tree_random_start/ --burnin_samples 1000 --adapt_samples 25000 --samples 50000 --thin 5 --alpha_star 0.44 --num_cores 4 --drop_species Homo_sapiens --randomize_psi_start_sd 0 --newick --tree Data/vertlife_mcc_dna_only_node_dated.tre
```
