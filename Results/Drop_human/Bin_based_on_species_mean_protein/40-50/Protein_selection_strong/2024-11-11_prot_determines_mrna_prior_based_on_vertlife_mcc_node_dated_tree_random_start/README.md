Model was run with bash command
```
Rscript --vanilla R_scripts/runMCMC_OUOU_parallel.R -i Data/Bin_based_on_species_mean_protein/10_groups/mean_se_rna_protein_40-50.tsv -o Results/Drop_human/Bin_based_on_species_mean_protein/40-50//Protein_selection_strong//2024-11-11_prot_determines_mrna_prior_based_on_vertlife_mcc_node_dated_tree_random_start/ --burnin_samples 0 --adapt_samples 30000 --samples 50000 --thin 5 --alpha_star 0.44 --num_cores 4 --drop_species Homo_sapiens --randomize_psi_start_sd 0 --newick --tree Data/vertlife_mcc_dna_only_node_dated.tre
```
