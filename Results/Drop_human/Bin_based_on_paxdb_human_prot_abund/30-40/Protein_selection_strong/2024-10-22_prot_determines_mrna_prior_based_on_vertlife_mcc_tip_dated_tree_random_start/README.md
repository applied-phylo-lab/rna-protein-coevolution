Model was run with bash command
```
Rscript --vanilla R_scripts/runMCMC_OUOU_parallel.R -i Data/Bin_based_on_paxdb_human_prot_abund/10_groups/mean_se_rna_protein_30-40.tsv -o Results/Drop_human/Bin_based_on_paxdb_human_prot_abund/30-40/Protein_selection_strong//2024-10-22_prot_determines_mrna_prior_based_on_vertlife_mcc_tip_dated_tree_random_start/ --burnin_samples 1000 --adapt_samples 25000 --samples 50000 --thin 5 --alpha_star 0.44 --num_cores 4 --drop_species Homo_sapiens --randomize_psi_start_sd 0 --newick --tree Data/vertlife_mcc_dna_only_fossil_bd.tre
```