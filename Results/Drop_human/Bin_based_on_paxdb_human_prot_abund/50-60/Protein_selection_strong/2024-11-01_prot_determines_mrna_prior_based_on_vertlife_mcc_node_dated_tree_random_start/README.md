Model was run with bash command
```
Rscript --vanilla R_scripts/runMCMC_OUOU_parallel.R -i Data/Bin_based_on_paxdb_human_prot_abund/10_groups/mean_se_rna_protein_50-60.tsv -o Results/Drop_human/Bin_based_on_paxdb_human_prot_abund/50-60//Protein_selection_strong//2024-11-01_prot_determines_mrna_prior_based_on_vertlife_mcc_node_dated_tree_random_start/ --burnin_samples 0 --adapt_samples 30000 --samples 50000 --thin 5 --alpha_star 0.44 --num_cores 4 --drop_species Homo_sapiens --randomize_psi_start_sd 0 --newick --tree Data/vertlife_mcc_dna_only_node_dated.tre
```
