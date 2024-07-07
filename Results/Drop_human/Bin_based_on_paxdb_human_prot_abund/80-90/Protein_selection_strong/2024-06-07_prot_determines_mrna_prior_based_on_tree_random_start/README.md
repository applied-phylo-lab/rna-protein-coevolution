Model was run with bash command
```
Rscript --vanilla R_scripts/runMCMC_OUOU_parallel.R -i Data/Bin_based_on_paxdb_human_prot_abund/10_groups/mean_se_rna_protein_80-90.tsv -o Results/Drop_human/Bin_based_on_paxdb_human_prot_abund/80-90//Protein_selection_strong//2024-06-07_prot_determines_mrna_prior_based_on_tree_random_start/ --burnin_samples 1000 --adapt_samples 25000 --samples 50000 --thin 5 --alpha_star 0.44 --num_cores 4 --drop_species Homo_sapiens --randomize_psi_start_sd 0 --tree Data/tree_11sp_noGpig.nex
```
