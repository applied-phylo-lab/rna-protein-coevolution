Model was run with bash command
```
Rscript --vanilla R_scripts/runMCMC_OUOU_parallel.R -i Model_testing/Simulated_data/Parameters_same_100_loci/simulated_data_99 -o Model_testing/Results/Test_DIC/Truth_C_eq_0/100_loci/C_neq_0_reverse_trait_fit/Simulated_data_99 --burnin_samples 1000 --adapt_samples 10000 --samples 20000 --no_std_err --reverse_trait_order --sigma_off_diagonal --num_cores 1 --thin 5 --randomize_psi_start_sd 0 --alpha_star 0.44 --tree Data/tree_11sp_noGpig.nex
```