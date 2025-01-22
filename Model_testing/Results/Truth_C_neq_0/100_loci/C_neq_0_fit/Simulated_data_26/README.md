Model was run with bash command
```
Rscript --vanilla R_scripts/runMCMC_OUOU_parallel.R -i Model_testing/Simulated_data/Parameters_diff_w_c_100_loci_0.5_sdlog/simulated_data_26 -o Model_testing/Results/Truth_C_neq_0/C_neq_0_fit/Simulated_data_26 --burnin_samples 0 --adapt_samples 0 --samples 50000 --prev_adapt_run Model_testing/Results/Truth_C_neq_0/C_neq_0_fit/Simulated_data_26/adaptive_mcmc.Rda --no_std_err --sigma_off_diagonal --num_cores 1 --thin 5 --randomize_psi_start_sd 0 --alpha_star 0.44 --tree Data/tree_11sp_noGpig.nex
```
