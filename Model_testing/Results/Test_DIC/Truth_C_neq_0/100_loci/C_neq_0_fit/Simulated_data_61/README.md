Model was run with bash command
```
Rscript --vanilla R_scripts/runMCMC_OUOU_parallel.R -i Model_testing/Simulated_data/Parameters_same_w_c_100_loci/simulated_data_61 -o Model_testing/Results/Test_DIC/Truth_C_neq_0/100_loci/C_neq_0_fit/Simulated_data_61 --burnin_samples 1000 --adapt_samples 10000 --samples 20000 --no_std_err --num_cores 1 --thin 5 --sigma_off_diagonal --randomize_psi_start_sd 0 --alpha_star 0.44 --tree Data/tree_11sp_noGpig.nex
```
