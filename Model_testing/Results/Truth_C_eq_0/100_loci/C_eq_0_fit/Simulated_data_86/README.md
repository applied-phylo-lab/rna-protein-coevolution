Model was run with bash command
```
Rscript --vanilla R_scripts/runMCMC_OUOU_parallel.R -i Model_testing/Simulated_data/Parameters_diff_wo_c_100_loci_0.5_sdlog/simulated_data_86 -o Model_testing/Results/Truth_C_eq_0/C_eq_0_fit/Simulated_data_86 --burnin_samples 0 --adapt_samples 50000 --samples 50000 --no_std_err --num_cores 1 --thin 5 --randomize_psi_start_sd 0 --alpha_star 0.44 --tree Data/tree_11sp_noGpig.nex
```
