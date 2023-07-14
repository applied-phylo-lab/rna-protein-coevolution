#!/bin/bash

INPUT="Model_testing/Simulated_data/Parameters_diff_50_loci_0.5_sdlog/simulated_data_$1"
OUTPUT="Model_testing/Results/Parameters_diff_50_loci_0.5_sdlog/Simulated_data_$1"
#nohup Rscript --vanilla R_scripts/runMCMC_OUOU_full_model.R -i ${INPUT} -o ${OUTPUT} --adapt_samples 20000 --samples 1 --thin 1 --alpha_star 0.35 --no_std_err --tree Data/tree_11sp_noGpig.nex &> Model_testing/Results/Full_model_parameters_same_50_loci_better_values/simulated_data_$1.out &
nohup Rscript --vanilla R_scripts/runMCMC_OUOU.R -i ${INPUT} -o ${OUTPUT} --adapt_samples 20000 --samples 1 --thin 5 --alpha_star 0.44 --tree Data/tree_11sp_noGpig.nex &> Model_testing/Results/Parameters_diff_50_loci_0.5_sdlog/simulated_data_$1.out &

# OUTPUT="Model_testing/Results/Parameters_diff_reverse_trait_order/Simulated_data_$1"
# nohup Rscript --vanilla R_scripts/runMCMC_OUOU.R -i ${INPUT} -o ${OUTPUT} --adapt_samples 50000 --samples 25000 --thin 5 --reverse_trait_order --alpha_star 0.44 --tree Data/tree_11sp_noGpig.nex &>  Model_testing/Results/Parameters_diff_reverse_trait_order/simulated_data_$1.out & 
wait