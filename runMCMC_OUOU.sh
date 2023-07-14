#!/bin/bash

INPUT="Data/mean_se_rna_protein.tsv"
#OUTPUT="Results/Real_data_prot_determines_mrna_prior_based_on_tree/"
#nohup Rscript --vanilla R_scripts/runMCMC_OUOU_real_data.R -i ${INPUT} -o ${OUTPUT} --adapt_samples 10000 --samples 50000 --thin 5 --alpha_star 0.44 --tree Data/tree_11sp_noGpig.nex &> Results/real_data_prot_determines_mrna_prior_based_on_tree.out &
#nohup Rscript --vanilla R_scripts/runMCMC_OUOU_real_data.R -i ${INPUT} -o ${OUTPUT} --adapt_samples 0 --samples 50000 --thin 5 --alpha_star 0.44 --prev_rwm_run ${OUTPUT}/rwm_mcmc_1.Rda --rwm_run_number 2 --tree Data/tree_11sp_noGpig.nex &> Results/real_data_prot_determines_mrna_rwm_2.out &


#OUTPUT="Results/Real_data_mrna_determines_prot_prior_based_on_tree/"
#nohup Rscript --vanilla R_scripts/runMCMC_OUOU_real_data.R -i ${INPUT} -o ${OUTPUT} --adapt_samples 10000 --samples 50000 --thin 5 --alpha_star 0.44 --reverse_trait_order --tree Data/tree_11sp_noGpig.nex &> Results/real_data_mrna_determines_prot_prior_based_on_tree.out &
#nohup Rscript --vanilla R_scripts/runMCMC_OUOU_real_data.R -i ${INPUT} -o ${OUTPUT} --adapt_samples 0 --samples 50000 --thin 5 --alpha_star 0.44 --prev_rwm_run ${OUTPUT}/rwm_mcmc_1.Rda --rwm_run_number 2 --reverse_trait_order --tree Data/tree_11sp_noGpig.nex &> Results/real_data_mrna_determines_prot_rwm_2.out &

OUTPUT="Results/Real_data_mrna_determines_prot_off_diagonal_sigma_prior_based_on_tree_run_2/"
nohup Rscript --vanilla R_scripts/runMCMC_OUOU_real_data.R -i ${INPUT} -o ${OUTPUT} --adapt_samples 10000 --samples 50000 --thin 5 --alpha_star 0.44 --sigma_off_diagonal --reverse_trait_order --tree Data/tree_11sp_noGpig.nex &> Results/real_data_mrna_determines_prot_off_diagonal_sigma_prior_based_on_tree_run_2.out &
#nohup Rscript --vanilla R_scripts/runMCMC_OUOU_real_data.R -i ${INPUT} -o ${OUTPUT} --adapt_samples 0 --samples 50000 --thin 5 --alpha_star 0.44 --prev_rwm_run ${OUTPUT}/rwm_mcmc_1.Rda --rwm_run_number 2 --reverse_trait_order --tree Data/tree_11sp_noGpig.nex &> Results/real_data_mrna_determines_prot_rwm_2.out &


#OUTPUT="Results/Real_data_independent_prior_based_on_tree/"
#nohup Rscript --vanilla R_scripts/runMCMC_OUOU_real_data.R -i ${INPUT} -o ${OUTPUT} --adapt_samples 10000 --samples 50000 --thin 5 --alpha_star 0.44 --independent_evo --tree Data/tree_11sp_noGpig.nex &> Results/real_data_independent_prior_based_on_tree.out &
#nohup Rscript --vanilla R_scripts/runMCMC_OUOU_real_data.R -i ${INPUT} -o ${OUTPUT} --adapt_samples 0 --samples 50000 --thin 5 --alpha_star 0.44 --prev_rwm_run ${OUTPUT}/rwm_mcmc_1.Rda --rwm_run_number 2 --independent_evo --tree Data/tree_11sp_noGpig.nex &> Results/real_data_independent_rwm_2.out &



# OUTPUT="Results/Real_data_prot_determines_mrna_omit_root/"
# nohup Rscript --vanilla R_scripts/runMCMC_OUOU_real_data.R -i ${INPUT} -o ${OUTPUT} --adapt_samples 10000 --samples 50000 --thin 5 --alpha_star 0.44 --omit_root --tree Data/tree_11sp_noGpig.nex &> Results/real_data_prot_determines_mrna_omit_root.out &
# 
# OUTPUT="Results/Real_data_mrna_determines_prot_omit_root/"
# nohup Rscript --vanilla R_scripts/runMCMC_OUOU_real_data.R -i ${INPUT} -o ${OUTPUT} --adapt_samples 10000 --samples 50000 --thin 5 --alpha_star 0.44 --omit_root --reverse_trait_order --tree Data/tree_11sp_noGpig.nex &> Results/real_data_mrna_determines_prot_omit_root.out & 
# 
# OUTPUT="Results/Real_data_independent_omit_root/"
# nohup Rscript --vanilla R_scripts/runMCMC_OUOU_real_data.R -i ${INPUT} -o ${OUTPUT} --adapt_samples 10000 --samples 50000 --thin 5 --alpha_star 0.44 --omit_root --independent_evo --tree Data/tree_11sp_noGpig.nex &> Results/real_data_independent_omit_root.out &


# OUTPUT="Results/Real_data_prot_determines_mrna_w_schur_decomp/"
# nohup Rscript --vanilla R_scripts/runMCMC_OUOU_real_data.R -i ${INPUT} -o ${OUTPUT} --adapt_samples 0 --samples 50000 --thin 5 --alpha_star 0.44 --prev_adapt_run ${OUTPUT}/adaptive_mcmc.Rda --tree Data/tree_11sp_noGpig.nex &> Results/real_data_prot_determines_mrna_w_schur_decomp.out &
# 
# OUTPUT="Results/Real_data_mrna_determines_prot_w_schur_decomp/"
# nohup Rscript --vanilla R_scripts/runMCMC_OUOU_real_data.R -i ${INPUT} -o ${OUTPUT} --adapt_samples 0 --samples 50000 --thin 5 --alpha_star 0.44 --prev_adapt_run ${OUTPUT}/adaptive_mcmc.Rda --reverse_trait_order --tree Data/tree_11sp_noGpig.nex &> Results/real_data_mrna_determines_prot_w_schur_decomp.out &
# 
# OUTPUT="Results/Real_data_independent_w_schur_decomp/"
# nohup Rscript --vanilla R_scripts/runMCMC_OUOU_real_data.R -i ${INPUT} -o ${OUTPUT} --adapt_samples 0 --samples 50000 --thin 5 --alpha_star 0.44 --prev_adapt_run ${OUTPUT}/adaptive_mcmc.Rda --independent_evo --tree Data/tree_11sp_noGpig.nex &> Results/real_data_independent_w_schur_decomp.out &

#
# OUTPUT="Results/Real_data_prot_determines_mrna_w_schur_decomp_omit_root/"
# nohup Rscript --vanilla R_scripts/runMCMC_OUOU_real_data.R -i ${INPUT} -o ${OUTPUT} --adapt_samples 10000 --samples 50000 --thin 5 --alpha_star 0.44 --omit_root --tree Data/tree_11sp_noGpig.nex &> Results/real_data_prot_determines_mrna_w_schur_decomp_omit_root.out &
# 
# OUTPUT="Results/Real_data_mrna_determines_prot_w_schur_decomp_omit_root/"
# nohup Rscript --vanilla R_scripts/runMCMC_OUOU_real_data.R -i ${INPUT} -o ${OUTPUT} --adapt_samples 10000 --samples 50000 --thin 5 --alpha_star 0.44 --omit_root --reverse_trait_order --tree Data/tree_11sp_noGpig.nex &> Results/real_data_mrna_determines_prot_w_schur_decomp_omit_root.out & 
# 
# OUTPUT="Results/Real_data_independent_w_schur_decomp_omit_root/"
# nohup Rscript --vanilla R_scripts/runMCMC_OUOU_real_data.R -i ${INPUT} -o ${OUTPUT} --adapt_samples 10000 --samples 50000 --thin 5 --alpha_star 0.44 --omit_root --independent_evo --tree Data/tree_11sp_noGpig.nex &> Results/real_data_independent_w_schur_decomp_omit_root.out &


wait