library(data.table)
library(dplyr)
library(stringr)
library(R.utils)
library(TwoSampleMR)
#devtools::install_github("xue-hr/MRcML")
library(MRcML)
library(ggplot2)

harmonised_file <- as.data.frame(fread("/path/to/home/TRS_project/MR/output/harmonised_reverse/PD_LBD_EUR_2021_harmonised.txt")) %>% filter(mr_keep == "TRUE")
harmonised_file$outcome <- "LDB"
# Run methods available in two sample MR
twosamplemr_methods <- mr(harmonised_file, method_list = c("mr_ivw","mr_weighted_median","mr_egger_regression","mr_penalised_weighted_median"))

# Run MRcML and extract the output
cML_output <- MRcML::mr_cML(harmonised_file$beta.exposure,
                                   harmonised_file$beta.outcome,
                                   harmonised_file$se.exposure,
                                   harmonised_file$se.outcome,
                                   n = 65942,
                                   random_start = 100,
                                   random_seed = 1)

cML_results <- data.frame(id.exposure = harmonised_file$id.exposure[1], id.outcome = harmonised_file$id.outcome[1], outcome = "LBD", exposure = "PD", method = "Constrained Maximum Likelihood",nsnp = (nrow(harmonised_file)) - (length(cML_output$BIC_invalid)), b = (cML_output$MA_BIC_theta), se = (cML_output$MA_BIC_se), pval = (cML_output$MA_BIC_p))

# Bind MRcML with the other MR methods and write
twosamplemr_methods <- rbind(twosamplemr_methods,cML_results)
write.table(twosamplemr_methods, file = paste0("/path/to/home/TRS_project/MR/output/IVW_MR_results/PD_LBS_MR_sens_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# Obtain Isq G-X to assess measurement error and suitability for MR-Egger (want > 0.9)
Isq_G_X <- Isq(abs(harmonised_file$beta.exposure),harmonised_file$se.exposure)

# Test for pleioptropy using MR-Egger intercept, add Isq G-X, write file
pleio <- mr_pleiotropy_test(harmonised_file)
pleio$Isq_G_X <- Isq_G_X
write.table(pleio, file = paste0("/path/to/home/TRS_project/MR/output/IVW_MR_results/PD_LBS_MR_pleio_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Test heterogeneity and write
het <- mr_heterogeneity(harmonised_file)
write.table(het, file = paste0("/path/to/home/TRS_project/MR/output/IVW_MR_results/PD_LBS_MR_het_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Perform LOO, write and plot
LOO <- mr_leaveoneout(harmonised_file)
write.table(LOO, file = paste0("/path/to/home/TRS_project/MR/output/IVW_MR_results/PD_LBS_MR_LOO_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

LOO_plot <- mr_leaveoneout_plot(LOO)
ggsave(LOO_plot[[1]], filename = paste0("/path/to/home/TRS_project/MR/output/IVW_MR_results/PD_LBD_MR_leave_plot.pdf"), width = 7, height = 14)

# Run single snp and plot funnel
single_snp <- mr_singlesnp(harmonised_file)

funnel_plot <- mr_funnel_plot(single_snp)
ggsave(funnel_plot[[1]], filename = paste0("/path/to/home/TRS_project/MR/output/IVW_MR_results/PD_LBD_MR_funnel_plot.pdf"), width = 7, height = 7)

# Finally plot scatter 
scatter_plot <- mr_scatter_plot(twosamplemr_methods, harmonised_file)
ggsave(scatter_plot[[1]], filename = paste0("/path/to/home/TRS_project/MR/output/IVW_MR_results/PD_LBD_MR_scatter_plot.pdf"), width = 7, height = 7)


