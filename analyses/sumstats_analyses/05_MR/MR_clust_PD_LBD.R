remotes::install_github("cnfoley/mrclust", build_vignettes = FALSE)
install.packages("remotes")
library(mrclust)
library(ggplot2)
library(dplyr)
library(data.table)
#remotes::install_github("phenoscanner/phenoscanner")
install.packages('ieugwasr')
library(ieugwasr)

# Load SNP locations
snp_db <- SNPlocs.Hsapiens.dbSNP155.GRCh38

# In the script we aim to further examine the causal effect of PD on LBD using clustering MR
setwd("/path/to/home/TRS_project/MR/output/IVW_MR_results/")

# Read in the harmonised file
harmonised_file <- as.data.frame(fread("/path/to/home/TRS_project/MR/output/harmonised_reverse/PD_LBD_EUR_2021_harmonised.txt")) %>% filter(mr_keep == "TRUE")

# Asign the effect sizes and se and define ratio estimate for MR-Clust
# Exposure
bx = harmonised_file$beta.exposure
bxse = harmonised_file$se.exposure

# Outcome
by = harmonised_file$beta.outcome
byse = harmonised_file$se.outcome

# Ratio
ratio_est = by/bx
ratio_est_se = byse/abs(bx)

# Instrument names
snp_names = harmonised_file$SNP


# Set seed for MR-Clust
set.seed(1)

# Run MR-Clust
res_em <- mr_clust_em(theta = ratio_est, theta_se = ratio_est_se, bx = bx,
                      by = by, bxse = bxse, byse = byse, obs_names = snp_names)

# Write the results from the best clustering
write.table(res_em$results$best, file = paste0("/path/to/home/TRS_project/MR/output/IVW_MR_results/PD_LBD_MR_clust_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

res_em$results$best
# Extract high quality clusters (variants asigned to clusters with a probability >0.8)
res80 <- mrclust::pr_clust(dta = res_em$results$best, prob = 0.8, min_obs =  4) # There are 11 variants in cluster 1 and 7 in the null cluster
mrclust::pr_clust(dta = res_em$results$best, prob = 0.8, min_obs =  4) 
# Get clustered high prob snps for annotation
harmonised_file_80 <- harmonised_file %>% filter(SNP %in% res80$observation)

harmonised_file_80 <- harmonised_file_80 %>% select(SNP, chr.exposure, pos.exposure, effect_allele.exposure, other_allele.exposure, beta.exposure, se.exposure, pval.exposure, samplesize.exposure)


harmonised_file_80_clust_1 <- harmonised_file_80 %>% filter(SNP %in% res80$observation[res80$cluster_class == 1])

# Write the results from the best clustering
write.table(harmonised_file_80, file = paste0("/path/to/home/TRS_project/MR/output/IVW_MR_results/PD_LBD_MR_clust_prob_80_SNPs.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# Write the results from the best clustering
write.table(harmonised_file_80_clust_1, file = paste0("/path/to/home/TRS_project/MR/output/IVW_MR_results/PD_LBD_MR_clust_prob_80_SNPs_clust_1.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# Plot these high quality clusters
# Asign snp names, betas and se
keep80 = which(snp_names %in% res80$observation)
bx80   = bx[keep80]
bxse80 = bxse[keep80]
by80   = by[keep80]
byse80 = byse[keep80]
snp_names80 = snp_names[keep80]

# Plot with ggplot
plot.sbp.pr80 = two_stage_plot(res = res80, bx = bx80, by = by80, bxse = bxse80,
                               byse = byse80, obs_names = snp_names80) + 
  ggplot2::xlim(0, max(abs(bx80) + 2*bxse80)) + 
  ggplot2::xlab("SNP effect on PD") + 
  ggplot2::ylab("SNP effect on LBD") + 
  ggplot2::ggtitle("")

# Save plot
ggsave(plot.sbp.pr80, filename = paste0("/path/to/home/TRS_project/MR/output/IVW_MR_results/PD_LBD_MR_Clust_plot.pdf"), width = 8, height = 7)

# Perform PHEWAS on the variants in cluster 1 (non-null cluster)
cluster_1 <- res80 %>% filter(cluster_class == 1)
cluster_null <- res80 %>% filter(cluster_class == 'Null')
background_1 <- harmonised_file %>% filter(!(SNP %in% cluster_1$observation), mr_keep == "TRUE")
background_null <- harmonised_file %>% filter(!(SNP %in% cluster_null$observation), mr_keep == "TRUE")

# Instead of "YOUR_API_TOKEN add your OpenGWAS API"
phewas_res <- phewas(variants = cluster_1$observation, pval = 1e-5, opengwas_jwt = "YOUR_API_TOKEN")
phewas_res_null <- phewas(variants = cluster_null$observation, pval = 1e-5, opengwas_jwt = "YOUR_API_TOKEN")
phewas_background_1 <- phewas(variants = background_1$SNP, pval = 1e-5, opengwas_jwt = "YOUR_API_TOKEN")
phewas_res_null <- phewas(variants = background_null$SNP, pval = 1e-5, opengwas_jwt = "YOUR_API_TOKEN")

# Write the results from the best clustering
write.table(phewas_res, file = paste0("/path/to/home/TRS_project/MR/output/IVW_MR_results/PD_LBD_MR_clust_PHEWAS_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(phewas_res_null, file = paste0("/path/to/home/TRS_project/MR/output/IVW_MR_results/PD_LBD_MR_clust_null_cluster_PHEWAS_results.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(phewas_background_1, file = paste0("/path/to/home/TRS_project/MR/output/IVW_MR_results/PD_LBD_MR_clust_PHEWAS_results_background.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(phewas_res_null, file = paste0("/path/to/home/TRS_project/MR/output/IVW_MR_results/PD_LBD_MR_clust_null_cluster_PHEWAS_results_background.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Create some files for ANNOVAR
cluster_1_anno <- harmonised_file %>% 
  filter(SNP %in% cluster_1$observation)

# Duplicate the position column manually
cluster_1_anno <- cluster_1_anno %>%
  mutate(pos.exposure_2 = pos.exposure) %>%
  dplyr::select(chr.exposure, pos.exposure, pos.exposure_2, other_allele.exposure, effect_allele.exposure, SNP)

write.table(cluster_1_anno, file = paste0("/path/to/home/TRS_project/MR/annotation_follow_up/input/PD_LBD_MR_clust_cluster_1_snps.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


cluster_null_anno <- harmonised_file %>% filter(SNP %in% cluster_null$observation)
# Duplicate the position column manually
cluster_null_anno <- cluster_null_anno %>%
  mutate(pos.exposure_2 = pos.exposure) %>%
  dplyr::select(chr.exposure, pos.exposure, pos.exposure_2, other_allele.exposure, effect_allele.exposure, SNP)

write.table(cluster_null_anno, file = paste0("/path/to/home/TRS_project/MR/annotation_follow_up/input/PD_LBD_MR_clust_cluster_null_snps.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

nrow(background_1)

