library(remotes)
library(data.table)
library(dplyr)
#library(tidyr)
library(stringr)
library(R.utils)
library(ieugwasr)
#remotes::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
#remotes::install_github("CBIIT/LDlinkR")
library(LDlinkR)

## This script extracts proxies where needed, harmonises them and performs IVW-MR for the effect of the exposure on PD
# We can define a function to get proxies from LDLinkR based on code from Shea Andrews (here: https://andrewslabucsf.github.io/MR-tutorial/scripts/mr_harmonization.html)

# First we define the munge_proxies function for later use
munge_proxies <- function(LDLink_file, outcome, outcome_clump){
  LDLink_file_path <- LDLink_file
  proxy_snps <- readr::read_tsv(LDLink_file_path, skip = 1, col_names = F) %>%
    rename(id = X1, func = X2, proxy_snp = X3, coord = X4, alleles = X5, maf = X6, 
           distance = X7, dprime = X8, rsq = X9, correlated_alleles = X10, FORGEdb = X11, RegulomeDB = X12) %>%
    tidyr::separate(coord, c('chr', 'pos'), sep = ":") %>%
    mutate(snp = ifelse(id == 1, proxy_snp, NA), 
           chr = str_replace(chr, 'chr', ""), 
           chr = as.numeric(chr), 
           pos = as.numeric(pos)) %>%
    tidyr::fill(snp, .direction = 'down') %>%
    relocate(snp, .before = proxy_snp) %>%
    dplyr::select(-id, -func, -FORGEdb, -RegulomeDB) %>%
    filter(rsq >= 0.8)
  
  # Munge proxy snp and outcome data
  proxy_outcome <- left_join(
    proxy_snps, outcome, by = c("proxy_snp" = "SNP")
  ) %>%
    tidyr::separate(correlated_alleles, c("target_a1.outcome", "proxy_a1.outcome", 
                                          "target_a2.outcome", "proxy_a2.outcome"), sep = ",|=") %>%
    filter(!is.na(chr.outcome)) %>%
    arrange(snp, -rsq, abs(distance)) %>%
    group_by(snp) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      proxy.outcome = TRUE,
      target_snp.outcome = snp,
      proxy_snp.outcome = proxy_snp, 
    ) %>% 
    mutate(
      new_effect_allele.outcome = case_when(
        proxy_a1.outcome == effect_allele.outcome & proxy_a2.outcome == other_allele.outcome ~ target_a1.outcome,
        proxy_a2.outcome == effect_allele.outcome & proxy_a1.outcome == other_allele.outcome ~ target_a2.outcome,
        TRUE ~ NA_character_
      ), 
      new_other_allele.outcome = case_when(
        proxy_a1.outcome == effect_allele.outcome & proxy_a2.outcome == other_allele.outcome ~ target_a2.outcome,
        proxy_a2.outcome == effect_allele.outcome & proxy_a1.outcome == other_allele.outcome ~ target_a1.outcome,
        TRUE ~ NA_character_
      ), 
      effect_allele.outcome = new_effect_allele.outcome, 
      other_allele.outcome = new_other_allele.outcome
    ) %>%
    dplyr::select(-proxy_snp, -chr, -pos, -alleles, -maf, -distance, -rsq, -dprime,  
                  -new_effect_allele.outcome, -new_other_allele.outcome) %>%
    relocate(target_a1.outcome, proxy_a1.outcome, target_a2.outcome, proxy_a2.outcome, .after = proxy_snp.outcome) %>%
    rename(SNP = snp) %>%
    relocate(SNP, .after = samplesize.outcome)
  
  # Merge outcome and proxy outcomes
  outcome_dat <- bind_rows(
    outcome_clump, proxy_outcome
  ) %>% 
    arrange(chr.outcome, pos.outcome)
  
  outcome_dat
}

# Set the working directory
setwd("/path/to/home/TRS_project/MR/instruments/clumped_instruments")

# Read in PD GWAS and set as data frame and give phenotype column
PD_GWAS <- fread("/path/to/home/TRS_project/MR/input/outcome/Parkinsons_disease_GP2_CLINICAL_ONLY_EUR_2025_alleles_checked_1KG.txt")
PD_GWAS <- as.data.frame(PD_GWAS)
PD_GWAS$Phenotype <- "PD"

# We can format the whole PD GWAS now for later as it is always the outcome
PD_GWAS_formatted <- format_data(
  PD_GWAS,
  type = "outcome",
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "FRQ",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  samplesize_col = "N",
  chr_col = "CHR",
  pos_col = "BP",
  min_pval = 1e-300,
  log_pval = FALSE
)

# Define a list of the files that contain SNPs that need proxies
need_proxy <- list.files(pattern = "_needs_proxies.txt")
print(need_proxy)
#TEST:
#need_proxy <- need_proxy [1]

mr_results_all <- data.frame()

for (i in need_proxy) {
  
  need_proxy_instruments <- fread(i, header = FALSE)
  exposure_name <- str_remove(i, "_instruments_list_needs_proxies.txt")
  
  if (nrow(need_proxy_instruments) == 0) {
    
    instruments <- fread(paste0(exposure_name, "_instruments_list.txt"), header = F)
    
    exposure_file <- fread(
      paste0("/path/to/home/TRS_project/MR/input/exposures/",
             exposure_name, "_alleles_checked_1KG.txt")
    ) %>%
      filter(snp %in% instruments$V1)
    
    exposure_file <- as.data.frame(exposure_file)
    
    exposure_file$Phenotype <- exposure_name
    
    exposure_formatted <- format_data(
      exposure_file,
      type = "exposure",
      header = TRUE,
      snp_col = "snp",
      beta_col = "b",
      se_col = "se",
      eaf_col = "freq",
      effect_allele_col = "A1",
      other_allele_col = "A2",
      pval_col = "p",
      samplesize_col = "n",
      min_pval = 1e-300,
      log_pval = FALSE
    )
    
    F_stats <- exposure_formatted$beta.exposure^2/exposure_formatted$se.exposure^2
    min_F_stat <- min(F_stats)
    max_F_stat <- max(F_stats)
    mean_F_stat <- mean(F_stats)
    
    
    print(paste0("The minimum F statistic is ",min_F_stat))
    print(paste0("The maximum F statistic is ",max_F_stat))
    print(paste0("The mean F statistic is ",mean_F_stat))
    
    
    
    outcome_data <- format_data(
      PD_GWAS,
      type = "outcome",
      snps = exposure_file$snp,
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      eaf_col = "FRQ",
      effect_allele_col = "A1",
      other_allele_col = "A2",
      pval_col = "P",
      samplesize_col = "N",
      chr_col = "CHR",
      pos_col = "BP",
      min_pval = 1e-300,
      log_pval = FALSE
    )
    
    harmonised_data <- harmonise_data(
      exposure_formatted,
      outcome_data,
      action = 2
    )
    
    write.table(harmonised_data, file = paste0("/path/to/home/TRS_project/MR/output/harmonised/",exposure_name,"_PD_harmonised.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    
    mr_results <- mr(harmonised_data, method_list = "mr_ivw")
    
    
    mr_results$min_F_stat <- min_F_stat
    mr_results$max_F_stat <- max_F_stat
    mr_results$mean_F_stat <- mean_F_stat
    
    print(mr_results)
    
    mr_results_all <- rbind(mr_results_all,mr_results)
    
  }else{
    
    instruments <- fread(paste0(exposure_name, "_instruments_list.txt"), header = F)
    
    exposure_file <- fread(
      paste0("/path/to/home/TRS_project/MR/input/exposures/",
             exposure_name, "_alleles_checked_1KG.txt")
    ) %>%
      filter(snp %in% instruments$V1)
    
    exposure_file <- as.data.frame(exposure_file)
    
    exposure_file$Phenotype <- exposure_name
    
    exposure_formatted <- format_data(
      exposure_file,
      type = "exposure",
      header = TRUE,
      snp_col = "snp",
      beta_col = "b",
      se_col = "se",
      eaf_col = "freq",
      effect_allele_col = "A1",
      other_allele_col = "A2",
      pval_col = "p",
      samplesize_col = "n",
      min_pval = 1e-300,
      log_pval = FALSE
    )
    
    F_stats <- exposure_formatted$beta.exposure^2/exposure_formatted$se.exposure^2
    min_F_stat <- min(F_stats)
    max_F_stat <- max(F_stats)
    mean_F_stat <- mean(F_stats)
    
    
    print(paste0("The minimum F statistic is ",min_F_stat))
    print(paste0("The maximum F statistic is ",max_F_stat))
    print(paste0("The mean F statistic is ",mean_F_stat))
    
    
    
    matched_snps_list <- fread(paste0(exposure_name, "_instruments_list_matched.txt"), header = F)
    matched_snps <- PD_GWAS_formatted %>% filter(SNP %in% matched_snps_list$V1)
    
    LDproxy_batch(need_proxy_instruments$V1, 
                  pop = "EUR",             # Match population ancestries
                  r2d = "r2", 
                  token = 'f3d93542a719', 
                  append = TRUE,           # We appended the results of each LDlink query to a single file
                  genome_build = "grch37") # Select genome build based on summary stats
    
    system("mv combined_query_snp_list_grch37.txt exposure_outcome_proxy_snps.txt")
    
    outcome_data_w_proxies <- munge_proxies("exposure_outcome_proxy_snps.txt", PD_GWAS_formatted, matched_snps)
    
    harmonised_data <- harmonise_data(
      exposure_formatted,
      outcome_data_w_proxies,
      action = 2
    )
    
    write.table(harmonised_data, file = paste0("/path/to/home/TRS_project/MR/output/harmonised/",exposure_name,"_PD_harmonised.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    
    
    mr_results <- mr(harmonised_data, method_list = "mr_ivw")
    
    mr_results$min_F_stat <- min_F_stat
    mr_results$max_F_stat <- max_F_stat
    mr_results$mean_F_stat <- mean_F_stat
    
    print(mr_results)
    
    mr_results_all <- rbind(mr_results_all,mr_results)
    
    
  }
  
}

mr_results_all <- mr_results_all %>% arrange(pval)

write.table(mr_results_all, file = paste0("/path/to/home/TRS_project/MR/output/IVW_MR_results/all_sig_TRS_exposures_PD_IVW.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



