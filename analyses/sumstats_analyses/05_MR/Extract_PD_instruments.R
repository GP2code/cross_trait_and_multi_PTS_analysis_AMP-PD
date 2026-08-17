library(remotes)
library(data.table)
library(dplyr)
#library(tidyr)
library(stringr)
library(R.utils)
library(ieugwasr)
#remotes::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)

# This script creates a a list of instruments in the risk factors and a list of snps that are shared between the PD outcome GWAS and the risk factor
# We can use these output together in extract proxy instruments if they are missing from the outcome dataset

#Set working directory
setwd("/path/to/home/TRS_project/MR/instruments/matched_snps/")

PD_GWAS <- fread("/path/to/home/TRS_project/MR/input/outcome/Parkinsons_disease_GP2_CLINICAL_ONLY_EUR_2025_alleles_checked_1KG.txt")

PD_SNPs <- fread("/path/to/home/TRS_project/MR/input/outcome/Parkinsons_disease_GP2_CLINICAL_ONLY_EUR_2025_alleles_checked_1KG_SNP_list.txt")

exposure_name <- "PD"

PD_GWAS$Phenotype <- exposure_name

PD_GWAS <- as.data.frame(PD_GWAS)

# Format using twosampleMR
PD_GWAS_formatted <- format_data(
  PD_GWAS,
  type = "exposure",
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

if (min(PD_GWAS_formatted$pval.exposure)<= 5e-8){
  PD_GWAS_clumped <- PD_GWAS_formatted %>%
    rename(rsid = SNP,
           pval = pval.exposure) %>%
    ieugwasr::ld_clump(clump_r2 = 0.001,
                       clump_p = 5e-8,
                       clump_kb = 10000,
                       plink_bin = genetics.binaRies::get_plink_binary(), 
                       bfile = "/path/to/home/TRS_project/MR/ref/EUR")}
  
  
  PD_GWAS_clumped <- PD_GWAS_clumped %>% rename(SNP = rsid, pval.exposure = pval)
  
  PD_GWAS_clumped <- PD_GWAS_clumped %>% select(SNP)
  sink("/path/to/home/TRS_project/MR/instruments/clumped_instruments_reverse/PD_MR_instruments_log.txt", append = TRUE)
  print("===============================================")
  print("===============================================")
  
  print(paste0("There are ",nrow(PD_GWAS_clumped)," total instruments for ",exposure_name))
  print(paste0("Writing full instruments list for ",exposure_name))
  
  write.table(PD_GWAS_clumped, file = paste0("/path/to/home/TRS_project/MR/instruments/clumped_instruments_reverse/",exposure_name,"_instruments_list.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  outcome_matched_snps <- list.files(pattern="_PD_SNPs.txt")
  # TEST: outcome_matched_snps <- outcome_matched_snps[1]
  #print(outcome_matched_snps) 
  
  for(i in outcome_matched_snps){
    # read in GWAS files
    outcome_matched_snps_file <-fread(i, header = F)
    
    outcome_name <- str_remove(i, pattern = "_PD_SNPs.txt")
    
    
    PD_GWAS_matched <- PD_GWAS_clumped %>% filter(SNP %in% outcome_matched_snps_file$V1)
    
    print(paste0("The number of PD instruments matched with ",outcome_name," is ",nrow(PD_GWAS_matched)))
    
    print(paste0("Writing matched instruments list for ",outcome_name))
    
    write.table(PD_GWAS_matched, file = paste0("/path/to/home/TRS_project/MR/instruments/clumped_instruments_reverse/PD_",outcome_name,"_instruments_list_matched.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    
    PD_GWAS_clumped_needs_proxies <- PD_GWAS_clumped %>% filter(!(SNP %in% outcome_matched_snps_file$V1))
    
    print(paste0("The number of PD instruments that need proxies for ",outcome_name," is ",nrow(PD_GWAS_clumped_needs_proxies)))  
    
    print(paste0("Writing instruments that need proxies for ",outcome_name))
    
    write.table(PD_GWAS_clumped_needs_proxies, file = paste0("/path/to/home/TRS_project/MR/instruments/clumped_instruments/PD_",outcome_name,"_instruments_list_needs_proxies.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    
    print("===============================================")
    print("===============================================")
  }
  
  sink()
    

