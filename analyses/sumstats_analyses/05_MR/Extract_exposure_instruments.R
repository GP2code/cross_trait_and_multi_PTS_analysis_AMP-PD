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
setwd("/path/to/home/TRS_project/MR/input/exposures/")

ref_snps <- fread("/path/to/home/TRS_project/MR/ref/EUR.bim")

PD_GWAS <- fread("/path/to/home/TRS_project/MR/input/outcome/Parkinsons_disease_GP2_CLINICAL_ONLY_EUR_2025_munged_GRCh37.csv")

# Extract the matching snps in order to obtain proxies later
  PD_GWAS <- PD_GWAS %>% filter(SNP %in% ref_snps$V2)

  # Merge to get reference alleles
  PD_GWAS <- PD_GWAS %>%
  left_join(ref_snps, by = c("SNP" = "V2"))

  # Only keep SNPs with matching alleles (either strand)
  PD_GWAS <- PD_GWAS %>% 
    filter(
      (A1 == V5 & A2 == V6) | 
      (A1 == V6 & A2 == V5)
    )

  write.table(PD_GWAS, file = paste0("/path/to/home/TRS_project/MR/input/outcome/Parkinsons_disease_GP2_CLINICAL_ONLY_EUR_2025_alleles_checked_1KG.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    
  PD_SNPs <- PD_GWAS %>% select(SNP)

  write.table(PD_SNPs, file = paste0("/path/to/home/TRS_project/MR/input/outcome/Parkinsons_disease_GP2_CLINICAL_ONLY_EUR_2025_alleles_checked_1KG_SNP_list.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    


# List the exposure files
exposure_files <- list.files(pattern="_munged_GRCh37_SMR_format.csv")
#exposure_files <- exposure_files[1]
print(exposure_files) 

sink("/path/to/home/TRS_project/MR/instruments/clumped_instruments/MR_instruments_log.txt", append = TRUE)
for(i in exposure_files){
  # read in GWAS files
  print("===============================================")
  print("===============================================")
  exposure_file<-fread(i)
  
  exposure_name <- str_remove(i, pattern = "_munged_GRCh37_SMR_format.csv")
  
  
  # Extract the matching snps in order to obtain proxies later
  exposure_file <- exposure_file %>% filter(snp %in% ref_snps$V2)

  # Merge to get reference alleles
  exposure_file <- exposure_file %>%
  left_join(ref_snps, by = c("snp" = "V2"))

  # Only keep SNPs with matching alleles (either strand)
  exposure_file <- exposure_file %>% 
    filter(
      (A1 == V5 & A2 == V6) | 
      (A1 == V6 & A2 == V5)
    )
  
  write.table(exposure_file, file = paste0("/path/to/home/TRS_project/MR/input/exposures/",exposure_name,"_alleles_checked_1KG.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        
  exposure_snps <- exposure_file %>% select(snp)
  
  
  print(paste0("Running for ...",exposure_name))
  print(paste0("The number of SNPs on in the GWAS also in the reference panel before joining is ",nrow(exposure_snps)))
  
  matched_snps <- exposure_snps %>% filter(snp %in% PD_SNPs$SNP)
  
  print(paste0("The number of rows after matching with the PD GWAS is ",nrow(matched_snps)))
  
  print(paste0("Writing matched SNPs..."))
  
  write.table(matched_snps, file = paste0("/path/to/home/TRS_project/MR/instruments/matched_snps/",exposure_name,"_PD_SNPs.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # We can reduce the size of the exposure file for clumping
  exposure_file <- exposure_file %>% filter(p<=5e-6)
  
  exposure_file$Phenotype <- exposure_name
  
  exposure_file <- as.data.frame(exposure_file)
  
  # Next we need to extract the instruments
  print(paste0("Extracting IVs for ",exposure_name))
  
  
  # Format using twosampleMR
  exposure_data <- format_data(
    exposure_file,
    type = "exposure",
    snps = NULL,
    header = TRUE,
    phenotype_col = "Phenotype",
    snp_col = "snp",
    beta_col = "b",
    se_col = "se",
    eaf_col = "freq",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "p",
    samplesize_col = "n",
    min_pval = 1e-300,
    log_pval = FALSE)
  
  
  # Clump the exposure data
  if (min(exposure_data$pval.exposure)<= 5e-8){
    exposure_data_clumped <- exposure_data %>%
      rename(rsid = SNP,
             pval = pval.exposure) %>%
      ieugwasr::ld_clump(clump_r2 = 0.001,
                         clump_p = 5e-8,
                         clump_kb = 10000,
                         plink_bin = genetics.binaRies::get_plink_binary(), 
                         bfile = "/path/to/home/TRS_project/MR/ref/EUR")
    
    print(paste0(exposure_name," has ",nrow(exposure_data_clumped)," IVs at 5e-8, proceeding..."))
    
  }else{
    
    print(paste0(exposure_name," has no possible IVs at 5e-8, using 5e-6..."))

    exposure_data_clumped <- exposure_data %>%
      rename(rsid = SNP,
             pval = pval.exposure) %>%
      ieugwasr::ld_clump(clump_r2 = 0.001,
                         clump_p = 5e-6,
                         clump_kb = 10000,
                         plink_bin = genetics.binaRies::get_plink_binary(), 
                         bfile = "/path/to/home/TRS_project/MR/ref/EUR")
    }
  
  
  exposure_data_clumped <- exposure_data_clumped %>% rename(SNP = rsid, pval.exposure = pval)
  
  
  #print(nrow(exposure_data_clumped))
  
  # Check for sufficient insturments
  if (nrow(exposure_data_clumped)< 5){
    
    print(paste0(exposure_name," has only ",nrow(exposure_data_clumped)," IVs at 5e-8, using 5e-6"))

    
    exposure_data_clumped <- exposure_data %>%
      rename(rsid = SNP,
             pval = pval.exposure) %>%
      ieugwasr::ld_clump(clump_r2 = 0.001,
                         clump_p = 5e-6,
                         clump_kb = 10000,
                         plink_bin = genetics.binaRies::get_plink_binary(), 
                         bfile = "/path/to/home/TRS_project/MR/ref/EUR")
    
    
    exposure_data_clumped <- exposure_data_clumped %>% rename(SNP = rsid, pval.exposure = pval)
    
  }
  
  
        print(paste0("The total number of instruments selected for ",exposure_name," is ",nrow(exposure_data_clumped)))
        
        
        ## We may need to get proxies for any missing instruments
        exposure_data_clumped <- exposure_data_clumped %>% select(SNP)
        
        print(paste0("Writing full instruments list for ",exposure_name))
        
        write.table(exposure_data_clumped, file = paste0("/path/to/home/TRS_project/MR/instruments/clumped_instruments/",exposure_name,"_instruments_list.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


        
        exposure_data_matched <- exposure_data_clumped %>% filter(SNP %in% PD_SNPs$SNP)

        print(paste0("The number of instruments matched with the PD GWAS for ",exposure_name," is ",nrow(exposure_data_matched)))

        print(paste0("Writing matched instruments list for ",exposure_name))

        write.table(exposure_data_matched, file = paste0("/path/to/home/TRS_project/MR/instruments/clumped_instruments/",exposure_name,"_instruments_list_matched.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


        exposure_data_clumped_needs_proxies <- exposure_data_clumped %>% filter(!(SNP %in% PD_SNPs$SNP))

        print(paste0("The number of instruments that need proxies for ",exposure_name," is ",nrow(exposure_data_clumped_needs_proxies)))  

        print(paste0("Writing instruments that need proxies for ",exposure_name))

        write.table(exposure_data_clumped_needs_proxies, file = paste0("/path/to/home/TRS_project/MR/instruments/clumped_instruments/",exposure_name,"_instruments_list_needs_proxies.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

        
          print("===============================================")
          print("===============================================")
}

sink()
    

