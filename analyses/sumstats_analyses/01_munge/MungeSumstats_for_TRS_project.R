# Munge Summary Statistics for TRS 
## Clean summary statistics using the MungeSumstats package in R ##

# if (!require("BiocManager")) install.packages("BiocManager")
# BiocManager::install("MungeSumstats")
# BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
# BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
# BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37")
# BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")

library(MungeSumstats)
library(data.table)
library(dplyr)
library(stringr)
library(R.utils)


# To avoid mungesumstats getting confused, we rename so the effect and non-effect are unambiguous

possible_colnames <- c(CHR = "#chrom", CHR = "chromosome", CHR = "chr", CHR = "CHROM", CHR = "#CHR", CHR = "Chr", CHR = "#CHROM", CHR = "Chromsome",
  BP = "pos", BP = "base_pair_location", BP = "Bp", BP = "PosB38", BP = "Position_hg38", BP = "POS", BP = "POS_b37", BP = "base_pair_position", BP = "Position",
  SNP = "snp", SNP = "rsids", SNP = "variant_id", SNP = "marker", SNP = "rsid", SNP = "ID", SNP = "SNPID", SNP = "RSID", SNP = "rsID", SNP = "markername", SNP = "rs_id", SNP = "MarkerName",
  EFFECT_ALLELE = "A1", EFFECT_ALLELE = "alt", EFFECT_ALLELE = "effect_allele", EFFECT_ALLELE = "Allele1", EFFECT_ALLELE = "a1", EFFECT_ALLELE = "EA", EFFECT_ALLELE = "Tested_Allele", EFFECT_ALLELE = "ALLELE1", EFFECT_ALLELE = "Effect_allele", EFFECT_ALLELE = "ALT", EFFECT_ALLELE = "EffectAllele",
  NON_EFFECT_ALLELE = "ref", NON_EFFECT_ALLELE = "A2", NON_EFFECT_ALLELE = "other_allele", NON_EFFECT_ALLELE = "Allele2", NON_EFFECT_ALLELE = "a2", NON_EFFECT_ALLELE = "OA", NON_EFFECT_ALLELE = "Other_Allele", NON_EFFECT_ALLELE = "ALLELE0", NON_EFFECT_ALLELE = "NEA", NON_EFFECT_ALLELE = "Non_Effect_allele", NON_EFFECT_ALLELE = "REF", NON_EFFECT_ALLELE = "NonEffectAllele",
  FRQ = "af_alt", FRQ = "effect_allele_frequency", FRQ = "AF", FRQ = "EAF_A1", FRQ = "all_meta_AF", FRQ = "EAFrq", FRQ = "HRC_FRQ_A1", FRQ = "Freq_Tested_Allele", FRQ = "A1FREQ", FRQ = "Freq1", FRQ = "POOLED_ALT_AF", FRQ = "freq", FRQ = "FREQ", FRQ = "EAF",
  BETA = "beta", BETA = "beta1", BETA = "inv_var_meta_beta", BETA = "Effect", BETA = "EFFECT_SIZE", BETA = "b", BETA = "Beta",
  SE = "sebeta", SE = "standard_error", SE = "se", SE = "inv_var_meta_sebeta", SE = "StdErr",
  Z = "Zscore",
  P = "pval", P = "p_value", P = "pValue", P = "Pval", P = "PVAL", P = "inv_var_meta_p", P = "P_BOLT_LMM", P = "P-value", P = "P.value", P = "pvalue", P = "p",
  N = "Ntotal", N_eff = "Net", N = "n", N = "NMISS", N = "samplesize", N = "NTOT", N = "sample_size",
  INFO = "IMPINFO", INFO = "UK_info", INFO = "info_score")


# Set working directory

setwd("/path/to/home/TRS_project/sumstats/raw/")

data("sumstatsColHeaders")

# Test: files_of_interest <- list.files(pattern="AD_no_UKB_GRCh37.txt")
files_of_interest <- list.files(pattern="GRCh38|GRCh37")

print(files_of_interest)

# test: files_of_interest <- files_of_interest[1]

# Start loop for munging
for(i in files_of_interest){

# Obtain trait name for later and also the cleaned trait name to check if the output file already exists
trait_name <- str_remove(i, pattern = ".gz")
trait_name <- str_remove(trait_name, pattern = ".tsv")
trait_name <- str_remove(trait_name, pattern = ".txt")
trait_name <- str_remove(trait_name, pattern = ".csv")
output_trait_name <- str_remove(trait_name, pattern = "_GRCh37")
output_trait_name <- str_remove(output_trait_name, pattern = "_GRCh38")

output_file_path <- paste0("/path/to/home/TRS_project/sumstats/munged/",output_trait_name,"_munged_GRCh37.csv")

# Check if the file exists, if so the skip to next raw GWAS file for processing
if (file.exists(output_file_path)) {
  print(paste0("Munged file already exists: ",output_file_path))
  next  # Use only inside a loop
} else {
  print(paste0("File not found. Proceeding with processing for ",trait_name))
}


# If there is no munged file, read in raw GWAS file
df<-fread(i, header=T)

# For reference later if needed
print(paste0("There are ",nrow(df)," rows for ",trait_name," before pre-munge filtering.."))

#Rename columns
df <- df %>%
   dplyr::rename(any_of(possible_colnames))

# Capitalise alleles
df$EFFECT_ALLELE <- toupper(df$EFFECT_ALLELE)
df$NON_EFFECT_ALLELE <- toupper(df$NON_EFFECT_ALLELE)

  # Convert to numeric if necessary
  df$BETA <- suppressWarnings(as.numeric(df$BETA))
  df$SE   <- suppressWarnings(as.numeric(df$SE))
  df$P <- suppressWarnings(as.numeric(df$P))

# Filter out non-bialleleic snp
df <- df %>% filter(EFFECT_ALLELE %in% c("A","C","T","G"))
df <- df %>% filter(NON_EFFECT_ALLELE %in% c("A","C","T","G"))


# Filter to retain only MAF ≥ 0.01
if ("FRQ" %in% colnames(df)) {
  df$FRQ <- suppressWarnings(as.numeric(df$FRQ))
  df <- df %>% filter(FRQ >= 0.01 & FRQ <= 0.99)
} 


# Create Z column if BETA and SE are already available and Z is not
if ("BETA" %in% colnames(df) && 
    "SE" %in% colnames(df) && 
    !("Z" %in% colnames(df))) {



  df$Z <- df$BETA / df$SE
}

# Restrict to only autosomes
if ("CHR" %in% colnames(df)){
df <- df %>% filter(CHR != "X")
df <- df %>% filter(CHR != "Y")
df <- df %>% filter(CHR != "XY")
df <- df %>% filter(CHR != "23")
df <- df %>% filter(CHR != "24")
}

# Select columns that are required for downstream analyses, removing redundant columns and reducing files size
cols_of_interest <- c("SNP","CHR","BP","EFFECT_ALLELE","NON_EFFECT_ALLELE","FRQ","BETA","SE","Z","P","N","Neff","INFO")

# Find columns that are both desired and present in the df
available_cols <- intersect(cols_of_interest, names(df))

# Subset
df <- df[, ..available_cols]

print(paste0("There are ",nrow(df)," rows for ",trait_name," after pre-munge filtering..."))

# Check the trait name contains GRCh38, if so convert to GRCh37 (which is used by the majority of downstream software)
if (grepl("GRCh38", trait_name)) {

  trait_name <- str_remove(trait_name, pattern = "_GRCh38")

  format_sumstats(
  path = df,
  ref_genome = "GRCh38", # Set reference genome build for the sumstats
  dbSNP = 155,
  convert_ref_genome = "GRCh37", # Convert to GRCh37
  convert_small_p = TRUE, # For p-values outside of R's range, these are converted to 0
  convert_large_p = TRUE, # For p-value over 1, these are converted to 1
  convert_neg_p = TRUE, # For p-values less that 0 (i.e. negative) these are converted to 0
  compute_z = FALSE, # Whether to convert beta and se to z-score (not required, set as FALSE)
  force_new_z = FALSE, # For if z column already exists (not required, set as FALSE)
  compute_n = 0L, # Compute missing N for SNPs (not required)
  convert_n_int = TRUE, # If N is not an integer, this is rounded
  impute_beta = FALSE, # Impute missing effect value if beta column missing (not required)
  impute_se = FALSE, # Impute missing effect standard error value if column missing (not required) 
  analysis_trait = NULL, # For if multiple traits are contained within the sumstats (not required)
  INFO_filter = 0.7, # If INFO column is present, filter on specified value
  FRQ_filter = 0, # No filtering applied if set to 0. Instead filter in MAF column seperately.
  pos_se = TRUE, # Check all standard errors are positive, if not remove those that arent
  effect_columns_nonzero = TRUE, # Removes SNPs with 0 beta effects 
  N_std = 5, # Remove SNPs with >5 SDs above mean N
  N_dropNA = FALSE, # Drop rows where N is missing (FALSE as not all sumstats will have N here)
  chr_style = 'Ensembl', # Drop chr or CHR prefix if it is present in CHR column (i.e. chr1 etc.)
  on_ref_genome = TRUE, # Checks all SNPs are on the reference genome and imputes if missing
  strand_ambig_filter = FALSE, # Remove strand ambigous SNPs
  allele_flip_check = TRUE, # Check reference allele (A1 here) requires flipping
  allele_flip_drop = TRUE, # If neither allele matches reference genome, they are dropped
  allele_flip_z = TRUE, # Flip the Z score column along with the allele
  allele_flip_frq = TRUE, # Flip the effect allele column
  bi_allelic_filter = FALSE, # Remove non-bialleic SNPs
  flip_frq_as_biallelic = TRUE, # Flip multialleleic snps as though their are bialeleic (should be fine as have already removed EFFECT_ALLELE and NON_EFFECT_ALLELE columns with multiple allele values)
  snp_ids_are_rs_ids = TRUE, # Defines whether SNPs IDs should be inferred as RSIDs
  remove_multi_rs_snp = TRUE, # Remove SNPs with more than 1 RSID (eg. rs123_rs234)
  frq_is_maf = TRUE, # Set as true, stops renaming of column if major allele freq is inferred
  indels = TRUE, # Exclude indels
  sort_coordinates = TRUE, # Sort the coordinates of the resulting summary statistics
  nThread = 1, # Number of threads for analysis
  save_path = paste0("/path/to/home/TRS_project/sumstats/munged/",trait_name,"_munged_GRCh37.csv"),
  write_vcf = FALSE, # Write vcf format file (set to FALSE)
  tabix_index = FALSE, # Related to vcf command above, ignore
  return_data = TRUE, # Return data in a specified format (see next command)
  return_format = "data.table", # Return format as data table
  ldsc_format = FALSE, # If requiring format for use straight away in ldsc
  log_folder_ind = FALSE, # For if log of filtered out snps are required
  log_mungesumstats_msgs = TRUE, # For a log of all messages from the above commands
  log_folder = paste0("/path/to/home/TRS_project/sumstats/munged/munge_logs/",trait_name,"_munged_GRCh37.log"),
  imputation_ind = FALSE, # If columns should be added to the new sumstats specifying filter steps
  force_new = FALSE, # Force a new outout if one already exists (redundant)
  mapping_file = sumstatsColHeaders
)

# else, munge on GRCh37
} else {

  trait_name <- str_remove(trait_name, pattern = "_GRCh37")

 format_sumstats(
  path = df,
  ref_genome = "GRCh37", # Set reference genome build for the sumstats
  dbSNP = 155,
  convert_ref_genome = NULL, # No conversion needed
  convert_small_p = TRUE, # For p-values outside of R's range, these are converted to 0
  convert_large_p = TRUE, # For p-value over 1, these are converted to 1
  convert_neg_p = TRUE, # For p-values less that 0 (i.e. negative) these are converted to 0
  compute_z = FALSE, # Whether to convert beta and se to z-score (not required, set as FALSE)
  force_new_z = FALSE, # For if z column already exists (not required, set as FALSE)
  compute_n = 0L, # Compute missing N for SNPs (not required, set as 0)
  convert_n_int = TRUE, # If N is not an integer, this is rounded
  impute_beta = FALSE, # Impute missing effect value if beta column missing (not required)
  impute_se = FALSE, # Impute missing effect standard error value if column missing (not required) 
  analysis_trait = NULL, # For if multiple traits are contained within the sumstats (not required)
  INFO_filter = 0.7, # If INFO column is present, filter on specified value
  FRQ_filter = 0, # No filtering applied if set to 0. Instead filter in MAF column seperately.
  pos_se = TRUE, # Check all standard errors are positive, if not remove those that arent
  effect_columns_nonzero = TRUE, # Removes SNPs with 0 beta effects 
  N_std = 5, # Remove SNPs with >5 SDs above mean N
  N_dropNA = FALSE, # Drop rows where N is missing (FALSE as not all sumstats will have N here)
  chr_style = 'Ensembl', # Drop chr or CHR prefix if it is present in CHR column (i.e. chr1 etc.)
  on_ref_genome = TRUE, # Checks all SNPs are on the reference genome and imputes if missing
  strand_ambig_filter = FALSE, # Remove strand ambigous SNPs
  allele_flip_check = TRUE, # Check reference allele (A1 here) requires flipping
  allele_flip_drop = TRUE, # If neither allele matches reference genome, they are dropped
  allele_flip_z = TRUE, # Flip the Z score column along with the allele
  allele_flip_frq = TRUE, # Flip the effect allele column
  bi_allelic_filter = FALSE, # Remove non-bialleic SNPs
  flip_frq_as_biallelic = TRUE, # Flip multialleleic snps as though their are bialeleic (should be fine as have already removed EFFECT_ALLELE and NON_EFFECT_ALLELE columns with multiple allele values)
  snp_ids_are_rs_ids = TRUE, # Defines whether SNPs IDs should be inferred as RSIDs
  remove_multi_rs_snp = TRUE, # Remove SNPs with more than 1 RSID (eg. rs123_rs234)
  frq_is_maf = TRUE, # Set as true, stops renaming of column if major allele freq is inferred
  indels = TRUE, # Exclude indels
  sort_coordinates = TRUE, # Sort the coordinates of the resulting summary statistics
  nThread = 1, # Number of threads for analysis
  save_path = paste0("/path/to/home/TRS_project/sumstats/munged/",trait_name,"_munged_GRCh37.csv"),
  write_vcf = FALSE, # Write vcf format file (set to FALSE)
  tabix_index = FALSE, # Related to vcf command above, ignore
  return_data = TRUE, # Return data in a specified format (see next command)
  return_format = "data.table", # Return format as data table
  ldsc_format = FALSE, # If requiring format for use straight away in ldsc
  log_folder_ind = FALSE, # For if log of filtered out snps are required
  log_mungesumstats_msgs = TRUE, # For a log of all messages from the above commands
  log_folder = paste0("/path/to/home/TRS_project/sumstats/munged/munge_logs/",trait_name,"_munged_GRCh37.log"),
  imputation_ind = FALSE, # If columns should be added to the new sumstats specifying filter steps
  force_new = FALSE, # Force a new outout if one already exists (redundant)
  mapping_file = sumstatsColHeaders
)


}

}

