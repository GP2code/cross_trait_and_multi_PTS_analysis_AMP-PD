library(data.table)
library(stringr)
library(dplyr)
library(ieugwasr)
#remotes::install_github("MRCIEU/genetics.binaRies")

setwd("/path/to/home/TRS_project/SMR/post_process/post_process_input")

# Output directory (make sure it exists)
output_dir_main <- "/path/to/home/TRS_project/SMR/post_process/post_process_output"


# List all .fusion files
files <- list.files(pattern = ".txt")

print(files)

for (i in files){
  df <- fread(i)
  
  name <- str_remove(i, "_results.txt")
  
  
  if (grepl("single_SNP", i)){
    df_smallest_p <- df %>%
      group_by(topSNP) %>%
      slice_min(order_by = p_SMR, n = 1)
    
    p_col <- "p_SMR"
    
  }else{
    
    df_smallest_p <- df %>%
      group_by(topSNP) %>%
      slice_min(order_by = p_SMR_multi, n = 1)
    
    p_col <- "p_SMR_multi"
  }
  
  
  SMR_results_clumped <- df_smallest_p %>%
    rename(rsid = topSNP,
           pval = !!sym(p_col)) %>%
    ieugwasr::ld_clump(clump_r2 = 0.1,
                       clump_p = 1,
                       clump_kb = 250,
                       plink_bin = genetics.binaRies::get_plink_binary(), 
                       bfile = "/path/to/home/TRS_project/ref/g1000_eur/g1000_eur")
  
  
  
  SMR_results_clumped <- as.data.table(SMR_results_clumped)
  
  print(SMR_results_clumped)
  
  sink("/path/to/home/TRS_project/SMR/post_process/post_process_output/n_genes.txt", append = TRUE)
  print(paste0("The number of rows for ",name," is ",nrow(df)))
  print(paste0("The number of rows for ",name," after retaining top gene for each SNP is ",nrow(df_smallest_p)))
  print(paste0("The number of rows for ",name," after retaining clumping SNPs is ",nrow(SMR_results_clumped)))
  print("Writing clumped SMR results...")
  sink()
  
  # Write outputs
  write.table(SMR_results_clumped, file = file.path(output_dir_main, paste0(name, "_post_process_clumped.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
}


