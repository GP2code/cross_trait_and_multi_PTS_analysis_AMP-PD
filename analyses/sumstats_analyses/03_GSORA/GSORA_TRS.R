#install.packages(c("data.table", "dplyr", "stringr"))
#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))

library(data.table)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

select <- dplyr::select

setwd("/path/to/home/TRS_project/TRS/FDR_scores_TWAS_and_SMR_post_process/TWAS_and_SMR_FDR_post_process_scorefiles")

pass_list <- fread("/path/to/home/TRS_project/gene_set_analysis/input/ALL_FUSION_SMR_FDR_results_meta_analysed_post_process.txt") %>% 
  filter(threshold_pass == "PASS")

pass_list$predictor <- paste0(pass_list$predictor, "_scorefile.txt")

out_dir <- "/path/to/home/TRS_project/gene_set_analysis/output/"

for (i in pass_list$predictor) {
  
  df <- fread(i)
  
  if (nrow(df) < 3) next
  
  pheno_name <- str_remove(i, "_scorefile.txt")
  
  # --- enrichGO (ORA) ---
  ora_GO <- tryCatch(
    enrichGO(
      gene          = df$GENE,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENSEMBL",
      ont           = "ALL",
      pAdjustMethod = "none",
      pvalueCutoff  = 1,
      qvalueCutoff  = 1,
      minGSSize     = 3,
      maxGSSize     = 500,
      readable      = TRUE
    ),
    error = function(e) { message("enrichGO failed for ", pheno_name, ": ", e$message); NULL }
  )
  
  if (!is.null(ora_GO)) {
    tmp <- as.data.frame(ora_GO)
    tmp <- tmp %>% dplyr::select(-qvalue, -p.adjust)
    rownames(tmp) <- NULL
    if (nrow(tmp) > 0) {
      tmp$FDR <- p.adjust(tmp$pvalue, method = "BH")
      tmp$pheno <- pheno_name
      tmp <- tmp %>% filter(pvalue < 0.05)
      print(head(tmp))
      write.table(tmp,
                  file      = paste0(out_dir, pheno_name, "_enrichment_ORA_GO.txt"),
                  sep       = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    }
  }
}


