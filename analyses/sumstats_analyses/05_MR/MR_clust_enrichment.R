library(data.table)
library(dplyr)

# Load data

df_1 <- fread("/path/to/home/AMP-PD/results/MR/PD_LBD_MR_clust_null_cluster_PHEWAS_results.txt")
df_1$cluster <- 1

nrow(df_1)

df_background <- fread("/path/to/home/AMP-PD/results/MR/PD_LBD_MR_clust_null_cluster_PHEWAS_results_background.txt")
df_background$cluster <- 0

harmonised_file <- fread("/path/to/home/AMP-PD/scripts/6.MR/PD_LBD_EUR_2021_harmonised.txt") %>% filter(mr_keep == "TRUE") %>% select(SNP)

df <- rbind(df_1, df_background)

# =========================================================
# 2. CLEAN (unique SNP–trait–cluster associations)
# =========================================================

df <- df %>%
  distinct(rsid, trait, cluster)

# =========================================================
# 3. DEFINE UNIVERSES
# =========================================================

all_snps <- 
  unique(harmonised_file$SNP)

cluster_snps <- df %>%
  filter(cluster == 1) %>%
  pull(rsid) %>%
  unique()

background_snps <- setdiff(all_snps, cluster_snps)

n_cluster <- length(cluster_snps)
n_background <- length(background_snps)

# =========================================================
# 4. PER-TRAIT ENRICHMENT
# =========================================================

traits <- unique(df_1$trait)

results <- lapply(traits, function(tr) {
  
  df_tr <- df %>% filter(trait == tr)
  
  trait_snps <- unique(df_tr$rsid)
  
  # -----------------------------
  # counts
  # -----------------------------
  a <- length(intersect(trait_snps, cluster_snps))      # cluster hits
  c <- length(intersect(trait_snps, background_snps))   # background hits
  
  b <- n_cluster - a
  d <- n_background - c
  
  # -----------------------------
  # proportions (Table 4 style)
  # -----------------------------
  tp_rate <- a / n_cluster
  fp_rate <- c / n_background
  
  # -----------------------------
  # contingency table
  # -----------------------------
  mat <- matrix(
    c(a, b,
      c, d),
    nrow = 2,
    byrow = TRUE
  )
  
  # Fisher exact test (enrichment in cluster)
  ft <- fisher.test(mat, alternative = "greater")
  
  data.frame(
    trait = tr,
    
    true_positives = paste0(a, "/", n_cluster),
    false_positives = paste0(c, "/", n_background),
    
    tp_rate = tp_rate,
    fp_rate = fp_rate,
    
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value
  )
})

res_df <- bind_rows(results)


res_df <- res_df %>%
  arrange(p_value)

# =========================================================
# 7. OUTPUT
# =========================================================

res_df

write.table(res_df,"/path/to/home/AMP-PD/results/MR/MR_clust_cluster_null_pheWAS_enrichment.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
