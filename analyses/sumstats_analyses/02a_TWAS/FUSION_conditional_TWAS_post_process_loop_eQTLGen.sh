#!/bin/bash
#$ -N FUSION_for_TRS_project_post_process
#$ -pe smp 1
#$ -l h_vmem=30G
#$ -l h_rt=12:0:0
#$ -cwd
#$ -j y
#$ -o /path/to/home/TRS_project/TWAS/FUSION/output/logs/FUSION_for_TRS_project_post_process.log

# Change to the FUSION directory
cd /path/to/home/TRS_project/TWAS/FUSION/fusion_twas-master || exit 1

module load R/4.4.1

# Define paths
SIG_DIR="/path/to/home/TRS_project/TWAS/FUSION/post_process_dir/post_process_input"
SUMSTATS_DIR="/path/to/home/TRS_project/sumstats/FUSION_format"
LD_DIR="/path/to/home/TRS_project/TWAS/FUSION/fusion_twas-master/LDREF"
OUTPUT_DIR="/path/to/home/TRS_project/TWAS/FUSION/post_process_dir/post_process_output"

for SIG_FILE in "${SIG_DIR}"/*.txt; do
    
    # Extract clean name
    SUMSTAT_NAME=$(basename "$SIG_FILE" \
    | sed -E 's/_FDR_fusion_combined_results\.txt$/' \
    | sed -E 's/_COLOC/g')

    TRS_NAME=$(basename "$SIG_FILE" \
    | sed -E 's/_combined_results\.txt$/')

    echo "Processing: $SUMSTAT_NAME"

    for CHR in {1..22}; do
        OUTPUT_FILE="${OUTPUT_DIR}/eQTLGen.eQTL_chr${CHR}_${TRS_NAME}_post_process.fusion"

            Rscript FUSION.post_process.R \
                --sumstats "${SUMSTATS_DIR}/${SUMSTAT_NAME}_munged_GRCh37_FUSION_format_no_ambig.txt" \
                --input "${SIG_FILE}" \
                --chr "${CHR}" \
                --ref_ld_chr "${LD_DIR}/1000G.EUR." \
                --out "${OUTPUT_FILE}" \
                --zthresh 2 \
                --max_r2 0.95 \
                --min_r2 0.1 \
                --locus_win 100000 \
                --verbose 2
    done
done






# submit: qsub /path/to/home/TRS_project/sumstats/scripts/FUSION_post_process_loop_eQTLGen.sh

## NOTE: For the coloc function you need to install coloc package version 3.2.1 
# remotes::install_version("coloc", version = "3.2-1", repos = "http:/cran.us.r-project.org")

# Here we set zthresh to 2 to stop genes being skipped due to low z-score, as genes for scores were selected at FDR significance, as opposed to Bonferroni


