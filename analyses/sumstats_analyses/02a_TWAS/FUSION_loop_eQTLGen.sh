#!/bin/bash
#$ -N FUSION_for_TRS_project
#$ -pe smp 10
#$ -l h_vmem=10G
#$ -l h_rt=72:0:0
#$ -cwd
#$ -j y
#$ -o /path/to/home/TRS_project/TWAS/FUSION/output/logs/FUSION_for_TRS_project.log

# Change to the FUSION directory
cd /path/to/home/TRS_project/TWAS/FUSION/fusion_twas-master || exit 1

module load R

# Define paths
SUMSTATS_DIR="/path/to/home/TRS_project/sumstats/FUSION_format"
WEIGHTS_DIR="/path/to/home/TRS_project/TWAS/FUSION/fusion_twas-master/weights"
LD_DIR="/path/to/home/TRS_project/TWAS/FUSION/fusion_twas-master/LDREF"
OUTPUT_DIR="/path/to/home/TRS_project/TWAS/FUSION/output"
N_FILE="/path/to/home/TRS_project/sumstats/FUSION_format/NFILE/FUSION_GWASN.csv"

# Loop over all .csv summary statistic files
for SUMSTAT_FILE in "${SUMSTATS_DIR}"/*_no_ambig.txt; do
    SUMSTAT_NAME=$(basename "${SUMSTAT_FILE}" _munged_GRCh37_FUSION_format_no_ambig.txt)

    # Look up GWAS sample size for this file
    GWASN=$(awk -F',' -v file="$(basename "$SUMSTAT_NAME")" 'NR>1 && $1 == file {print $2}' "$N_FILE")

    
    if [ -z "$GWASN" ]; then
        echo "Sample size not found for $(basename "$SUMSTAT_NAME") in $N_FILE. Exiting."
        exit 1
    fi

    for CHR in {1..22}; do
        OUTPUT_FILE="${OUTPUT_DIR}/eQTLGen.eQTL_chr${CHR}_${SUMSTAT_NAME}.fusion"

        if [ -f "${OUTPUT_FILE}" ]; then
            echo "Output file ${OUTPUT_FILE} already exists. Skipping..."
        else
            echo "Running twas for ${SUMSTAT_NAME} with N = ${GWASN}..."

            Rscript FUSION.assoc_test.R \
                --sumstats "${SUMSTAT_FILE}" \
                --weights "${WEIGHTS_DIR}/eQTLGen.eQTL_h2_0.01_cv.performance_updated.pos" \
                --weights_dir "${WEIGHTS_DIR}" \
                --chr "${CHR}" \
                --ref_ld_chr "${LD_DIR}/1000G.EUR." \
                --GWASN "${GWASN}" \
                --coloc_P 0.01 \
                --PANELN "${WEIGHTS_DIR}/eQTLGen_PANEL_N.txt" \
                --out "${OUTPUT_FILE}"
        fi
    done
done





# submit: qsub /path/to/home/TRS_project/sumstats/scripts/FUSION_loop_eQTLGen.sh

## NOTE: For the coloc function you need to install coloc package version 3.2.1 
# remotes::install_version("coloc", version = "3.2-1", repos = "http://cran.us.r-project.org")


