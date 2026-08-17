#!/bin/bash
#$ -N SMR_for_TRS_project
#$ -pe smp 5
#$ -l h_vmem=10G
#$ -l h_rt=24:0:0
#$ -cwd
#$ -j y
#$ -o /path/to/home/TRS_project/SMR/logs/SMR_for_TRS_project.log

# Define paths
BFILE="/path/to/home/TRS_project/ref/g1000_eur/g1000_eur"
SUMSTAT_FOLDER="/path/to/home/TRS_project/sumstats/SMR_format"
EQTLs="/path/to/home/TRS_project/SMR/eQTLs/cis-eQTLs-full_eQTLGen"
OUTPUT_PATH="/path/to/home/TRS_project/SMR/output"


# Loop over all .txt files in the summary statistics folder
for SUMSTAT_FILE in "${SUMSTAT_FOLDER}"/*.csv; do
    SUMSTAT_NAME=$(basename "${SUMSTAT_FILE}" _munged_GRCh37_SMR_format.csv)
    OUTPUT_NAME="${OUTPUT_PATH}/cis-eQTLs-full_eQTLGen_${SUMSTAT_NAME}"

    # Check if output file already exists
    if [ -f "${OUTPUT_NAME}.msmr" ]; then
        echo "Output file ${OUTPUT_NAME}.msmr already exists. Skipping..."
    else
        echo "Processing ${SUMSTAT_FILE}..."
        mkdir -p "${OUTPUT_PATH}"
        /path/to/home/TRS_project/SMR/smr_software/smr \
            --bfile "${BFILE}" \
            --gwas-summary "${SUMSTAT_FILE}" \
            --beqtl-summary "${EQTLs}" \
            --smr-multi \
            --ld-multi-snp 0.1 \
            --peqtl-smr 5e-8 \
            --ld-upper-limit 0.9 \
            --ld-lower-limit 0.05 \
            --peqtl-heidi 1.57e-3 \
            --heidi-min-m 3 \
            --heidi-max-m 20 \
            --maf 0.01 \
            --heidi-mtd 1 \
            --diff-freq 0.2 \
            --out "${OUTPUT_NAME}" \
            --diff-freq-prop 0.05 \
            --thread-num 5
    fi
done

# submit: qsub /path/to/home/TRS_project/sumstats/scripts/SMR_loop_eQTLGen.sh