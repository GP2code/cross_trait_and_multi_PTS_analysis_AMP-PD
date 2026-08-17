#!/bin/bash
#$ -N SBayesRC_calc_hapmap_sumstats_array
#$ -t 1-18
#$ -pe smp 1
#$ -l h_vmem=50G
#$ -l h_rt=6:0:0
#$ -cwd
#$ -j y
#$ -o /path/to/home/TRS_project/PGS/scorefiles/SBayesRC_calc_hapmap_sumstats_array.$TASK_ID.log

set -euo pipefail

# Define directories
MAIN_DIR="/path/to/home/TRS_project"
SUMSTATS_DIR="${MAIN_DIR}/PGS/sumstats/imputed"
LD_DIR="/path/to/home/omic_score_project/PGS/ref"
OUTPUT_DIR="${MAIN_DIR}/PGS/scorefiles"
GCTB_BINARY="/path/to/home/TRS_project/general_software/gctb"
ANNO_DIR="/path/to/home/omic_score_project/PGS/annotations"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Create array of all .txt files
mapfile -t SUMSTAT_FILES < <(find "${SUMSTATS_DIR}" -maxdepth 1 -name "*_hapmap3.imputed.ma" | sort)

# Check if files exist
if [ ${#SUMSTAT_FILES[@]} -eq 0 ]; then
    echo "ERROR: No imputed _hapmap3.imputed.ma files found in $SUMSTATS_DIR" >&2
    exit 1
fi

# Get the current file based on array task ID (SGE_TASK_ID is 1-indexed)
ARRAY_INDEX=$((SGE_TASK_ID - 1))

# Check if task ID is within bounds
if [ "$ARRAY_INDEX" -ge ${#SUMSTAT_FILES[@]} ]; then
    echo "Task $SGE_TASK_ID: No file assigned (only ${#SUMSTAT_FILES[@]} files available)"
    exit 0
fi

SUMSTAT_FILE="${SUMSTAT_FILES[$ARRAY_INDEX]}"
SUMSTAT_NAME=$(basename "${SUMSTAT_FILE}" .imputed.ma)
OUTPUT_NAME="${OUTPUT_DIR}/${SUMSTAT_NAME}_SBayesRC_calc"

echo "Task $SGE_TASK_ID processing: $SUMSTAT_NAME"

# Check if output file already exists
if [ -f "${OUTPUT_NAME}.snpRes" ]; then
    echo "✓ Output file ${OUTPUT_NAME}.snpRes already exists. Skipping..."
    exit 0
fi

echo "Calculating weights for ${SUMSTAT_NAME}..."

# Run GCTB
if "$GCTB_BINARY" \
    --gwas-summary "$SUMSTAT_FILE" \
    --ldm-eigen "$LD_DIR"/ukbEUR_HM3 \
    --annot "$ANNO_DIR"/annot_baseline2.2.txt \
    --sbayes RC \
    --thread 1 \
    --out "$OUTPUT_NAME"; then
    echo "✓ Successfully processed ${SUMSTAT_NAME}"
    exit 0
else
    echo "✗ ERROR: Failed to process ${SUMSTAT_NAME}" >&2
    exit 1
fi


# qsub /path/to/home/TRS_project/sumstats/scripts/SBayesRC_calc_hapmap_sumstats_array.sh
