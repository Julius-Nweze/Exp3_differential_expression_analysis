#!/bin/bash
#SBATCH --account=def-yergeaue-ab_cpu
#SBATCH --time=00-01:00:00              # Duration (D-HH:MM:SS)
#SBATCH --mem=4G                        # Total memory for the job
#SBATCH --mail-type=END,FAIL            # Notifications for job completion or failure
#SBATCH --mail-user=julius-eyiuche.nweze@inrs.ca  # Your email for notifications
#===============================================================================
# Title          : 9_run_deseq2.sh
# Description    : Runs DESeq2 differential expression analysis using R
# Author         : Julius Eyiuche Nweze
# Date           : 03/03/2025
# Version        : v1.0
# Usage          : sbatch 9_run_deseq2.sh
#===============================================================================

# Exit immediately if a command exits with a non-zero status
set -euo pipefail

# Define base and output directories
BASE_DIR="/home/juli24/projects/def-yergeaue-ab/juli24/Exp3"
OUTPUT_DIR="${BASE_DIR}/DEG_analysis"

# Create output directory if it doesn’t exist
mkdir -p "$OUTPUT_DIR"

# Input file paths
INPUT_FILE="${BASE_DIR}/merged_gene_abundance.tsv"
R_SCRIPT="${BASE_DIR}/MT_bacteria.R"
OUTPUT_FILE="${OUTPUT_DIR}/MT_DEG_bacteria.csv"

# Log info
echo "==========================================="
echo "Starting DESeq2 differential expression job"
echo "Base directory  : $BASE_DIR"
echo "Input file      : $INPUT_FILE"
echo "R script        : $R_SCRIPT"
echo "Output file     : $OUTPUT_FILE"
echo "Job started at  : $(date)"
echo "==========================================="

# Run DESeq2 R script
Rscript "$R_SCRIPT" "$INPUT_FILE" "$OUTPUT_FILE"

echo "==========================================="
echo "Job completed successfully at: $(date)"
echo "Results saved in: $OUTPUT_FILE"
echo "==========================================="

