#!/bin/bash
#SBATCH --account=def-yergeaue-ab       # Replace with your own account
#SBATCH --mem-per-cpu=10G               # Memory allocation per CPU
#SBATCH --time=00-00:10:00               # Job time limit (DD-HH:MM)
#SBATCH --cpus-per-task=10               # Number of CPU cores (adjust as needed)
#SBATCH --job-name=NA-HMM      # Job name for easier identification
#SBATCH --output=NA-HMM-%j.out # Standard output file
#SBATCH --error=NA-HMM-%j.err  # Standard error file
#===============================================================================
# Title          : 11a_NA_deg.sh
# Description    : Search for genes
# Author         : Julius Eyiuche Nweze
# Date           : 03/03/2025
# Version        : V.0
# Usage          : sbatch 11a_NA_deg.sh
#===============================================================================

# Load necessary modules
module load StdEnv/2023 hmmer/3.4 

# Define input and output paths
HMM_FILE="/project/6004719/projects/GROW/MT_2024-10-06/bacteria/Julius/NA_DEG/NA_deg.hmm"
FASTA_FILE="/project/6004719/projects/GROW/MT_2024-10-06/bacteria/Julius/hmm_index/Contigs_renamed.faa"
OUTPUT_DIR="/project/6004719/projects/GROW/MT_2024-10-06/bacteria/Julius/NA_DEG"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run hmmsearch with Stockholm alignment output
hmmsearch -E 1e-5 -A "$OUTPUT_DIR/NA_genes.sto" "$HMM_FILE" "$FASTA_FILE"

# Run hmmsearch with tabular output
hmmsearch -E 1e-5 --tblout "$OUTPUT_DIR/NA_genes.txt" "$HMM_FILE" "$FASTA_FILE"

# Convert Stockholm alignment to FASTA
esl-reformat fasta "$OUTPUT_DIR/NA_genes.sto" > "$OUTPUT_DIR/NA_genes.fa"
