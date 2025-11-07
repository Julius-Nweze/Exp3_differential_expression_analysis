#!/bin/bash
#SBATCH --account=def-yergeaue-ab           # Replace with your own account
#SBATCH --mem-per-cpu=10G                   # Memory allocation per CPU
#SBATCH --time=00:20:00                     # Job time limit (HH:MM:SS)
#SBATCH --cpus-per-task=10                  # Number of CPU cores
#SBATCH --job-name=CANT-HMM                 # Job name for easier identification
#SBATCH --output=CANT-HMM-%j.out      # Standard output log
#SBATCH --error=CANT-HMM-%j.err       # Standard error log
#===============================================================================
# Title          : 12a_CANT_hmm.sh
# Description    : Search for genes
# Author         : Julius Eyiuche Nweze
# Date           : 03/03/2025
# Version        : V.0
# Usage          : sbatch 12a_CANT_hmm.sh
#===============================================================================

# Load required modules
module load StdEnv/2023 hmmer/3.4

# Define input and output paths
HMM_FILE="/project/6004719/projects/GROW/MT_2024-10-06/bacteria/Julius/CANT/CANT-HYD.hmm"
FASTA_FILE="/project/6004719/projects/GROW/MT_2024-10-06/bacteria/Julius/hmm_index/Contigs_renamed.faa"
OUTPUT_DIR="/project/6004719/projects/GROW/MT_2024-10-06/bacteria/Julius/CANT"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run hmmsearch with Stockholm alignment output
hmmsearch -E 1e-5 -A "$OUTPUT_DIR/CANT.sto" "$HMM_FILE" "$FASTA_FILE"

# Run hmmsearch with tabular output
hmmsearch -E 1e-5 --tblout "$OUTPUT_DIR/CANT.txt" "$HMM_FILE" "$FASTA_FILE"

# Convert Stockholm alignment to FASTA
esl-reformat fasta "$OUTPUT_DIR/CANT.sto" > "$OUTPUT_DIR/CANT_genes.fa"

