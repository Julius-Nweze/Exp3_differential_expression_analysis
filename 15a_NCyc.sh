#!/bin/bash
#SBATCH --account=def-yergeaue-ab         # Replace with your own account
#SBATCH --mem-per-cpu=20G                 # Memory allocation per CPU
#SBATCH --time=01-00:00:00                 # Job time limit (DD-HH:MM)
#SBATCH --cpus-per-task=10                # Number of CPU cores (adjust as needed)
#SBATCH --job-name=NCyc-db        # Job name for easier identification
#SBATCH --output=NCyc-db-%j.out   # Standard output file
#SBATCH --error=NCyc-db-%j.err    # Standard error file
#===============================================================================
# Title          : 14a_NCyc.sh
# Description    : A simple loop to generate SingleM script for each sample read
# Author         : Julius Eyiuche Nweze
# Date           : 03/03/2025
# Version        : V.0
# Usage          : sbatch 14a_NCyc.sh
#===============================================================================

# Load required modules
module load StdEnv/2023 gcc/12.3 blast+/2.14.1

# Path to the predicted proteins from contigs (FASTA)
ASS="/home/juli24/projects/def-yergeaue-ab/juli24/Exp3/Julius/blast_index/Contigs_renamed.faa"

# Path to the phosphorus cycling protein database (FASTA)
PHOS="/home/juli24/projects/def-yergeaue-ab/juli24/Exp3/Julius/NCyc/NCyc_95.faa"

# Output files and directories
OUTDIR="/home/juli24/projects/def-yergeaue-ab/juli24/Exp3/Julius/NCyc/Result"
# OUTPUT_INDEX="/home/juli24/projects/def-yergeaue-ab/juli24/Exp3/Julius/blast_index/Contigs_renamed.faa"
OUTPUT="$OUTDIR/NCyc.output.txt"
OUTPUT2="$OUTDIR/NCyc.output2.txt"
OUTPUT_GENE="$OUTDIR/NCyc_gene.faa"

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Run BLASTP pipeline
# makeblastdb -in "$ASS" -out "$OUTPUT_INDEX" -dbtype prot -parse_seqids

# Full blastp output with annotations
blastp -query "$PHOS" -db "$ASS" -out "$OUTPUT" \
  -outfmt '6 qseqid sseqid pident qcovus length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames stitle qcovs' \
  -evalue 1e-5 &&

# Extract only subject sequence IDs
blastp -query "$PHOS" -db "$ASS" -out "$OUTPUT2" \
  -outfmt '6 sseqid' -evalue 1e-5 &&

# Retrieve matching protein sequences
blastdbcmd -db "$ASS" -entry_batch "$OUTPUT2" -out "$OUTPUT_GENE" -outfmt '%f'
