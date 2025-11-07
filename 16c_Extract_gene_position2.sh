#!/bin/bash
#SBATCH --account=def-yergeaue-ab         # Replace with your own account
#SBATCH --mem-per-cpu=50G                 # Memory allocation per CPU
#SBATCH --time=2-00:00:00                 # Job time limit (DD-HH:MM)
#SBATCH --cpus-per-task=10                # Number of CPU cores
#SBATCH --job-name=anvio-contig-db        # Job name
#SBATCH --output=anvio-contig-db-%j.out   # Standard output file
#SBATCH --error=anvio-contig-db-%j.err    # Standard error file

#===========================================================================================================
# Script to filter and extract phosphorus-cycling genes from BLAST output (PCyCDB) based on identity (≥30%)
# and alignment length (≥25 aa). Outputs include .bed files and FASTA of filtered hits with updated headers.
#===========================================================================================================

# Load necessary modules
module load StdEnv/2023 gcc/12.3 seqtk/1.4

# Set the input paths
INPUT="/project/6004719/projects/GROW/MT_2024-10-06/bacteria"
BLAST_RESULT="$INPUT/Julius/PCyCDB/Result/PCyCDB.output.txt"
CONTIGS="$INPUT/assembly/Contigs.fasta"

# Set the output directory and files
OUTDIR="$INPUT/Julius/PCyCDB/Result2"
OUTPUT_FILTERED="$OUTDIR/PCyCDB.filtered.output.txt"
BED_RAW="$OUTDIR/PCyCDB.phosphate_genes_filtered_ID.bed"
BED_FINAL="$OUTDIR/PCyCDB.phosphate_genes_filtered.bed"
GENES_FASTA="$OUTDIR/PCyCDB_phospho_genes_filtered.fasta"

# Create output directory if not already there
mkdir -p "$OUTDIR"

# Step 1: Filter BLAST results by identity ≥30% and alignment length ≥25 aa
awk '$3 >= 30 && $5 >= 25' "$BLAST_RESULT" > "$OUTPUT_FILTERED"

# Step 2: Check if any hits passed the filter
if [[ ! -s "$OUTPUT_FILTERED" ]]; then
    echo "Warning: No genes passed the filtering threshold. Exiting script."
    exit 1
fi

# Step 3: Extract original gene IDs with positions
cut -f2,10,11 "$OUTPUT_FILTERED" > "$BED_RAW"

# Step 4: Modify gene IDs to match contig names (strip suffixes)
awk -F'\t' '{
    orig_id = $1
    split(orig_id, a, "_")
    mod_id = a[1]
    for (i = 2; i < length(a); i++) {
        mod_id = mod_id "_" a[i]
    }
    print orig_id, mod_id, $2, $3
}' OFS='\t' "$BED_RAW" > "$BED_FINAL"

# Step 5: Extract sequences using seqtk from modified IDs
seqtk subseq "$CONTIGS" <(cut -f2 "$BED_FINAL") > "$GENES_FASTA"

# Step 6: Rename FASTA headers to include original ID, modified ID, and positions
TMP_FASTA="$GENES_FASTA.tmp"

awk 'FNR==NR {
    map[$2] = $1 "|" $2 "|" $3 "|" $4
    next
}
(/^>/) {
    seq_id = substr($0, 2)
    print ">" map[seq_id]
    next
}
{ print }' "$BED_FINAL" "$GENES_FASTA" > "$TMP_FASTA" && mv "$TMP_FASTA" "$GENES_FASTA"

echo "✅ Gene filtering and extraction complete."
echo "🔍 Output FASTA: $GENES_FASTA with headers: >original_ID|modified_ID|start|end"
