#!/bin/bash

# Create output directory if it doesn't exist
mkdir -p New_scripts

# Loop through all .sh files in Script directory
for file in Scripts/*.sh; do
    filename=$(basename "$file")
    newfile="New_scripts/$filename"

    sed -e '/# Designate the output directory for OTU/ {
s|.*|# Input directory\
INPUT="/project/6004719/projects/GROW/MT_2024-10-06/bacteria/Julius/PCyCDB/Result2";|
}' \
        -e '/^OUT_otu=/d' \
        -e '/# Designate the output directory for Profile/ {
s|.*|# Output directory\
OUTPUT="/project/6004719/projects/GROW/MT_2024-10-06/bacteria/Julius/PCyCDB/TPM";|
}' \
        -e '/^OUT_profile=/d' \
        -e '/# Create the output directory if it doesn'\''t exist/,/mkdir -p "\$OUT_profile"/c\
# Create output directory if it doesn'\''t exist\
mkdir -p "$OUTPUT"' \
        -e '/# Run singlem with the selected files/,/otu_table.tsv"/c\
coverm contig -1 $R1s -2 $R2s -r $INPUT/PCyCDB_phospho_genes_filtered.fasta -o $OUTPUT/PCyCDB_genes.coverm.tsv -m tpm --min-read-percent-identity 0.95 -p bwa-mem --threads 5' \
        "$file" > "$newfile"

    echo "Updated: $file → $newfile"
done
