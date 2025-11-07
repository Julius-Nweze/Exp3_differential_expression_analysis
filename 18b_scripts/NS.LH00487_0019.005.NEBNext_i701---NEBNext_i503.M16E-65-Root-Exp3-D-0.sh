#!/bin/bash
#SBATCH --account=rrg-yergeaue
#SBATCH --time=5:00:00       # The duration in HH:MM:SS format of each task in the array
#SBATCH --cpus-per-task=10 
#SBATCH --job-name=SingleM
#SBATCH --mem=10G
#SBATCH --export=ALL
#SBATCH --mail-type=END,FAIL      # Notifications for job done & fail
#SBATCH --mail-user=julius-eyiuche.nweze@inrs.ca  # Set your email address

# Load necessary modules
module load apptainer-suid/1.1

# Define and export the sample name
sample_name=NS.LH00487_0019.005.NEBNext_i701---NEBNext_i503.M16E-65-Root-Exp3-D-0
echo "Processing sample: $sample_name"

# Input directory
INPUT="/project/6004719/projects/GROW/MT_2024-10-06/bacteria/Julius/PCyCDB/Result2";

# Output directory
OUTPUT="/project/6004719/projects/GROW/MT_2024-10-06/bacteria/Julius/PCyCDB/TPM";

# Define the R1 and R2 fastq files for the current sample
R1s=("/project/6004719/projects/GROW/MT_2024-10-06/bacteria/raw_reads/NS.LH00487_0019.005.NEBNext_i701---NEBNext_i503.M16E-65-Root-Exp3-D-0_R1.fastq.gz")
R2s=("/project/6004719/projects/GROW/MT_2024-10-06/bacteria/raw_reads/NS.LH00487_0019.005.NEBNext_i701---NEBNext_i503.M16E-65-Root-Exp3-D-0_R2.fastq.gz")

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT"


# Check if R1 and R2 files exist before running SingleM
if [[ -f "${R1s[0]}" && -f "${R2s[0]}" ]]; then
coverm contig -1 $R1s -2 $R2s -r $INPUT/PCyCDB_phospho_genes_filtered.fasta -o $OUTPUT/PCyCDB_genes.coverm.tsv -m tpm --min-read-percent-identity 0.95 -p bwa-mem --threads 5
else
    echo "Warning: R1 or R2 files not found for $sample_name"
fi
