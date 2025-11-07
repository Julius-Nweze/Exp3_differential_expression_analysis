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
sample_name=NS.LH00487_0019.005.NEBNext_i705---NEBNext_i501.M8G-379-Root-Exp3-D19
echo "Processing sample: $sample_name"

# Designate the output directory for OTU
OUT_otu=/project/6004719/projects/GROW/MT_2024-10-06/bacteria/Julius/SingleM/OTU

# Designate the output directory for Profile
OUT_profile=/project/6004719/projects/GROW/MT_2024-10-06/bacteria/Julius/SingleM/Profile

# Define the R1 and R2 fastq files for the current sample
R1s=("/project/6004719/projects/GROW/MT_2024-10-06/bacteria/raw_reads/NS.LH00487_0019.005.NEBNext_i705---NEBNext_i501.M8G-379-Root-Exp3-D19_R1.fastq.gz")
R2s=("/project/6004719/projects/GROW/MT_2024-10-06/bacteria/raw_reads/NS.LH00487_0019.005.NEBNext_i705---NEBNext_i501.M8G-379-Root-Exp3-D19_R2.fastq.gz")

# Create the output directory if it doesn't exist
mkdir -p "$OUT_otu"
mkdir -p "$OUT_profile"


# Check if R1 and R2 files exist before running SingleM
if [[ -f "${R1s[0]}" && -f "${R2s[0]}" ]]; then
    # Run singlem with the selected files
    singularity run /home/juli24/scratch/singlem_0.18.3.sif pipe -1 "${R1s[0]}" -2 "${R2s[0]}" \
        -p "$OUT_profile/${sample_name}.profile.tsv" --otu-table "$OUT_otu/NS.LH00487_0019.005.NEBNext_i705---NEBNext_i501.M8G-379-Root-Exp3-D19.otu_table.tsv"
else
    echo "Warning: R1 or R2 files not found for $sample_name"
fi
