#!/bin/bash
#===============================================================================
# Title          : 1_generate_script.sh
# Description    : A simple loop to generate SingleM script for each sample read
# Author         : Julius Eyiuche Nweze
# Date           : 03/03/2025
# Version        : V.0
# Usage          : sbatch 1_generate_script.sh
#===============================================================================
#SBATCH --account=rrg-yergeaue
#SBATCH --time=00:10:00          # Duration (HH:MM:SS)
#SBATCH --cpus-per-task=2
#SBATCH --job-name=SingleM_generate_script
#SBATCH --mem=5G
#SBATCH --export=ALL
#SBATCH --mail-type=END,FAIL    # Notifications for job completion/failure
#SBATCH --mail-user=julius-eyiuche.nweze@inrs.ca  # Your email address
#SBATCH --output=generate-scripts_%j.out
#SBATCH --error=generate-scripts_%j.err
#================================================================================

# Sample.txt should contain one directory name per line, e.g., NS.1871.001, NS.1871.002, etc.
sample_file="Sample.txt"

# Create the output directory if it doesn't exist
mkdir -p Scripts

# Loop over each line in Sample.txt
while read -r sample_name; do
    # Create a new script file with the sample_name in Scripts/SingleM
    output_script="Scripts/${sample_name}.sh"
    
    # Copy the template script contents
    cat <<EOF > "$output_script"
#!/bin/bash
#SBATCH --account=rrg-yergeaue
#SBATCH --time=5:00:00       # The duration in HH:MM:SS format of each task in the array
#SBATCH --cpus-per-task=10 
#SBATCH --job-name=SingleM
#SBATCH --mem=10G
#SBATCH --export=ALL
#SBATCH --mail-type=END,FAIL      # Notifications for job done & fail
#SBATCH --mail-user=julius-eyiuche.nweze@inrs.ca  # Set your email address
#SBATCH --output="${sample_name}_%j.out"
#SBATCH --error="${sample_name}_%j.err"

# Load necessary modules
module load apptainer-suid/1.1

# Define and export the sample name
sample_name=${sample_name}
echo "Processing sample: \$sample_name"

# Designate the output directory for OTU
OUT_otu=/project/6004719/projects/GROW/MT_2024-10-06/bacteria/Julius/SingleM/OTU

# Designate the output directory for Profile
OUT_profile=/project/6004719/projects/GROW/MT_2024-10-06/bacteria/Julius/SingleM/Profile

# Define the R1 and R2 fastq files for the current sample
R1s=("/project/6004719/projects/GROW/MT_2024-10-06/bacteria/raw_reads/${sample_name}_R1.fastq.gz")
R2s=("/project/6004719/projects/GROW/MT_2024-10-06/bacteria/raw_reads/${sample_name}_R2.fastq.gz")

# Create the output directory if it doesn't exist
mkdir -p "\$OUT_otu"
mkdir -p "\$OUT_profile"


# Check if R1 and R2 files exist before running SingleM
if [[ -f "\${R1s[0]}" && -f "\${R2s[0]}" ]]; then
    # Run singlem with the selected files
    singularity run /home/juli24/scratch/singlem_0.18.3.sif pipe -1 "\${R1s[0]}" -2 "\${R2s[0]}" \\
        -p "\$OUT_profile/\${sample_name}.profile.tsv" --otu-table "\$OUT_otu/${sample_name}.otu_table.tsv"
else
    echo "Warning: R1 or R2 files not found for \$sample_name"
fi
EOF

    # Make the generated script executable
    chmod +x "$output_script"
    echo "Generated script: $output_script"
    
done < "$sample_file"

