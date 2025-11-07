#!/bin/bash
#SBATCH --account=def-yergeaue-ab_cpu
#SBATCH --time=00-02:00:00       # The duration in D-HH:MM:SS format of each task
#SBATCH --mem=4G                 # Total memory for the job (adjust if needed)
#SBATCH --mail-type=END,FAIL     # Notifications for job done & fail
#SBATCH --mail-user=julius-eyiuche.nweze@inrs.ca  # Set your email address
#===============================================================================
# Title          : 5_run_singlem_summarise_OTU.sh
# Description    : A simple loop to generate SingleM script for each sample read
# Author         : Julius Eyiuche Nweze
# Date           : 03/03/2025
# Version        : V.0
# Usage          : sbatch 5_run_singlem_summarise_OTU.sh
#===============================================================================

# Load necessary modules
module load StdEnv/2023 apptainer/1.2.4

# Directory containing the profile files
PROFILE_DIR="/project/6004719/projects/GROW/MT_2024-10-06/bacteria/Julius/SingleM/OTU"

# Combine OTU tables using singlem summarise
singularity run /home/juli24/scratch/singlem_0.18.3.sif summarise \
    --input-otu-tables $PROFILE_DIR/*.otu_table.tsv \
    --output-otu-table $PROFILE_DIR/Combined.otu_table.csv

# Check if the command succeeded
if [ $? -eq 0 ]; then
    echo "Successfully combined OTU tables."
else
    echo "Error combining OTU tables." >&2
fi
