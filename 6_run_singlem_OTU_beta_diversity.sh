#!/bin/bash
#SBATCH --account=def-yergeaue-ab_cpu
#SBATCH --time=00-00:20:00       # The duration in D-HH:MM:SS format of each task
#SBATCH --mem=4G                 # Total memory for the job (adjust if needed)
#SBATCH --mail-type=END,FAIL      # Notifications for job done & fail
#SBATCH --mail-user=julius-eyiuche.nweze@inrs.ca  # Set your email address
#===============================================================================
# Title          : 6_run_singlem_OTU_beta_diversity.sh
# Description    : A simple loop to generate SingleM script for each sample read
# Author         : Julius Eyiuche Nweze
# Date           : 03/03/2025
# Version        : V.0
# Usage          : sbatch 6_run_singlem_OTU_beta_diversity.sh
#===============================================================================

# Load necessary modules
module load StdEnv/2023 apptainer/1.2.4

# Directory containing the profile files
PROFILE_DIR="/project/6004719/projects/GROW/MT_2024-10-06/bacteria/Julius/SingleM/OTU"

# Specify the beta diversity directory for conversion scripts
CONVERT_TO_EBD="/home/juli24/projects/rrg-yergeaue/juli24/R/SingleM/ExpressBetaDiversity/scripts"

# Specify the ExpressBetaDiversity directory
EXPRESS_BETA_DIVERSITY="/home/juli24/projects/rrg-yergeaue/juli24/R/SingleM/ExpressBetaDiversity/bin"

# Output directory for Beta Diversity results
BETA_DIR="$PROFILE_DIR/Beta_diversity"
mkdir -p "$BETA_DIR"

# Define output filenames
UNIFRAC_FILE="$BETA_DIR/otu_table-.S3.5.ribosomal_protein_S2_rpsB.unifrac.S3.12.ribosomal_L1.unifrac.S3.48.rplD_bact.unifrac.S3.58.tsf.unifrac.S3.13.ribosomal_S9.unifrac.S3.1.ribosomal_protein_L2_rplB.unifrac.S3.59.uS11_bact.unifrac"
EBD_FILE="$BETA_DIR/otu_table.ebd"
BD_OUTPUT="$BETA_DIR/Combined.otu_table_edited_BD"

# Bray-Curtis beta diversity calculation using SingleM summarise
singularity run /home/juli24/scratch/singlem_0.18.3.sif summarise \
    --input-otu-table "$FILTERED_FILE" \
    --unifrac-by-otu "$UNIFRAC_FILE"

# Check if SingleM summarise command succeeded
if [ $? -eq 0 ]; then
    echo "Successfully calculated Bray-Curtis beta diversity using SingleM."

    # Convert the Unifrac output to EBD format
    python "$CONVERT_TO_EBD/convertToEBD.py" "$UNIFRAC_FILE" "$EBD_FILE"

    # Check if the conversion was successful
    if [ $? -eq 0 ]; then
        echo "Successfully converted Unifrac output to EBD format."

        # Run ExpressBetaDiversity
        "$EXPRESS_BETA_DIVERSITY/ExpressBetaDiversity" -s "$EBD_FILE" -c Bray-Curtis -p "$BD_OUTPUT"

        # Check if ExpressBetaDiversity command succeeded
        if [ $? -eq 0 ]; then
            echo "Successfully calculated Express Beta Diversity."
        else
            echo "Error calculating Express Beta Diversity." >&2
        fi
    else
        echo "Error converting Unifrac output to EBD format." >&2
    fi
else
    echo "Error calculating Bray-Curtis beta diversity using SingleM." >&2
fi
