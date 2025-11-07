#!/bin/bash
#SBATCH --account=def-yergeaue-ab_cpu
#SBATCH --time=00-02:00:00       # The duration in D-HH:MM:SS format of each task in the array
#SBATCH --mem=4G                 # Total memory for the job (adjust if needed)
#SBATCH --mail-type=END,FAIL     # Notifications for job done & fail
#SBATCH --mail-user=julius-eyiuche.nweze@inrs.ca  # Set your email address
#===============================================================================
# Title          : 4_run_singlem_summarise_Profile.sh
# Description    : A simple loop to generate SingleM script for each sample read
# Author         : Julius Eyiuche Nweze
# Date           : 03/03/2025
# Version        : V.0
# Usage          : sbatch 4_run_singlem_summarise_Profile.sh
#===============================================================================

# Load necessary modules
module load StdEnv/2023 apptainer/1.2.4

# Directory containing the profile files
PROFILE_DIR="/project/6004719/projects/GROW/MT_2024-10-06/bacteria/Julius/SingleM/Profile"

# Loop through all profile files in the directory
for profile_file in "$PROFILE_DIR"/*.profile.tsv; do
    # Get the base name of the profile file (without the directory and extension)
    basename=$(basename "$profile_file" .profile.tsv)

    # Construct the output file paths
    species_by_site_output="$PROFILE_DIR/$basename.species_by_site_phylum.csv"
    profile_with_extras_output="$PROFILE_DIR/$basename.with_extras.tsv"
    species_by_site_level="phylum"
   
    # singlem summarise: Converting a taxonomic profile to a more convenient format
    # Run the SingleM summarise command within the Singularity container
    singularity run /home/juli24/scratch/singlem_0.18.3.sif summarise \
        --input-taxonomic-profile "$profile_file" \
        --output-species-by-site-relative-abundance "$species_by_site_output" \
        --output-species-by-site-level "$species_by_site_level"

    # Long form with extras
    singularity run /home/juli24/scratch/singlem_0.18.3.sif summarise \
        --input-taxonomic-profile "$profile_file" \
        --output-taxonomic-profile-with-extras "$profile_with_extras_output"

    # Check if the commands succeeded
    if [ $? -eq 0 ]; then
        echo "Successfully processed $profile_file"
    else
        echo "Error processing $profile_file" >&2
    fi

done
