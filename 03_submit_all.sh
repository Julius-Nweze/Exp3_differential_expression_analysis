#!/bin/bash
#===============================================================================
# Title          : 3_submit_all.sh
# Description    : Submit script to cluster
# Author         : Julius Eyiuche Nweze
# Date           : 03/03/2025
# Version        : V.0
# Usage          : sbatch 3_submit_all.sh
#===============================================================================

# Change to the directory containing the scripts
cd Scripts || exit

# Loop through each .sh script and submit it
for script in *.sh; do
    sbatch "$script"
done
