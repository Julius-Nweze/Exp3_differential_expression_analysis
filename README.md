# README.md

.. image:: https://zenodo.org/badge/1091796028.svg
  :target: https://doi.org/10.5281/zenodo.17553524


## Project Information
- **Title:** Microbial and Plant Gene Expression Response to Naphthenic Acids Fractional Components  
- **Author:** Julius Eyiuche Nweze  
- **Date:** 2025-11-06  
- **Email:** julipeale2001@gmail.com  

## Overview
This repository contains a collection of R, Bash, and HMM scripts for analyzing active microbial communities, as well as for identifying, mapping, and characterizing genes involved in naphthenic acid fractional component (NAFC) degradation, microbial functional pathways, and plant transcriptomic responses.
The project integrates metatranscriptomic workflows focusing on active microbial community profiling, hydrocarbon-degrading genes, differential gene expression, and biogeochemical cycling pathways that contribute to NAFC degradation.

## Table of Contents
1. [Setup](#setup)  
   - [Dependencies](#dependencies)  
   - [Installation](#installation)  
2. [Data](#data)  
   - [Input Files](#input-files)  
   - [Output Files](#output-files)  
3. [Analysis Workflow](#analysis-workflow)  
   - [Functional Gene Detection (HMM)](#functional-gene-detection-hmm)  
   - [Gene Mapping and Expression Analysis](#gene-mapping-and-expression-analysis)  
   - [Microbial and Plant DEG Analysis](#microbial-and-plant-deg-analysis)  
   - [Microbial Community Profiling](#microbial-community-profiling)  
4. [Scripts Overview](#scripts-overview)  
5. [Figures and Outputs](#figures-and-outputs)  
6. [Reproducibility](#reproducibility)  

## Setup

### Dependencies
- R (≥ 4.3.0)
- Python (optional, for parsing or data cleaning)
- bash / SLURM environment
- Required R packages: `tidyverse`, `DESeq2`, `ggplot2`, `vegan`, `rstatix`, `pheatmap`, `kableExtra`, etc.

### Installation
Clone this repository and ensure the execution permissions for shell scripts:
```bash
git clone https://github.com/julipeale2001/Naphthenic_Acid_DEG.git
cd Naphthenic_Acid_DEG
chmod +x *.sh
```

## Data

### Input Files
- `merged_gene_abundance.tsv` — gene abundance data across samples  
- `NA_deg.hmm`, `CANT-HYD.hmm`, `MCycDB.hmm` — HMM profiles for hydrocarbon degradation pathways  
- `Metadata.csv` — sample metadata  
- GTF files for plant gene annotation  

### Output Files
- DEG results (bacteria and plant) in CSV format  
- Pathway summaries and visualizations in RMarkdown reports  
- Normalized TPM tables for mapped genes  

## Analysis Workflow

### Functional Gene Detection (HMM)
Scripts such as `11a_NA_deg.sh`, `12a_CANT_hmm.sh`, `13a_HADEG.sh`,`15a_NCyc.sh`, `16a_Phospho_gene.sh`, `17a_Sulf.sh` and `14a_MCycDB.sh` perform HMM-based gene searches for hydrocarbon-degrading and nutrient cycling genes across metagenomic or transcriptomic datasets.

### Gene Mapping and Expression Analysis
Scripts (`18a_Gene_mapping_TPM.sh`, `18b_generate_script_Gene_mapping_TPM.sh`, `8_MT_bacteria.R` and `19b_MT_plant.R`) generate TPM-normalized gene expression matrices, followed by DESeq2-based differential expression analysis (`9_run_deseq2.sh`).

### Microbial and Plant DEG Analysis
RMarkdown files like `10_DEGs_Naphthenic_acid.Rmd` and `19d_DEGs_Plant_Naphthenic_acid.Rmd` analyze and visualize differentially expressed genes (DEGs) in microbial and plant systems under NA exposure.

### Microbial Community Profiling
Scripts such as `4_run_singlem_summarise_Profile.sh`, `5_run_singlem_summarise_OTU.sh`, `6_run_singlem_OTU_beta_diversity.sh` and `7_Microbial_community.Rmd` process and analyze microbial community composition using SingleM and diversity metrics.

## Scripts Overview
| Script/File | Description |
|--------------|-------------|
| `1_generate_script.sh` | Generates submission scripts in 2_Scripts folder for multiple samples analysis of microbial community in the raw reads |
| `2_Scripts` | A folder containing the generated read script for multiple samples analysis of microbial community in the raw reads |
| `3_submit_all.sh` | Submits all generated SLURM jobs to microbial community in the raw reads |
| `4_run_singlem_summarise_Profile.sh` | Summarizes SingleM Profile tables |
| `5_run_singlem_summarise_OTU.sh` | Summarizes SingleM OTU tables |
| `6_run_singlem_OTU_beta_diversity.sh` | Uses SingleM OTU tables to calculate beta diversity between samples |
| `7_Microbial_community.Rmd` | Visualizes active microbial community and beta diversity results |
| `8_MT_bacteria.R` | Processes metatranscriptomic bacterial gene data (differentially expressed genes) using the gene abundance |
| `9_run_deseq2.sh` | Runs DESeq2 differential expression pipeline using 9_run_deseq2.sh |
| `10_DEGs_Naphthenic_acid.Rmd` | Visualizes the number, taxonomic origins and genes that are differentially expressed |
| `11a_NA_deg.sh` | Searches NA-degrading genes using HMMs |
| `11b_NA_deg.hmm` | HMM profile of the previously identified NA-degrading genes |
| `11c_NA_HMM.Rmd` | Visualizes the abundance of previously identified NA-degrading genes |
| `12a_CANT_hmm.sh` | Searches for genes using CANT-HYD: a curated database of phylogeny-derived hidden markov models for annotation of marker genes involved in hydrocarbon degradation |
| `12b_CANT-HYD.hmm` | Hydrocarbon degradation marker database (CANT-HYD) from https://github.com/dgittins/CANT-HYD-HydrocarbonBiodegradation |
| `12c_CANT.Rmd` | Visualizes the abundance of hydrocarbon degradation marker genes from CANT-HYD database |
| `13a_HADEG.sh` | Searches for genes using HADEG: A Curated Hydrocarbon Aerobic Degradation Enzymes and Genes Database from https://github.com/jarojasva/HADEG |
| `13b_HADEG.Rmd` | Visualizes the abundance of hydrocarbon degradation marker genes from HADEG database |
| `14a_MCycDB.sh` | Searches for genes using MCycDB: a curated database for comprehensively profiling methane cycling processes of environmental microbiomes from https://github.com/qichao1984/MCycDB |
| `14b_MCycDB.Rmd` | Visualizes the abundance of methane cycling marker genes from MCycDB database |
| `15a_NCyc.sh` | Searches for genes using NCyc：a curated integrative database for fast and accurate metagenomic profiling of nitrogen cycling genes from https://github.com/qichao1984/NCyc |
| `15b_NCyc.Rmd` | Visualizes the abundance of nitrogen cycling marker genes from NCyc database |
| `16a_Phospho_gene.sh` | Searches for genes using Phosphorus-cycling-database (PCyCDB) from https://github.com/ZengJiaxiong/Phosphorus-cycling-database |
| `16b_Extract_gene_position.sh` | Extracts genomic information of candidate genes |
| `16c_Extract_gene_position2.sh` | Extracts genomic coordinates of candidate genes |
| `16d_PCyCDB.Rmd` | Visualizes the abundance of nitrogen cycling marker genes from PCyCDB database |
| `17a_Sulf.sh` | Visualizes the abundance of sulfur cycling marker genes from SulfDB database (https://github.com/qichao1984/SCycDB) |
| `17b_Sulf.Rmd` | Visualizes the abundance of sulfur cycling marker genes from SulfDB database |
| `18_generate_script_Gene_mapping_TPM.sh` | Generate scripts in 18c_scripts for gene-to-TPM mapping to reads |
| `18b_scripts` | A folder that contains the generated scripts for gene-to-TPM mapping to reads |
| `19a_Plant_extract_gene_info_gtf.sh` | Extracts plant genetic coordinate information from a GTF file |
| `19b_MT_plant.R` | Processes metatranscriptomic plant gene data (differentially expressed genes) using the plant gene abundance |
| `19c_DEGs_Plant_Naphthenic_acid.Rmd` | Performs plant DEG analysis and visualization |

## Figures and Outputs
All figures are generated automatically from RMarkdown scripts (`.Rmd`) and saved in the `Figures/` directory.   
The generated figures are modified using Inkscape, a free and open-source software vector graphics editor. 
Each `.Rmd` script produces interactive plots, bar charts, and heatmaps summarizing gene abundance, taxonomy, and expression.

## Reproducibility
- **R version:** 4.3.3  
- **Operating System:** Unix/Linux (HPC environment)  
- **Scheduler:** SLURM  
- **Citation:** Please cite this repository and the relevant HMM databases used in your analysis.

---
© 2025 Julius Eyiuche Nweze. All rights reserved.
