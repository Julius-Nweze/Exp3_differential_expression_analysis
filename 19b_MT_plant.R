#!/usr/bin/env Rscript
# Install required libraries ---------------------------------------------------
################################################################################
.libPaths(c('~/R/library', .libPaths())) # Uncomment if you have no write access to R path
# https://hbctraining.github.io/DGE_workshop/lessons/05_DGE_DESeq2_analysis2.html
repo <- "http://cran.wu.ac.at"
#lib.loc <- Sys.getenv("R_LIBS_USER")

# update.packages(
# #    lib.loc, 
#     repos = repo,
#     ask = FALSE
#)

.cran_libs <- c(
"devtools", # The tools required to efficiently manage other packages
"data.table", # Optimized for fast reading of large files (fread)
"readr", # For reading .csv files efficiently
"tibble", # For tibble manipulations (data frames with extra features)
"tidyverse", # Helps to transform and better present data. 
"ggplot2", # For creating declarative and customizable graphics
"magrittr", # Forward and backward “pipe”-like operator
"dplyr", # A grammar of data manipulation
"readxl", # designed for reading Excel files
"DESeq2" # For RNA-seq count data analysis and differential expression
) 

.inst <- .cran_libs %in% installed.packages()
if (any(!.inst)) {
   install.packages(.cran_libs[!.inst],
                    repos = repo
                    )
}

.bioc_libs <- c(
  #"multtest", #Resampling-based multiple hypothesis testing
)

.bioc_inst <- .bioc_libs %in% installed.packages()
if (any(!.bioc_inst)) {
   if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
   BiocManager::install(ask = F)  # upgrade bioC packages
   BiocManager::install(.bioc_libs[!.bioc_inst], ask = F)
}

.local_libs <- c()

.inst <- names(.local_libs) %in% installed.packages()
if (any(!.inst)) {
   install.packages(paste0("~/R/", .local_libs[!.inst]) ,repos = NULL, type = "source")
}

.github_libs <- c(
  "wilkelab/ggtext", # Improved text rendering support for 'ggplot2' 
  "ACCLAB/dabestr" # Data Analysis using Bootstrap-Coupled Estimation 
)

.github_lib_names <- stringr::str_replace(.github_libs, ".*/(.*)$", "\\1")

.github_inst <- .github_lib_names %in% installed.packages()
if (any(!.github_inst)) {
  devtools::install_github(.github_libs[!.github_inst],
                           
                           dependencies = TRUE)
}

# Load packages into session, and print package version
(loaded.libs <- sapply(c(.cran_libs, .bioc_libs, names(.local_libs), .github_lib_names), require, character.only = TRUE))
if (!all(loaded.libs)) {stop(paste("Package(s):", names(loaded.libs[loaded.libs == FALSE]), "could not be loaded"))}
sapply(c(.cran_libs, .bioc_libs, names(.local_libs), .github_lib_names), packageVersion)
################################################################################


# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

################################################################################# For RO sample


# Load the large gene_MT abundance in TSV file
# data <- fread("merged_gene_abundance.tsv", sep = "\t")
data <- fread(input_file, sep = "\t")

########################################################## Run for each water type: RO
# Load the metadata
metadata <- read_tsv("metadata.tsv")
metadata_full <- read_xlsx("metadata_Root.xlsx")

# Merge metadata
metadata <-  metadata %>% dplyr::left_join(metadata_full, by=c("SampleID", "Time"), multiple = "all", relationship = "many-to-many")

# Set the first column as row names in data
data <- data %>%
  column_to_rownames(var = colnames(data)[1])  # Set the first column as row names

# Set the first column as row names in metadata
metadata <- metadata %>% filter(Water_type == "RO") %>%
  column_to_rownames(var = colnames(metadata)[1])  # Set the first column as row names


# Find mismatches
mismatched_columns <- setdiff(colnames(data), rownames(metadata))
mismatched_rows <- setdiff(rownames(metadata), colnames(data))

# Print mismatches
if (length(mismatched_columns) > 0) {
  message("Columns in data not found in metadata row names: ", paste(mismatched_columns, collapse = ", "))
}

if (length(mismatched_rows) > 0) {
  message("Rows in metadata not found in data column names: ", paste(mismatched_rows, collapse = ", "))
}

# Subset data to keep only columns present in metadata row names
data <- data[, colnames(data) %in% rownames(metadata)]

# Reorder columns in data to match row names in metadata
data <- data[, match(rownames(metadata), colnames(data))]

# Set column names in data to match row names in metadata
colnames(data) <- rownames(metadata)

# Verify alignment
identical(colnames(data), rownames(metadata)) # Should return TRUE
################################################################################



# Construct DESEQDataSet Object---------------------------- 
################################################################################
# Make sure the metadata has the right structure and factor levels
# Ensure that Time and Sample_type are factors
metadata$Time <- factor(metadata$Time, levels = c("D-0", "D6", "D19", "D40", "D60"))



# Create the DESeq2 Dataset: In the design formula ~ Sample_type * Time, the DESeq2 model is structured to analyze the interaction effects of Sample_type and Time.
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadata,
                              design = ~ Time)

################################################################################


# Pre-filtering----------------------------
#################################################################################
# Here we perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples. 
# The count of 10 is a reasonable choice for bulk RNA-seq. A recommendation for the minimal number of samples is to specify the smallest group size, e.g. here there are 3 treated samples.
# If there are not discrete groups, one can use the minimal number of samples where non-zero counts would be considered interesting. 
# One can also omit this step entirely and just rely on the independent filtering procedures available in results(), either IHW or genefilter https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilt. 
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds2 <- dds[keep,]
################################################################################

# Run DESeq2 Analysis----------------------------
################################################################################
# Now you can run the DESeq2 analysis: Results tables are generated using the function results, which extracts a results table with log2 fold changes, p values and adjusted p values. 
# With no additional arguments to results, the log2 fold change and Wald test p value will be for the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable over the reference level (see previous note on factor levels). 
# However, the order of the variables of the design do not matter so long as the user specifies the comparison to build a results table for, using the name or contrast arguments of results.
dds_res <- DESeq(dds2)
# resultsNames(dds_res)
################################################################################

# Extract Results for Each Comparison----------------------------
################################################################################
# Extract interaction effects: Extract results with a specific alpha threshold (e.g., 0.05; typically set to 0.05 by default)


results_interaction_d6_RO <- results(dds_res, name = "Time_D6_vs_D.0")
results_interaction_d19_RO <- results(dds_res, name = "Time_D19_vs_D.0")
results_interaction_d40_RO <- results(dds_res, name = "Time_D40_vs_D.0")
results_interaction_d60_RO <- results(dds_res, name = "Time_D60_vs_D.0")


# Define your threshold values for significance
padj_threshold <- 0.05
log2fc_threshold <- 1

# Convert to a dataframe with filtering the DEGs for D6
Unfiltered_degs_d6_RO <- as.data.frame(results_interaction_d6_RO)
Unfiltered_degs_d6_RO$Gene <- rownames(Unfiltered_degs_d6_RO) # Add the Gene column from the row names
Unfiltered_degs_d6_RO$Time <- "D6 vs D-0" # Add a Time column to each filtered DEGs dataframe
Unfiltered_degs_d6_RO$Water_type <- "RO" # Add a Water_type column to each filtered DEGs dataframe

# Convert to a dataframe with filtering the DEGs for D19
Unfiltered_degs_d19_RO <- as.data.frame(results_interaction_d19_RO)
Unfiltered_degs_d19_RO$Gene <- rownames(Unfiltered_degs_d19_RO) # Add the Gene column from the row names
Unfiltered_degs_d19_RO$Time <- "D19 vs D-0" # Add a Time column to each filtered DEGs dataframe
Unfiltered_degs_d19_RO$Water_type <- "RO" # Add a Water_type column to each filtered DEGs dataframe

# Convert to a dataframe with filtering the DEGs for D40
Unfiltered_degs_d40_RO <- as.data.frame(results_interaction_d40_RO)
Unfiltered_degs_d40_RO$Gene <- rownames(Unfiltered_degs_d40_RO) # Add the Gene column from the row names
Unfiltered_degs_d40_RO$Time <- "D40 vs D-0" # Add a Time column to each filtered DEGs dataframe
Unfiltered_degs_d40_RO$Water_type <- "RO" # Add a Water_type column to each filtered DEGs dataframe

# Convert to a dataframe with filtering the DEGs for D60
Unfiltered_degs_d60_RO <- as.data.frame(results_interaction_d60_RO)
Unfiltered_degs_d60_RO$Gene <- rownames(Unfiltered_degs_d60_RO) # Add the Gene column from the row names
Unfiltered_degs_d60_RO$Time <- "D60 vs D-0" # Add a Time column to each filtered DEGs dataframe
Unfiltered_degs_d60_RO$Water_type <- "RO" # Add a Water_type column to each filtered DEGs dataframe

# Merge the data frames:
Unfiltered_merged_degs_RO <- rbind(Unfiltered_degs_d6_RO, Unfiltered_degs_d19_RO, Unfiltered_degs_d40_RO, Unfiltered_degs_d60_RO)
Unfiltered_merged_degs_RO$SeqType <- "Metatranscriptome" # Add a SeqType column to each filtered DEGs dataframe

# Filter for significant DEGs and rearrange columns
Filtered_merged_degs_RO <- subset(Unfiltered_merged_degs_RO, padj < padj_threshold & abs(log2FoldChange) > log2fc_threshold)
Filtered_merged_degs_RO <- Filtered_merged_degs_RO[, c("SeqType", "Gene", "Time", setdiff(colnames(Filtered_merged_degs_RO), c("SeqType", "Gene", "Time")))]

# Save the merged DEGs dataframe to a CSV
# write.csv(Filtered_merged_degs_RO, "Filtered_merged_degs_RO.csv")
write.csv(Filtered_merged_degs_RO, output_file, row.names = FALSE)



# For RO at D6 comparison
# Load the large gene_MT abundance in TSV file
# data <- fread("merged_gene_abundance.tsv", sep = "\t")
data <- fread(input_file, sep = "\t")

########################################################## Run for each water type: RO using D6
# Load the metadata
metadata <- read_tsv("metadata.tsv")
metadata_full <- read_xlsx("metadata_Root.xlsx")

# Merge metadata
metadata <-  metadata %>% dplyr::left_join(metadata_full, by=c("SampleID", "Time"), multiple = "all", relationship = "many-to-many")

# Set the first column as row names in data
data <- data %>%
  column_to_rownames(var = colnames(data)[1])  # Set the first column as row names

# Set the first column as row names in metadata
metadata <- metadata %>% filter(Water_type == "RO") %>%
  filter(Time %in% c("D6", "D19", "D40", "D60")) %>%
  
  column_to_rownames(var = colnames(metadata)[1])  # Set the first column as row names


# Find mismatches
mismatched_columns <- setdiff(colnames(data), rownames(metadata))
mismatched_rows <- setdiff(rownames(metadata), colnames(data))

# Print mismatches
if (length(mismatched_columns) > 0) {
  message("Columns in data not found in metadata row names: ", paste(mismatched_columns, collapse = ", "))
}

if (length(mismatched_rows) > 0) {
  message("Rows in metadata not found in data column names: ", paste(mismatched_rows, collapse = ", "))
}

# Subset data to keep only columns present in metadata row names
data <- data[, colnames(data) %in% rownames(metadata)]

# Reorder columns in data to match row names in metadata
data <- data[, match(rownames(metadata), colnames(data))]

# Set column names in data to match row names in metadata
colnames(data) <- rownames(metadata)

# Verify alignment
identical(colnames(data), rownames(metadata)) # Should return TRUE
################################################################################



# Construct DESEQDataSet Object---------------------------- 
################################################################################
# Make sure the metadata has the right structure and factor levels
# Ensure that Time and Sample_type are factors
metadata$Time <- factor(metadata$Time, levels = c("D6", "D19", "D40", "D60"))



# Create the DESeq2 Dataset: In the design formula ~ Sample_type * Time, the DESeq2 model is structured to analyze the interaction effects of Sample_type and Time.
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadata,
                              design = ~ Time)

################################################################################


# Pre-filtering----------------------------
#################################################################################
# Here we perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples. 
# The count of 10 is a reasonable choice for bulk RNA-seq. A recommendation for the minimal number of samples is to specify the smallest group size, e.g. here there are 3 treated samples.
# If there are not discrete groups, one can use the minimal number of samples where non-zero counts would be considered interesting. 
# One can also omit this step entirely and just rely on the independent filtering procedures available in results(), either IHW or genefilter https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilt. 
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds2 <- dds[keep,]
################################################################################

# Run DESeq2 Analysis----------------------------
################################################################################
# Now you can run the DESeq2 analysis: Results tables are generated using the function results, which extracts a results table with log2 fold changes, p values and adjusted p values. 
# With no additional arguments to results, the log2 fold change and Wald test p value will be for the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable over the reference level (see previous note on factor levels). 
# However, the order of the variables of the design do not matter so long as the user specifies the comparison to build a results table for, using the name or contrast arguments of results.
dds_res <- DESeq(dds2)
# resultsNames(dds_res)
################################################################################

# Extract Results for Each Comparison----------------------------
################################################################################
# Extract interaction effects: Extract results with a specific alpha threshold (e.g., 0.05; typically set to 0.05 by default)



results_interaction_d19_RO <- results(dds_res, name = "Time_D19_vs_D6")
results_interaction_d40_RO <- results(dds_res, name = "Time_D40_vs_D6")
results_interaction_d60_RO <- results(dds_res, name = "Time_D60_vs_D6")


# Define your threshold values for significance
padj_threshold <- 0.05
log2fc_threshold <- 1


# Convert to a dataframe with filtering the DEGs for D19
Unfiltered_degs_d19_RO <- as.data.frame(results_interaction_d19_RO)
Unfiltered_degs_d19_RO$Gene <- rownames(Unfiltered_degs_d19_RO) # Add the Gene column from the row names
Unfiltered_degs_d19_RO$Time <- "D19 vs D6" # Add a Time column to each filtered DEGs dataframe
Unfiltered_degs_d19_RO$Water_type <- "RO" # Add a Water_type column to each filtered DEGs dataframe

# Convert to a dataframe with filtering the DEGs for D40
Unfiltered_degs_d40_RO <- as.data.frame(results_interaction_d40_RO)
Unfiltered_degs_d40_RO$Gene <- rownames(Unfiltered_degs_d40_RO) # Add the Gene column from the row names
Unfiltered_degs_d40_RO$Time <- "D40 vs D6" # Add a Time column to each filtered DEGs dataframe
Unfiltered_degs_d40_RO$Water_type <- "RO" # Add a Water_type column to each filtered DEGs dataframe

# Convert to a dataframe with filtering the DEGs for D60
Unfiltered_degs_d60_RO <- as.data.frame(results_interaction_d60_RO)
Unfiltered_degs_d60_RO$Gene <- rownames(Unfiltered_degs_d60_RO) # Add the Gene column from the row names
Unfiltered_degs_d60_RO$Time <- "D60 vs D6" # Add a Time column to each filtered DEGs dataframe
Unfiltered_degs_d60_RO$Water_type <- "RO" # Add a Water_type column to each filtered DEGs dataframe

# Merge the data frames:
Unfiltered_merged_degs_RO <- rbind(Unfiltered_degs_d19_RO, Unfiltered_degs_d40_RO, Unfiltered_degs_d60_RO)
Unfiltered_merged_degs_RO$SeqType <- "Metatranscriptome" # Add a SeqType column to each filtered DEGs dataframe

# Filter for significant DEGs and rearrange columns
Filtered_merged_degs_RO <- subset(Unfiltered_merged_degs_RO, padj < padj_threshold & abs(log2FoldChange) > log2fc_threshold)
Filtered_merged_degs_RO_D6 <- Filtered_merged_degs_RO[, c("SeqType", "Gene", "Time", setdiff(colnames(Filtered_merged_degs_RO), c("SeqType", "Gene", "Time")))]

# Save the merged DEGs dataframe to a CSV
# write.csv(Filtered_merged_degs_RO_D6, "Filtered_merged_degs_RO_D6.csv")
write.csv(Filtered_merged_degs_RO_D6, output_file, row.names = FALSE)






# For OSPW sample


# Load the data: Compartment (Sediments), Sample_type (No_plant and Carex) and Water_type (OSPW)----------------------------
################################################################################
# Load the large gene_MT abundance in TSV file
# data <- fread("merged_gene_abundance.tsv", sep = "\t")
data <- fread(input_file, sep = "\t")
data <- data[, lapply(.SD, function(x) if (is.numeric(x)) x + 1 else x)]
########################################################## Run for each water type: OSPW
# Load the metadata
metadata <- read_tsv("metadata.tsv")
metadata_full <- read_xlsx("metadata_Root.xlsx")

# Merge metadata
metadata <-  metadata %>% dplyr::left_join(metadata_full, by=c("SampleID", "Time"), multiple = "all", relationship = "many-to-many")

# Set the first column as row names in data
data <- data %>%
  column_to_rownames(var = colnames(data)[1])  # Set the first column as row names

# Set the first column as row names in metadata
metadata <- metadata %>% filter(Water_type == "OSPW") %>%
  column_to_rownames(var = colnames(metadata)[1])  # Set the first column as row names


# Find mismatches
mismatched_columns <- setdiff(colnames(data), rownames(metadata))
mismatched_rows <- setdiff(rownames(metadata), colnames(data))

# Print mismatches
if (length(mismatched_columns) > 0) {
  message("Columns in data not found in metadata row names: ", paste(mismatched_columns, collapse = ", "))
}

if (length(mismatched_rows) > 0) {
  message("Rows in metadata not found in data column names: ", paste(mismatched_rows, collapse = ", "))
}

# Subset data to keep only columns present in metadata row names
data <- data[, colnames(data) %in% rownames(metadata)]

# Reorder columns in data to match row names in metadata
data <- data[, match(rownames(metadata), colnames(data))]

# Set column names in data to match row names in metadata
colnames(data) <- rownames(metadata)

# Verify alignment
identical(colnames(data), rownames(metadata)) # Should return TRUE
################################################################################



# Construct DESEQDataSet Object---------------------------- 
################################################################################
# Make sure the metadata has the right structure and factor levels
# Ensure that Time and Sample_type are factors
metadata$Time <- factor(metadata$Time, levels = c("D-0", "D6", "D19", "D40", "D60"))



# Create the DESeq2 Dataset: In the design formula ~ Sample_type * Time, the DESeq2 model is structured to analyze the interaction effects of Sample_type and Time.
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadata,
                              design = ~ Time)

################################################################################


# Pre-filtering----------------------------
#################################################################################
# Here we perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples. 
# The count of 10 is a reasonable choice for bulk RNA-seq. A recommendation for the minimal number of samples is to specify the smallest group size, e.g. here there are 3 treated samples.
# If there are not discrete groups, one can use the minimal number of samples where non-zero counts would be considered interesting. 
# One can also omit this step entirely and just rely on the independent filtering procedures available in results(), either IHW or genefilter https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilt. 
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds2 <- dds[keep,]
################################################################################

# Run DESeq2 Analysis----------------------------
################################################################################
# Now you can run the DESeq2 analysis: Results tables are generated using the function results, which extracts a results table with log2 fold changes, p values and adjusted p values. 
# With no additional arguments to results, the log2 fold change and Wald test p value will be for the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable over the reference level (see previous note on factor levels). 
# However, the order of the variables of the design do not matter so long as the user specifies the comparison to build a results table for, using the name or contrast arguments of results.
# I got the error below: estimating size factors. Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,: every gene contains at least one zero, cannot compute log geometric means

dds_res <- DESeq(dds2)
# resultsNames(dds_res)
################################################################################

# Extract Results for Each Comparison----------------------------
################################################################################
# Extract interaction effects: Extract results with a specific alpha threshold (e.g., 0.05; typically set to 0.05 by default)


results_interaction_d6_OSPW <- results(dds_res, name = "Time_D6_vs_D.0")
results_interaction_d19_OSPW <- results(dds_res, name = "Time_D19_vs_D.0")
results_interaction_d40_OSPW <- results(dds_res, name = "Time_D40_vs_D.0")
results_interaction_d60_OSPW <- results(dds_res, name = "Time_D60_vs_D.0")


# Define your threshold values for significance
padj_threshold <- 0.05
log2fc_threshold <- 1

# Convert to a dataframe with filtering the DEGs for D6
Unfiltered_degs_d6_OSPW <- as.data.frame(results_interaction_d6_OSPW)
Unfiltered_degs_d6_OSPW$Gene <- rownames(Unfiltered_degs_d6_OSPW) # Add the Gene column from the row names
Unfiltered_degs_d6_OSPW$Time <- "D6 vs D-0" # Add a Time column to each filtered DEGs dataframe
Unfiltered_degs_d6_OSPW$Water_type <- "OSPW" # Add a Water_type column to each filtered DEGs dataframe

# Convert to a dataframe with filtering the DEGs for D19
Unfiltered_degs_d19_OSPW <- as.data.frame(results_interaction_d19_OSPW)
Unfiltered_degs_d19_OSPW$Gene <- rownames(Unfiltered_degs_d19_OSPW) # Add the Gene column from the row names
Unfiltered_degs_d19_OSPW$Time <- "D19 vs D-0" # Add a Time column to each filtered DEGs dataframe
Unfiltered_degs_d19_OSPW$Water_type <- "OSPW" # Add a Water_type column to each filtered DEGs dataframe

# Convert to a dataframe with filtering the DEGs for D40
Unfiltered_degs_d40_OSPW <- as.data.frame(results_interaction_d40_OSPW)
Unfiltered_degs_d40_OSPW$Gene <- rownames(Unfiltered_degs_d40_OSPW) # Add the Gene column from the row names
Unfiltered_degs_d40_OSPW$Time <- "D40 vs D-0" # Add a Time column to each filtered DEGs dataframe
Unfiltered_degs_d40_OSPW$Water_type <- "OSPW" # Add a Water_type column to each filtered DEGs dataframe

# Convert to a dataframe with filtering the DEGs for D60
Unfiltered_degs_d60_OSPW <- as.data.frame(results_interaction_d60_OSPW)
Unfiltered_degs_d60_OSPW$Gene <- rownames(Unfiltered_degs_d60_OSPW) # Add the Gene column from the row names
Unfiltered_degs_d60_OSPW$Time <- "D60 vs D-0" # Add a Time column to each filtered DEGs dataframe
Unfiltered_degs_d60_OSPW$Water_type <- "OSPW" # Add a Water_type column to each filtered DEGs dataframe

# Merge the data frames:
Unfiltered_merged_degs_OSPW <- rbind(Unfiltered_degs_d6_OSPW, Unfiltered_degs_d19_OSPW, Unfiltered_degs_d40_OSPW, Unfiltered_degs_d60_OSPW)
Unfiltered_merged_degs_OSPW$SeqType <- "Metatranscriptome" # Add a SeqType column to each filtered DEGs dataframe

# Filter for significant DEGs and rearrange columns
Filtered_merged_degs_OSPW <- subset(Unfiltered_merged_degs_OSPW, padj < padj_threshold & abs(log2FoldChange) > log2fc_threshold)
Filtered_merged_degs_OSPW <- Filtered_merged_degs_OSPW[, c("SeqType", "Gene", "Time", setdiff(colnames(Filtered_merged_degs_OSPW), c("SeqType", "Gene", "Time")))]

# Save the merged DEGs dataframe to a CSV
# write.csv(Filtered_merged_degs_OSPW, "Filtered_merged_degs_OSPW.csv")
write.csv(Filtered_merged_degs_OSPW, output_file, row.names = FALSE)




# For OSPW at D6 comparison
################################################################################
# Load the large gene_MT abundance in TSV file
# data <- fread("merged_gene_abundance.tsv", sep = "\t")
data <- fread(input_file, sep = "\t")
data <- data[, lapply(.SD, function(x) if (is.numeric(x)) x + 1 else x)]
########################################################## Run for each water type: OSPW
# Load the metadata
metadata <- read_tsv("metadata.tsv")
metadata_full <- read_xlsx("metadata_Root.xlsx")

# Merge metadata
metadata <-  metadata %>% dplyr::left_join(metadata_full, by=c("SampleID", "Time"), multiple = "all", relationship = "many-to-many")

# Set the first column as row names in data
data <- data %>%
  column_to_rownames(var = colnames(data)[1])  # Set the first column as row names

# Set the first column as row names in metadata
metadata <- metadata %>% filter(Water_type == "OSPW") %>%
  filter(Time %in% c("D6", "D19", "D40", "D60")) %>%
  column_to_rownames(var = colnames(metadata)[1])  # Set the first column as row names


# Find mismatches
mismatched_columns <- setdiff(colnames(data), rownames(metadata))
mismatched_rows <- setdiff(rownames(metadata), colnames(data))

# Print mismatches
if (length(mismatched_columns) > 0) {
  message("Columns in data not found in metadata row names: ", paste(mismatched_columns, collapse = ", "))
}

if (length(mismatched_rows) > 0) {
  message("Rows in metadata not found in data column names: ", paste(mismatched_rows, collapse = ", "))
}

# Subset data to keep only columns present in metadata row names
data <- data[, colnames(data) %in% rownames(metadata)]

# Reorder columns in data to match row names in metadata
data <- data[, match(rownames(metadata), colnames(data))]

# Set column names in data to match row names in metadata
colnames(data) <- rownames(metadata)

# Verify alignment
identical(colnames(data), rownames(metadata)) # Should return TRUE
################################################################################



# Construct DESEQDataSet Object---------------------------- 
################################################################################
# Make sure the metadata has the right structure and factor levels
# Ensure that Time and Sample_type are factors
metadata$Time <- factor(metadata$Time, levels = c("D6", "D19", "D40", "D60"))



# Create the DESeq2 Dataset: In the design formula ~ Sample_type * Time, the DESeq2 model is structured to analyze the interaction effects of Sample_type and Time.
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadata,
                              design = ~ Time)

################################################################################


# Pre-filtering----------------------------
#################################################################################
# Here we perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples. 
# The count of 10 is a reasonable choice for bulk RNA-seq. A recommendation for the minimal number of samples is to specify the smallest group size, e.g. here there are 3 treated samples.
# If there are not discrete groups, one can use the minimal number of samples where non-zero counts would be considered interesting. 
# One can also omit this step entirely and just rely on the independent filtering procedures available in results(), either IHW or genefilter https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilt. 
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds2 <- dds[keep,]
################################################################################

# Run DESeq2 Analysis----------------------------
################################################################################
# Now you can run the DESeq2 analysis: Results tables are generated using the function results, which extracts a results table with log2 fold changes, p values and adjusted p values. 
# With no additional arguments to results, the log2 fold change and Wald test p value will be for the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable over the reference level (see previous note on factor levels). 
# However, the order of the variables of the design do not matter so long as the user specifies the comparison to build a results table for, using the name or contrast arguments of results.
dds_res <- DESeq(dds2)
# resultsNames(dds_res)
################################################################################

# Extract Results for Each Comparison----------------------------
################################################################################
# Extract interaction effects: Extract results with a specific alpha threshold (e.g., 0.05; typically set to 0.05 by default)

results_interaction_d19_OSPW <- results(dds_res, name = "Time_D19_vs_D6")
results_interaction_d40_OSPW <- results(dds_res, name = "Time_D40_vs_D6")
results_interaction_d60_OSPW <- results(dds_res, name = "Time_D60_vs_D6")


# Define your threshold values for significance
padj_threshold <- 0.05
log2fc_threshold <- 1


# Convert to a dataframe with filtering the DEGs for D19
Unfiltered_degs_d19_OSPW <- as.data.frame(results_interaction_d19_OSPW)
Unfiltered_degs_d19_OSPW$Gene <- rownames(Unfiltered_degs_d19_OSPW) # Add the Gene column from the row names
Unfiltered_degs_d19_OSPW$Time <- "D19 vs D6" # Add a Time column to each filtered DEGs dataframe
Unfiltered_degs_d19_OSPW$Water_type <- "OSPW" # Add a Water_type column to each filtered DEGs dataframe

# Convert to a dataframe with filtering the DEGs for D40
Unfiltered_degs_d40_OSPW <- as.data.frame(results_interaction_d40_OSPW)
Unfiltered_degs_d40_OSPW$Gene <- rownames(Unfiltered_degs_d40_OSPW) # Add the Gene column from the row names
Unfiltered_degs_d40_OSPW$Time <- "D40 vs D6" # Add a Time column to each filtered DEGs dataframe
Unfiltered_degs_d40_OSPW$Water_type <- "OSPW" # Add a Water_type column to each filtered DEGs dataframe

# Convert to a dataframe with filtering the DEGs for D60
Unfiltered_degs_d60_OSPW <- as.data.frame(results_interaction_d60_OSPW)
Unfiltered_degs_d60_OSPW$Gene <- rownames(Unfiltered_degs_d60_OSPW) # Add the Gene column from the row names
Unfiltered_degs_d60_OSPW$Time <- "D60 vs D6" # Add a Time column to each filtered DEGs dataframe
Unfiltered_degs_d60_OSPW$Water_type <- "OSPW" # Add a Water_type column to each filtered DEGs dataframe

# Merge the data frames:
Unfiltered_merged_degs_OSPW <- rbind(Unfiltered_degs_d19_OSPW, Unfiltered_degs_d40_OSPW, Unfiltered_degs_d60_OSPW)
Unfiltered_merged_degs_OSPW$SeqType <- "Metatranscriptome" # Add a SeqType column to each filtered DEGs dataframe

# Filter for significant DEGs and rearrange columns
Filtered_merged_degs_OSPW <- subset(Unfiltered_merged_degs_OSPW, padj < padj_threshold & abs(log2FoldChange) > log2fc_threshold)
Filtered_merged_degs_OSPW_D6 <- Filtered_merged_degs_OSPW[, c("SeqType", "Gene", "Time", setdiff(colnames(Filtered_merged_degs_OSPW), c("SeqType", "Gene", "Time")))]

# Save the merged DEGs dataframe to a CSV
# write.csv(Filtered_merged_degs_OSPW_D6, "Filtered_merged_degs_OSPW_D6.csv")
write.csv(Filtered_merged_degs_OSPW_D6, output_file, row.names = FALSE)













################################################ Water type and Time: OSPW vs RO using D-0 as the initial time point


################################################################################
# Load the large gene_MT abundance in TSV file
# data <- fread("merged_gene_abundance.tsv", sep = "\t")
data <- fread(input_file, sep = "\t")
data <- data[, lapply(.SD, function(x) if (is.numeric(x)) x + 1 else x)]

# Load the metadata
metadata <- read_tsv("metadata.tsv")
metadata_full <- read_xlsx("metadata_Root.xlsx")

# Merge metadata
metadata <-  metadata %>% dplyr::left_join(metadata_full, by=c("SampleID", "Time"), multiple = "all", relationship = "many-to-many")

# Set the first column as row names in data
data <- data %>%
  column_to_rownames(var = colnames(data)[1])  # Set the first column as row names

# Set the first column as row names in metadata
metadata <- metadata %>%
  column_to_rownames(var = colnames(metadata)[1])  # Set the first column as row names


# Find mismatches
mismatched_columns <- setdiff(colnames(data), rownames(metadata))
mismatched_rows <- setdiff(rownames(metadata), colnames(data))

# Print mismatches
if (length(mismatched_columns) > 0) {
  message("Columns in data not found in metadata row names: ", paste(mismatched_columns, collapse = ", "))
}

if (length(mismatched_rows) > 0) {
  message("Rows in metadata not found in data column names: ", paste(mismatched_rows, collapse = ", "))
}

# Subset data to keep only columns present in metadata row names
data <- data[, colnames(data) %in% rownames(metadata)]

# Reorder columns in data to match row names in metadata
data <- data[, match(rownames(metadata), colnames(data))]

# Set column names in data to match row names in metadata
colnames(data) <- rownames(metadata)

# Verify alignment
identical(colnames(data), rownames(metadata)) # Should return TRUE
################################################################################



# Construct DESEQDataSet Object---------------------------- 
################################################################################
# Make sure the metadata has the right structure and factor levels
# Ensure that Time and Sample_type are factors
metadata$Time <- factor(metadata$Time, levels = c("D-0", "D6", "D19", "D40", "D60"))
metadata$Water_type <- factor(metadata$Water_type, levels = c("RO", "OSPW"))
# metadata$Greenhouse <- factor(metadata$Greenhouse, levels = c("E", "G"))


# Create the DESeq2 Dataset: In the design formula ~ Sample_type * Time, the DESeq2 model is structured to analyze the interaction effects of Sample_type and Time.
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadata,
                              design = ~ Water_type * Time)
# The interaction terms indicate how the effect of time differs between the "Root" and "OSPW" sample types. 
# These terms are critical for understanding whether the presence of Root alters the time-dependent expression changes relative to the "OSPW" condition.
################################################################################


# Pre-filtering----------------------------
#################################################################################
# Here we perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples. 
# The count of 10 is a reasonable choice for bulk RNA-seq. A recommendation for the minimal number of samples is to specify the smallest group size, e.g. here there are 3 treated samples.
# If there are not discrete groups, one can use the minimal number of samples where non-zero counts would be considered interesting. 
# One can also omit this step entirely and just rely on the independent filtering procedures available in results(), either IHW or genefilter https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilt. 
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds2 <- dds[keep,]
################################################################################

# Run DESeq2 Analysis----------------------------
################################################################################
# Now you can run the DESeq2 analysis: Results tables are generated using the function results, which extracts a results table with log2 fold changes, p values and adjusted p values. 
# With no additional arguments to results, the log2 fold change and Wald test p value will be for the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable over the reference level (see previous note on factor levels). 
# However, the order of the variables of the design do not matter so long as the user specifies the comparison to build a results table for, using the name or contrast arguments of results.
dds_res <- DESeq(dds2)
# resultsNames(dds_res)
################################################################################

# Extract Results for Each Comparison----------------------------
################################################################################
# Extract interaction effects: Extract results with a specific alpha threshold (e.g., 0.05; typically set to 0.05 by default)
# Water_typeRO.TimeD6 indciates how the effect of root differs from OSPW specifically in terms of the change in gene expression between D6 and D-0
# The term "Water_typeRO.TimeD6" represents the interaction between Sample_type ("Root" vs "OSPW") and Time (D6 vs D-0).
# It tells you how the difference in gene expression from D-0 to D6 is affected by the presence of Root, relative to OSPW.
# If "Water_typeRO.TimeD6" has a positive log2 fold change, it suggests that Root increases the gene expression difference between D6 and D-0 compared to what happens in OSPW. 
# Conversely, a negative value would suggest that Root decreases this difference in gene expression compared to OSPW.

results_interaction_d6 <- results(dds_res, name = "Water_typeOSPW.TimeD6")
results_interaction_d19 <- results(dds_res, name = "Water_typeOSPW.TimeD19")
results_interaction_d40 <- results(dds_res, name = "Water_typeOSPW.TimeD40")
results_interaction_d60 <- results(dds_res, name = "Water_typeOSPW.TimeD60")
# results_interaction_d8 <- results(dds_res, alpha = alpha, contrast = c("Treatment", "NO", "HA"))




# Define your threshold values for significance
padj_threshold <- 0.05
log2fc_threshold <- 1

# Convert to a dataframe with filtering the DEGs for D6
Unfiltered_degs_d6 <- as.data.frame(results_interaction_d6)
Unfiltered_degs_d6$Gene <- rownames(Unfiltered_degs_d6) # Add the Gene column from the row names
Unfiltered_degs_d6$Time <- "D6 vs D-0" # Add a Time column to each filtered DEGs dataframe

# Convert to a dataframe with filtering the DEGs for D19
Unfiltered_degs_d19 <- as.data.frame(results_interaction_d19)
Unfiltered_degs_d19$Gene <- rownames(Unfiltered_degs_d19) # Add the Gene column from the row names
Unfiltered_degs_d19$Time <- "D19 vs D-0" # Add a Time column to each filtered DEGs dataframe

# Convert to a dataframe with filtering the DEGs for D40
Unfiltered_degs_d40 <- as.data.frame(results_interaction_d40)
Unfiltered_degs_d40$Gene <- rownames(Unfiltered_degs_d40) # Add the Gene column from the row names
Unfiltered_degs_d40$Time <- "D40 vs D-0" # Add a Time column to each filtered DEGs dataframe

# Convert to a dataframe with filtering the DEGs for D60
Unfiltered_degs_d60 <- as.data.frame(results_interaction_d60)
Unfiltered_degs_d60$Gene <- rownames(Unfiltered_degs_d60) # Add the Gene column from the row names
Unfiltered_degs_d60$Time <- "D60 vs D-0" # Add a Time column to each filtered DEGs dataframe

# Merge the data frames:
Unfiltered_merged_degs <- rbind(Unfiltered_degs_d6, Unfiltered_degs_d19, Unfiltered_degs_d40, Unfiltered_degs_d60)
Unfiltered_merged_degs$SeqType <- "Metatranscriptome" # Add a Time column to each filtered DEGs dataframe

# Filter for significant DEGs and rearrange columns
Filtered_merged_degs <- subset(Unfiltered_merged_degs, padj < padj_threshold & abs(log2FoldChange) > log2fc_threshold)
Filtered_merged_degs <- Filtered_merged_degs[, c("SeqType", "Gene", "Time", setdiff(colnames(Filtered_merged_degs), c("SeqType", "Gene", "Time")))]

# Save the merged DEGs dataframe to a CSV
# write.csv(Filtered_merged_degs, "Filtered_merged_degs_OSPWvsTime.csv")
write.csv(Filtered_merged_degs, output_file, row.names = FALSE)









################################################ Water type and Time: OSPW vs RO at each time point



################################################ Water type and Time (Filter time 0)


# Load the data: Compartment (Sediments), Sample_type (No_plant and Carex) and Water_type (OSPW)----------------------------
################################################################################
# Load the large gene_MT abundance in TSV file
# data <- fread("merged_gene_abundance.tsv", sep = "\t")
data <- fread(input_file, sep = "\t")
data <- data[, lapply(.SD, function(x) if (is.numeric(x)) x + 1 else x)]

# Load the metadata
metadata <- read_tsv("metadata.tsv")
metadata_full <- read_xlsx("metadata_Root.xlsx")

# Merge metadata
metadata <-  metadata %>% dplyr::left_join(metadata_full, by=c("SampleID", "Time"), multiple = "all", relationship = "many-to-many")

# Set the first column as row names in data
data <- data %>%
  column_to_rownames(var = colnames(data)[1])  # Set the first column as row names

# Set the first column as row names in metadata
metadata <- metadata %>% filter(Time == "D-0") %>%
  column_to_rownames(var = colnames(metadata)[1])  # Set the first column as row names


# Find mismatches
mismatched_columns <- setdiff(colnames(data), rownames(metadata))
mismatched_rows <- setdiff(rownames(metadata), colnames(data))

# Print mismatches
if (length(mismatched_columns) > 0) {
  message("Columns in data not found in metadata row names: ", paste(mismatched_columns, collapse = ", "))
}

if (length(mismatched_rows) > 0) {
  message("Rows in metadata not found in data column names: ", paste(mismatched_rows, collapse = ", "))
}

# Subset data to keep only columns present in metadata row names
data <- data[, colnames(data) %in% rownames(metadata)]

# Reorder columns in data to match row names in metadata
data <- data[, match(rownames(metadata), colnames(data))]

# Set column names in data to match row names in metadata
colnames(data) <- rownames(metadata)

# Verify alignment
identical(colnames(data), rownames(metadata)) # Should return TRUE
################################################################################



# Construct DESEQDataSet Object---------------------------- 
################################################################################
# Make sure the metadata has the right structure and factor levels
# Ensure that Time and Sample_type are factors
# metadata$Time <- factor(metadata$Time, levels = c("D-0", "D6", "D19", "D40", "D60"))
metadata$Water_type <- factor(metadata$Water_type, levels = c("RO", "OSPW"))
# metadata$Greenhouse <- factor(metadata$Greenhouse, levels = c("E", "G"))


# Create the DESeq2 Dataset: In the design formula ~ Sample_type * Time, the DESeq2 model is structured to analyze the interaction effects of Sample_type and Time.
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadata,
                              design = ~ Water_type)
# The interaction terms indicate how the effect of time differs between the "Root" and "OSPW" sample types. 
# These terms are critical for understanding whether the presence of Root alters the time-dependent expression changes relative to the "OSPW" condition.
################################################################################


# Pre-filtering----------------------------
#################################################################################
# Here we perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples. 
# The count of 10 is a reasonable choice for bulk RNA-seq. A recommendation for the minimal number of samples is to specify the smallest group size, e.g. here there are 3 treated samples.
# If there are not discrete groups, one can use the minimal number of samples where non-zero counts would be considered interesting. 
# One can also omit this step entirely and just rely on the independent filtering procedures available in results(), either IHW or genefilter https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilt. 
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds2 <- dds[keep,]
################################################################################

# Run DESeq2 Analysis----------------------------
################################################################################
# Now you can run the DESeq2 analysis: Results tables are generated using the function results, which extracts a results table with log2 fold changes, p values and adjusted p values. 
# With no additional arguments to results, the log2 fold change and Wald test p value will be for the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable over the reference level (see previous note on factor levels). 
# However, the order of the variables of the design do not matter so long as the user specifies the comparison to build a results table for, using the name or contrast arguments of results.
dds_res_d0 <- DESeq(dds2)
# resultsNames(dds_res_d0)
################################################################################


################################################ Water type and Time (Filter time 6)


# Load the data: Compartment (Sediments), Sample_type (No_plant and Carex) and Water_type (OSPW)----------------------------
################################################################################
# Load the large gene_MT abundance in TSV file
# data <- fread("merged_gene_abundance.tsv", sep = "\t")
data <- fread(input_file, sep = "\t")
data <- data[, lapply(.SD, function(x) if (is.numeric(x)) x + 1 else x)]

# Load the metadata
metadata <- read_tsv("metadata.tsv")
metadata_full <- read_xlsx("metadata_Root.xlsx")

# Merge metadata
metadata <-  metadata %>% dplyr::left_join(metadata_full, by=c("SampleID", "Time"), multiple = "all", relationship = "many-to-many")

# Set the first column as row names in data
data <- data %>%
  column_to_rownames(var = colnames(data)[1])  # Set the first column as row names

# Set the first column as row names in metadata
metadata <- metadata %>% filter(Time == "D6") %>%
  column_to_rownames(var = colnames(metadata)[1])  # Set the first column as row names


# Find mismatches
mismatched_columns <- setdiff(colnames(data), rownames(metadata))
mismatched_rows <- setdiff(rownames(metadata), colnames(data))

# Print mismatches
if (length(mismatched_columns) > 0) {
  message("Columns in data not found in metadata row names: ", paste(mismatched_columns, collapse = ", "))
}

if (length(mismatched_rows) > 0) {
  message("Rows in metadata not found in data column names: ", paste(mismatched_rows, collapse = ", "))
}

# Subset data to keep only columns present in metadata row names
data <- data[, colnames(data) %in% rownames(metadata)]

# Reorder columns in data to match row names in metadata
data <- data[, match(rownames(metadata), colnames(data))]

# Set column names in data to match row names in metadata
colnames(data) <- rownames(metadata)

# Verify alignment
identical(colnames(data), rownames(metadata)) # Should return TRUE
################################################################################



# Construct DESEQDataSet Object---------------------------- 
################################################################################
# Make sure the metadata has the right structure and factor levels
# Ensure that Time and Sample_type are factors
# metadata$Time <- factor(metadata$Time, levels = c("D-0", "D6", "D19", "D40", "D60"))
metadata$Water_type <- factor(metadata$Water_type, levels = c("RO", "OSPW"))
# metadata$Greenhouse <- factor(metadata$Greenhouse, levels = c("E", "G"))


# Create the DESeq2 Dataset: In the design formula ~ Sample_type * Time, the DESeq2 model is structured to analyze the interaction effects of Sample_type and Time.
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadata,
                              design = ~ Water_type)
# The interaction terms indicate how the effect of time differs between the "Root" and "OSPW" sample types. 
# These terms are critical for understanding whether the presence of Root alters the time-dependent expression changes relative to the "OSPW" condition.
################################################################################


# Pre-filtering----------------------------
#################################################################################
# Here we perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples. 
# The count of 10 is a reasonable choice for bulk RNA-seq. A recommendation for the minimal number of samples is to specify the smallest group size, e.g. here there are 3 treated samples.
# If there are not discrete groups, one can use the minimal number of samples where non-zero counts would be considered interesting. 
# One can also omit this step entirely and just rely on the independent filtering procedures available in results(), either IHW or genefilter https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilt. 
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds2 <- dds[keep,]
################################################################################

# Run DESeq2 Analysis----------------------------
################################################################################
# Now you can run the DESeq2 analysis: Results tables are generated using the function results, which extracts a results table with log2 fold changes, p values and adjusted p values. 
# With no additional arguments to results, the log2 fold change and Wald test p value will be for the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable over the reference level (see previous note on factor levels). 
# However, the order of the variables of the design do not matter so long as the user specifies the comparison to build a results table for, using the name or contrast arguments of results.
dds_res_d6 <- DESeq(dds2)
# resultsNames(dds_res_d6)
################################################################################







################################################ Water type and Time (Filter time 19)

# Load the data: Compartment (Sediments), Sample_type (No_plant and Carex) and Water_type (OSPW)----------------------------
################################################################################
# Load the large gene_MT abundance in TSV file
# data <- fread("merged_gene_abundance.tsv", sep = "\t")
data <- fread(input_file, sep = "\t")
data <- data[, lapply(.SD, function(x) if (is.numeric(x)) x + 1 else x)]

# Load the metadata
metadata <- read_tsv("metadata.tsv")
metadata_full <- read_xlsx("metadata_Root.xlsx")

# Merge metadata
metadata <-  metadata %>% dplyr::left_join(metadata_full, by=c("SampleID", "Time"), multiple = "all", relationship = "many-to-many")

# Set the first column as row names in data
data <- data %>%
  column_to_rownames(var = colnames(data)[1])  # Set the first column as row names

# Set the first column as row names in metadata
metadata <- metadata %>% filter(Time == "D19") %>%
  column_to_rownames(var = colnames(metadata)[1])  # Set the first column as row names


# Find mismatches
mismatched_columns <- setdiff(colnames(data), rownames(metadata))
mismatched_rows <- setdiff(rownames(metadata), colnames(data))

# Print mismatches
if (length(mismatched_columns) > 0) {
  message("Columns in data not found in metadata row names: ", paste(mismatched_columns, collapse = ", "))
}

if (length(mismatched_rows) > 0) {
  message("Rows in metadata not found in data column names: ", paste(mismatched_rows, collapse = ", "))
}

# Subset data to keep only columns present in metadata row names
data <- data[, colnames(data) %in% rownames(metadata)]

# Reorder columns in data to match row names in metadata
data <- data[, match(rownames(metadata), colnames(data))]

# Set column names in data to match row names in metadata
colnames(data) <- rownames(metadata)

# Verify alignment
identical(colnames(data), rownames(metadata)) # Should return TRUE
################################################################################



# Construct DESEQDataSet Object---------------------------- 
################################################################################
# Make sure the metadata has the right structure and factor levels
# Ensure that Time and Sample_type are factors
# metadata$Time <- factor(metadata$Time, levels = c("D-0", "D6", "D19", "D40", "D60"))
metadata$Water_type <- factor(metadata$Water_type, levels = c("RO", "OSPW"))
# metadata$Greenhouse <- factor(metadata$Greenhouse, levels = c("E", "G"))


# Create the DESeq2 Dataset: In the design formula ~ Sample_type * Time, the DESeq2 model is structured to analyze the interaction effects of Sample_type and Time.
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadata,
                              design = ~ Water_type)
# The interaction terms indicate how the effect of time differs between the "Root" and "OSPW" sample types. 
# These terms are critical for understanding whether the presence of Root alters the time-dependent expression changes relative to the "OSPW" condition.
################################################################################


# Pre-filtering----------------------------
#################################################################################
# Here we perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples. 
# The count of 10 is a reasonable choice for bulk RNA-seq. A recommendation for the minimal number of samples is to specify the smallest group size, e.g. here there are 3 treated samples.
# If there are not discrete groups, one can use the minimal number of samples where non-zero counts would be considered interesting. 
# One can also omit this step entirely and just rely on the independent filtering procedures available in results(), either IHW or genefilter https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilt. 
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds2 <- dds[keep,]
################################################################################

# Run DESeq2 Analysis----------------------------
################################################################################
# Now you can run the DESeq2 analysis: Results tables are generated using the function results, which extracts a results table with log2 fold changes, p values and adjusted p values. 
# With no additional arguments to results, the log2 fold change and Wald test p value will be for the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable over the reference level (see previous note on factor levels). 
# However, the order of the variables of the design do not matter so long as the user specifies the comparison to build a results table for, using the name or contrast arguments of results.
dds_res_d19 <- DESeq(dds2)
# resultsNames(dds_res)




################################################ Water type and Time (Filter time 40)

# Load the data: Compartment (Sediments), Sample_type (No_plant and Carex) and Water_type (OSPW)----------------------------
################################################################################
# Load the large gene_MT abundance in TSV file
# data <- fread("merged_gene_abundance.tsv", sep = "\t")
data <- fread(input_file, sep = "\t")
data <- data[, lapply(.SD, function(x) if (is.numeric(x)) x + 1 else x)]

# Load the metadata
metadata <- read_tsv("metadata.tsv")
metadata_full <- read_xlsx("metadata_Root.xlsx")

# Merge metadata
metadata <-  metadata %>% dplyr::left_join(metadata_full, by=c("SampleID", "Time"), multiple = "all", relationship = "many-to-many")

# Set the first column as row names in data
data <- data %>%
  column_to_rownames(var = colnames(data)[1])  # Set the first column as row names

# Set the first column as row names in metadata
metadata <- metadata %>% filter(Time == "D40") %>%
  column_to_rownames(var = colnames(metadata)[1])  # Set the first column as row names


# Find mismatches
mismatched_columns <- setdiff(colnames(data), rownames(metadata))
mismatched_rows <- setdiff(rownames(metadata), colnames(data))

# Print mismatches
if (length(mismatched_columns) > 0) {
  message("Columns in data not found in metadata row names: ", paste(mismatched_columns, collapse = ", "))
}

if (length(mismatched_rows) > 0) {
  message("Rows in metadata not found in data column names: ", paste(mismatched_rows, collapse = ", "))
}

# Subset data to keep only columns present in metadata row names
data <- data[, colnames(data) %in% rownames(metadata)]

# Reorder columns in data to match row names in metadata
data <- data[, match(rownames(metadata), colnames(data))]

# Set column names in data to match row names in metadata
colnames(data) <- rownames(metadata)

# Verify alignment
identical(colnames(data), rownames(metadata)) # Should return TRUE
################################################################################



# Construct DESEQDataSet Object---------------------------- 
################################################################################
# Make sure the metadata has the right structure and factor levels
# Ensure that Time and Sample_type are factors
# metadata$Time <- factor(metadata$Time, levels = c("D-0", "D6", "D19", "D40", "D60"))
metadata$Water_type <- factor(metadata$Water_type, levels = c("RO", "OSPW"))
# metadata$Greenhouse <- factor(metadata$Greenhouse, levels = c("E", "G"))


# Create the DESeq2 Dataset: In the design formula ~ Sample_type * Time, the DESeq2 model is structured to analyze the interaction effects of Sample_type and Time.
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadata,
                              design = ~ Water_type)
# The interaction terms indicate how the effect of time differs between the "Root" and "OSPW" sample types. 
# These terms are critical for understanding whether the presence of Root alters the time-dependent expression changes relative to the "OSPW" condition.
################################################################################


# Pre-filtering----------------------------
#################################################################################
# Here we perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples. 
# The count of 10 is a reasonable choice for bulk RNA-seq. A recommendation for the minimal number of samples is to specify the smallest group size, e.g. here there are 3 treated samples.
# If there are not discrete groups, one can use the minimal number of samples where non-zero counts would be considered interesting. 
# One can also omit this step entirely and just rely on the independent filtering procedures available in results(), either IHW or genefilter https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilt. 
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds2 <- dds[keep,]
################################################################################

# Run DESeq2 Analysis----------------------------
################################################################################
# Now you can run the DESeq2 analysis: Results tables are generated using the function results, which extracts a results table with log2 fold changes, p values and adjusted p values. 
# With no additional arguments to results, the log2 fold change and Wald test p value will be for the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable over the reference level (see previous note on factor levels). 
# However, the order of the variables of the design do not matter so long as the user specifies the comparison to build a results table for, using the name or contrast arguments of results.
dds_res_d40 <- DESeq(dds2)
# resultsNames(dds_res)





################################################ Water type and Time (Filter time 60)

# Load the data: Compartment (Sediments), Sample_type (No_plant and Carex) and Water_type (OSPW)----------------------------
################################################################################
# Load the large gene_MT abundance in TSV file
# data <- fread("merged_gene_abundance.tsv", sep = "\t")
data <- fread(input_file, sep = "\t")
data <- data[, lapply(.SD, function(x) if (is.numeric(x)) x + 1 else x)]

# Load the metadata
metadata <- read_tsv("metadata.tsv")
metadata_full <- read_xlsx("metadata_Root.xlsx")

# Merge metadata
metadata <-  metadata %>% dplyr::left_join(metadata_full, by=c("SampleID", "Time"), multiple = "all", relationship = "many-to-many")

# Set the first column as row names in data
data <- data %>%
  column_to_rownames(var = colnames(data)[1])  # Set the first column as row names

# Set the first column as row names in metadata
metadata <- metadata %>% filter(Time == "D60") %>%
  column_to_rownames(var = colnames(metadata)[1])  # Set the first column as row names


# Find mismatches
mismatched_columns <- setdiff(colnames(data), rownames(metadata))
mismatched_rows <- setdiff(rownames(metadata), colnames(data))

# Print mismatches
if (length(mismatched_columns) > 0) {
  message("Columns in data not found in metadata row names: ", paste(mismatched_columns, collapse = ", "))
}

if (length(mismatched_rows) > 0) {
  message("Rows in metadata not found in data column names: ", paste(mismatched_rows, collapse = ", "))
}

# Subset data to keep only columns present in metadata row names
data <- data[, colnames(data) %in% rownames(metadata)]

# Reorder columns in data to match row names in metadata
data <- data[, match(rownames(metadata), colnames(data))]

# Set column names in data to match row names in metadata
colnames(data) <- rownames(metadata)

# Verify alignment
identical(colnames(data), rownames(metadata)) # Should return TRUE
################################################################################



# Construct DESEQDataSet Object---------------------------- 
################################################################################
# Make sure the metadata has the right structure and factor levels
# Ensure that Time and Sample_type are factors
# metadata$Time <- factor(metadata$Time, levels = c("D-0", "D6", "D19", "D40", "D60"))
metadata$Water_type <- factor(metadata$Water_type, levels = c("RO", "OSPW"))
# metadata$Greenhouse <- factor(metadata$Greenhouse, levels = c("E", "G"))


# Create the DESeq2 Dataset: In the design formula ~ Sample_type * Time, the DESeq2 model is structured to analyze the interaction effects of Sample_type and Time.
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadata,
                              design = ~ Water_type)
# The interaction terms indicate how the effect of time differs between the "Root" and "OSPW" sample types. 
# These terms are critical for understanding whether the presence of Root alters the time-dependent expression changes relative to the "OSPW" condition.
################################################################################


# Pre-filtering----------------------------
#################################################################################
# Here we perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples. 
# The count of 10 is a reasonable choice for bulk RNA-seq. A recommendation for the minimal number of samples is to specify the smallest group size, e.g. here there are 3 treated samples.
# If there are not discrete groups, one can use the minimal number of samples where non-zero counts would be considered interesting. 
# One can also omit this step entirely and just rely on the independent filtering procedures available in results(), either IHW or genefilter https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilt. 
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds2 <- dds[keep,]
################################################################################

# Run DESeq2 Analysis----------------------------
################################################################################
# Now you can run the DESeq2 analysis: Results tables are generated using the function results, which extracts a results table with log2 fold changes, p values and adjusted p values. 
# With no additional arguments to results, the log2 fold change and Wald test p value will be for the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable over the reference level (see previous note on factor levels). 
# However, the order of the variables of the design do not matter so long as the user specifies the comparison to build a results table for, using the name or contrast arguments of results.
dds_res_d60 <- DESeq(dds2)
# resultsNames(dds_res)








# Extract Results for Each Comparison----------------------------
################################################################################

results_interaction_d0 <- results(dds_res_d0, name = "Water_type_OSPW_vs_RO")
results_interaction_d6 <- results(dds_res_d6, name = "Water_type_OSPW_vs_RO")
results_interaction_d19 <- results(dds_res_d19, name = "Water_type_OSPW_vs_RO")
results_interaction_d40 <- results(dds_res_d40, name = "Water_type_OSPW_vs_RO")
results_interaction_d60 <- results(dds_res_d60, name = "Water_type_OSPW_vs_RO")


# Define your threshold values for significance
padj_threshold <- 0.05
log2fc_threshold <- 1

# Convert to a dataframe with filtering the DEGs for D0
Unfiltered_degs_d0 <- as.data.frame(results_interaction_d0)
Unfiltered_degs_d0$Gene <- rownames(Unfiltered_degs_d0) # Add the Gene column from the row names
Unfiltered_degs_d0$Time <- "D0" # Add a Time column to each filtered DEGs dataframe

# Convert to a dataframe with filtering the DEGs for D6
Unfiltered_degs_d6 <- as.data.frame(results_interaction_d6)
Unfiltered_degs_d6$Gene <- rownames(Unfiltered_degs_d6) # Add the Gene column from the row names
Unfiltered_degs_d6$Time <- "D6" # Add a Time column to each filtered DEGs dataframe

# Convert to a dataframe with filtering the DEGs for D19
Unfiltered_degs_d19 <- as.data.frame(results_interaction_d19)
Unfiltered_degs_d19$Gene <- rownames(Unfiltered_degs_d19) # Add the Gene column from the row names
Unfiltered_degs_d19$Time <- "D19" # Add a Time column to each filtered DEGs dataframe

# Convert to a dataframe with filtering the DEGs for D40
Unfiltered_degs_d40 <- as.data.frame(results_interaction_d40)
Unfiltered_degs_d40$Gene <- rownames(Unfiltered_degs_d40) # Add the Gene column from the row names
Unfiltered_degs_d40$Time <- "D40" # Add a Time column to each filtered DEGs dataframe

# Convert to a dataframe with filtering the DEGs for D60
Unfiltered_degs_d60 <- as.data.frame(results_interaction_d60)
Unfiltered_degs_d60$Gene <- rownames(Unfiltered_degs_d60) # Add the Gene column from the row names
Unfiltered_degs_d60$Time <- "D60" # Add a Time column to each filtered DEGs dataframe

# Merge the data frames:
Unfiltered_merged_degs <- rbind(Unfiltered_degs_d0, Unfiltered_degs_d6, Unfiltered_degs_d19, Unfiltered_degs_d40, Unfiltered_degs_d60)
Unfiltered_merged_degs$SeqType <- "Metatranscriptome" 
Unfiltered_merged_degs$Treatment <- "Water_type_OSPW_vs_RO_at_each_Time" 

# Filter for significant DEGs and rearrange columns
Filtered_merged_degs <- subset(Unfiltered_merged_degs, padj < padj_threshold & abs(log2FoldChange) > log2fc_threshold)
Filtered_merged_degs <- Filtered_merged_degs[, c("SeqType", "Gene", "Time", setdiff(colnames(Filtered_merged_degs), c("SeqType", "Gene", "Time")))]

# Save the merged DEGs dataframe to a CSV
# write.csv(Filtered_merged_degs, "Filtered_merged_degs_Water_type_OSPW_vs_RO_NoTimeIN.csv")
write.csv(Filtered_merged_degs, output_file, row.names = FALSE)
