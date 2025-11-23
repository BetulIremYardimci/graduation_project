#!/usr/bin/env Rscript

library(TCGAbiolinks)

# Load functions
source("lib/tcga_dataset_extraction.R")
source("lib/differential_expression.R")
source("lib/visualization.R")

# 1. Download data
cat("Downloading TCGA data...\n")
query <- download_tcga_data(cancer_type = "BRCA")

# 2. Prepare data
cat("Preparing data...\n")
data <- prepare_tcga_data(query)

# 3. Run differential expression
cat("Running differential expression analysis...\n")
results <- run_differential_expression(data)

# 4. Visualize
cat("Creating plots...\n")
create_plots(data, results)

cat("Analysis complete!\n")