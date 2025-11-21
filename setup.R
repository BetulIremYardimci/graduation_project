# setup.R - TCGA Histone Analysis Paketleri

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Bioconductor paketleri
bioc_packages <- c(
  "TCGAbiolinks",           # TCGA data download
  "SummarizedExperiment",   # Data structures
  "DESeq2",                 # Differential expression
  "edgeR",                  # Alternative DE analysis
  "limma",                  # Linear models
  "biomaRt",                # Gene annotations
  "ComplexHeatmap"          # Heatmaps
)

# CRAN paketleri
cran_packages <- c(
  "survival",               # Survival analysis
  "survminer",              # Survival plots
  "ggplot2",                # Visualizations
  "dplyr",                  # Data manipulation
  "tidyr",                  # Data tidying
  "yaml",                   # YAML config reader
  "pheatmap",               # Heatmaps
  "ggpubr"                  # Publication-ready plots
)

# Bioconductor
BiocManager::install(bioc_packages, update = TRUE, ask = FALSE)

# CRAN
install.packages(cran_packages)

cat("\n✓ All packages downloaded!\n")