# ============================================================================
# TCGA Pan-Cancer Histone Variant Analysis - Package Setup
# ============================================================================
# This script installs all required R packages for the analysis pipeline
# Run this once before starting the analysis
# ============================================================================

cat("\n")
cat("========================================================================\n")
cat("  TCGA Pan-Cancer Analysis - Package Installation\n")
cat("========================================================================\n\n")

# Check R version
r_version <- paste(R.version$major, R.version$minor, sep = ".")
cat("R version:", r_version, "\n")

if (as.numeric(R.version$major) < 4 || 
    (as.numeric(R.version$major) == 4 && as.numeric(R.version$minor) < 3)) {
  warning("R version 4.3.0 or higher is recommended!")
}

# ============================================================================
# INSTALL BiocManager (if not already installed)
# ============================================================================

cat("\n[1/3] Installing BiocManager...\n")

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  cat("  ✓ BiocManager installed\n")
} else {
  cat("  ✓ BiocManager already installed\n")
}

library(BiocManager)

# ============================================================================
# INSTALL Bioconductor Packages
# ============================================================================

cat("\n[2/3] Installing Bioconductor packages...\n")
cat("  (This may take 10-15 minutes)\n\n")

bioc_packages <- c(
  "TCGAbiolinks",           # TCGA data download and processing
  "SummarizedExperiment",   # Data container for genomics
  "DESeq2",                 # Differential expression analysis
  "edgeR",                  # Alternative DE analysis
  "limma",                  # Linear models for microarray/RNA-seq
  "biomaRt",                # Gene annotation
  "ComplexHeatmap",         # Advanced heatmaps
  "sva"                     # Batch effect correction (ComBat)
)

for (pkg in bioc_packages) {
  cat("  Installing", pkg, "...\n")
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

cat("\n  ✓ All Bioconductor packages installed\n")

# ============================================================================
# INSTALL CRAN Packages
# ============================================================================

cat("\n[3/3] Installing CRAN packages...\n\n")

cran_packages <- c(
  "survival",               # Survival analysis
  "survminer",              # Survival plots
  "ggplot2",                # Data visualization
  "dplyr",                  # Data manipulation
  "tidyr",                  # Data tidying
  "yaml",                   # YAML configuration files
  "pheatmap",               # Pretty heatmaps
  "ggpubr",                 # Publication-ready plots
  "cowplot",                # Multi-panel figures
  "RColorBrewer",           # Color palettes
  "viridis",                # Color scales
  "reshape2",               # Data reshaping
  "stringr",                # String operations
  "readr"                   # Fast data reading
)

for (pkg in cran_packages) {
  cat("  Installing", pkg, "...\n")
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, quiet = TRUE)
  }
}

cat("\n  ✓ All CRAN packages installed\n")

# ============================================================================
# VERIFY Installation
# ============================================================================

cat("\n")
cat("========================================================================\n")
cat("  Verifying Installation\n")
cat("========================================================================\n\n")

all_packages <- c(bioc_packages, cran_packages)
installation_status <- sapply(all_packages, function(pkg) {
  require(pkg, character.only = TRUE, quietly = TRUE)
})

if (all(installation_status)) {
  cat("✓ SUCCESS: All packages installed correctly!\n\n")
  cat("Installed packages:\n")
  for (pkg in all_packages) {
    version <- packageVersion(pkg)
    cat("  •", pkg, "-", as.character(version), "\n")
  }
} else {
  failed_packages <- names(installation_status[!installation_status])
  cat("✗ WARNING: Some packages failed to install:\n")
  for (pkg in failed_packages) {
    cat("  •", pkg, "\n")
  }
  cat("\nPlease install these packages manually:\n")
  for (pkg in failed_packages) {
    if (pkg %in% bioc_packages) {
      cat('  BiocManager::install("', pkg, '")\n', sep = "")
    } else {
      cat('  install.packages("', pkg, '")\n', sep = "")
    }
  }
}

# ============================================================================
# SAVE Session Info
# ============================================================================

cat("\n")
cat("========================================================================\n")
cat("  Saving Session Information\n")
cat("========================================================================\n\n")

sink("session_info.txt")
cat("TCGA Pan-Cancer Analysis - Package Installation\n")
cat("Date:", as.character(Sys.time()), "\n\n")
cat("========================================================================\n\n")
sessionInfo()
sink()

cat("Session information saved to: session_info.txt\n")

# ============================================================================
# INITIALIZE renv (Optional but Recommended)
# ============================================================================

cat("\n")
cat("========================================================================\n")
cat("  Optional: Initialize renv for Reproducibility\n")
cat("========================================================================\n\n")

cat("Would you like to initialize renv for package version control?\n")
cat("This creates an isolated, reproducible package environment.\n")
cat("(Recommended for publication-quality research)\n\n")

cat("To initialize renv later, run:\n")
cat('  install.packages("renv")\n')
cat('  renv::init()\n')
cat('  renv::snapshot()\n\n')

# ============================================================================
# NEXT STEPS
# ============================================================================

cat("\n")
cat("========================================================================\n")
cat("  Setup Complete! Next Steps:\n")
cat("========================================================================\n\n")

cat("1. Verify your config.yaml file has correct settings\n")
cat("2. Ensure you have stable internet connection\n")
cat("3. Check you have sufficient storage (~50 GB for 15 cancers)\n")
cat("4. Run the main analysis:\n")
cat('   source("histone_main.R")\n\n')

cat("For testing with small sample size:\n")
cat('   • Edit config.yaml: use_full_dataset: false\n')
cat('   • This will use only 100 samples per cancer\n\n')

cat("========================================================================\n")
cat("  Installation successful! Ready to analyze TCGA data.\n")
cat("========================================================================\n\n")

# ============================================================================
# SYSTEM INFORMATION
# ============================================================================

cat("System Information:\n")
cat("  OS:", Sys.info()["sysname"], Sys.info()["release"], "\n")
cat("  R version:", r_version, "\n")
cat("  Platform:", R.version$platform, "\n")
cat("  Memory:", round(memory.size()/1024, 1), "GB available\n")
cat("\n")