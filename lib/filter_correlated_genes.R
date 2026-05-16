# ============================================================================
# PRACTICAL SCRIPT: Filter Your Correlated Genes (FIXED PDF VERSION)
# ============================================================================
# Use this AFTER running main_single_cancer.R
# This is a SEPARATE, STANDALONE script
# ============================================================================

library(dplyr)

cat("=================================================================\n")
cat("FILTERING CORRELATED GENES\n")
cat("=================================================================\n\n")

# ============================================================================
# CONFIGURATION - CHANGE THESE FOR YOUR CANCER
# ============================================================================

CANCER <- "LUAD"  # Change this to match your current cancer
RESULTS_DIR <- paste0("results/", CANCER)

# Check if results directory exists
if (!dir.exists(RESULTS_DIR)) {
  stop("ERROR: Results directory not found: ", RESULTS_DIR, "\n",
       "       Make sure you've run main_single_cancer.R first!")
}

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading data from:", RESULTS_DIR, "\n")

# Check required files
required_files <- c(
  paste0(RESULTS_DIR, "/correlated_genes_summary.csv"),
  paste0(RESULTS_DIR, "/correlation_results.rds"),
  paste0(RESULTS_DIR, "/deseq2_full_results.csv")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("ERROR: Missing required files:\n  ",
       paste(missing_files, collapse = "\n  "))
}

# Load data
cat("  Loading correlated_genes_summary.csv...\n")
correlated_genes <- read.csv(paste0(RESULTS_DIR, "/correlated_genes_summary.csv"))

cat("  Loading correlation_results.rds...\n")
correlation_results <- readRDS(paste0(RESULTS_DIR, "/correlation_results.rds"))

cat("  Loading deseq2_full_results.csv...\n")
deseq2_full <- read.csv(paste0(RESULTS_DIR, "/deseq2_full_results.csv"))

cat("\nData loaded successfully!\n")
cat("Original correlated genes:", nrow(correlated_genes), "\n")
cat("Full transcriptome genes:", nrow(deseq2_full), "\n\n")


# ============================================================================
# FILTERING STRATEGIES
# ============================================================================

cat("=================================================================\n")
cat("APPLYING FILTERING STRATEGIES\n")
cat("=================================================================\n\n")

# --------------------------------------------------------------------
# OPTION 1: Stricter Correlation Threshold
# --------------------------------------------------------------------
cat("OPTION 1: Stricter Correlation (|r| > 0.5)\n")
cat("--------------------------------------------\n")

filtered_1 <- correlated_genes %>%
  filter(max_correlation > 0.5)

cat("Result:", nrow(filtered_1), "genes\n")
cat("Reduction:", round((1 - nrow(filtered_1)/nrow(correlated_genes)) * 100, 1), "%\n\n")


# --------------------------------------------------------------------
# OPTION 2: Multiple Histone Correlations
# --------------------------------------------------------------------
cat("OPTION 2: Multiple Histone Correlations (2+ histones)\n")
cat("------------------------------------------------------\n")

filtered_2 <- correlated_genes %>%
  filter(n_histone_correlated >= 2)

cat("Result:", nrow(filtered_2), "genes\n")
cat("Reduction:", round((1 - nrow(filtered_2)/nrow(correlated_genes)) * 100, 1), "%\n\n")


# --------------------------------------------------------------------
# OPTION 3: Combined Strict
# --------------------------------------------------------------------
cat("OPTION 3: Combined (Correlation + Multiple Histones)\n")
cat("-----------------------------------------------------\n")

filtered_3 <- correlated_genes %>%
  filter(
    max_correlation > 0.5,
    n_histone_correlated >= 2
  )

cat("Result:", nrow(filtered_3), "genes\n")
cat("Reduction:", round((1 - nrow(filtered_3)/nrow(correlated_genes)) * 100, 1), "%\n\n")


# --------------------------------------------------------------------
# OPTION 4a: Also Differentially Expressed
# --------------------------------------------------------------------
cat("OPTION 4a: Correlated + DE (padj < 0.05)\n")
cat("-----------------------------------------\n")

correlated_with_DE <- correlated_genes %>%
  inner_join(
    deseq2_full %>%
      select(ensembl_id, baseMean, log2FoldChange, padj),
    by = c("gene_id" = "ensembl_id")
  )

filtered_4a <- correlated_with_DE %>%
  filter(
    padj < 0.05,
    baseMean > 50
  )

cat("Result:", nrow(filtered_4a), "genes\n")
cat("Reduction:", round((1 - nrow(filtered_4a)/nrow(correlated_genes)) * 100, 1), "%\n\n")


# --------------------------------------------------------------------
# OPTION 4b: HIGH CONFIDENCE (RECOMMENDED!)
# --------------------------------------------------------------------
cat("OPTION 4b: HIGH CONFIDENCE ⭐ (RECOMMENDED)\n")
cat("-------------------------------------------\n")

filtered_4b <- correlated_with_DE %>%
  filter(
    max_correlation > 0.5,
    padj < 0.05,
    baseMean > 100,
    abs(log2FoldChange) > 1
  )

cat("Result:", nrow(filtered_4b), "genes\n")
cat("Reduction:", round((1 - nrow(filtered_4b)/nrow(correlated_genes)) * 100, 1), "%\n\n")


# --------------------------------------------------------------------
# OPTION 5: Top N per Histone
# --------------------------------------------------------------------
cat("OPTION 5: Top 50 per Histone Gene\n")
cat("----------------------------------\n")

top_per_histone <- data.frame()

for (histone_id in names(correlation_results$correlations)) {

  sig_genes <- correlation_results$correlations[[histone_id]]$significant

  if (nrow(sig_genes) > 0) {
    top_50 <- sig_genes %>%
      arrange(desc(abs(correlation))) %>%
      head(50) %>%
      mutate(histone_gene = histone_id)

    top_per_histone <- rbind(top_per_histone, top_50)
  }
}

top_unique <- top_per_histone %>%
  group_by(gene_id) %>%
  slice_max(abs(correlation), n = 1) %>%
  ungroup()

cat("Result:", nrow(top_unique), "genes\n")
cat("Reduction:", round((1 - nrow(top_unique)/nrow(correlated_genes)) * 100, 1), "%\n\n")


# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("=================================================================\n")
cat("SAVING FILTERED RESULTS\n")
cat("=================================================================\n\n")

write.csv(filtered_1,
          paste0(RESULTS_DIR, "/correlated_genes_filtered_r0.5.csv"),
          row.names = FALSE)
cat("✓ Saved: correlated_genes_filtered_r0.5.csv\n")

write.csv(filtered_3,
          paste0(RESULTS_DIR, "/correlated_genes_filtered_strict.csv"),
          row.names = FALSE)
cat("✓ Saved: correlated_genes_filtered_strict.csv\n")

write.csv(filtered_4a,
          paste0(RESULTS_DIR, "/correlated_genes_filtered_DE.csv"),
          row.names = FALSE)
cat("✓ Saved: correlated_genes_filtered_DE.csv\n")

write.csv(filtered_4b,
          paste0(RESULTS_DIR, "/correlated_genes_filtered_HIGH_CONFIDENCE.csv"),
          row.names = FALSE)
cat("✓ Saved: correlated_genes_filtered_HIGH_CONFIDENCE.csv\n")

write.csv(top_unique,
          paste0(RESULTS_DIR, "/correlated_genes_filtered_top50.csv"),
          row.names = FALSE)
cat("✓ Saved: correlated_genes_filtered_top50.csv\n\n")


# ============================================================================
# SUMMARY TABLE
# ============================================================================

cat("=================================================================\n")
cat("FILTERING SUMMARY\n")
cat("=================================================================\n\n")

summary_table <- data.frame(
  Filter = c(
    "Original",
    "Option 1: |r| > 0.5",
    "Option 2: 2+ histones",
    "Option 3: Combined",
    "Option 4a: + DE",
    "Option 4b: HIGH CONFIDENCE",
    "Option 5: Top 50/histone"
  ),
  N_Genes = c(
    nrow(correlated_genes),
    nrow(filtered_1),
    nrow(filtered_2),
    nrow(filtered_3),
    nrow(filtered_4a),
    nrow(filtered_4b),
    nrow(top_unique)
  ),
  Reduction = paste0(
    round(c(
      0,
      (1 - nrow(filtered_1)/nrow(correlated_genes)) * 100,
      (1 - nrow(filtered_2)/nrow(correlated_genes)) * 100,
      (1 - nrow(filtered_3)/nrow(correlated_genes)) * 100,
      (1 - nrow(filtered_4a)/nrow(correlated_genes)) * 100,
      (1 - nrow(filtered_4b)/nrow(correlated_genes)) * 100,
      (1 - nrow(top_unique)/nrow(correlated_genes)) * 100
    ), 1), "%"
  ),
  Recommendation = c(
    "Baseline",
    "Moderate",
    "Moderate",
    "Good",
    "Very Good",
    "Best ⭐",
    "Alternative"
  )
)

print(summary_table)

# Save summary
write.csv(summary_table,
          paste0(RESULTS_DIR, "/filtering_summary.csv"),
          row.names = FALSE)


# ============================================================================
# VISUALIZATION - FIXED VERSION
# ============================================================================

cat("\n=================================================================\n")
cat("CREATING COMPARISON PLOT\n")
cat("=================================================================\n\n")

# Load ggplot2
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  cat("Warning: ggplot2 not installed. Skipping plot.\n")
} else {

  library(ggplot2)

  # Prepare data for plotting
  plot_data <- summary_table[-1, ]  # Remove "Original" row

  tryCatch({

    # Create plot
    p <- ggplot(plot_data, aes(x = reorder(Filter, -N_Genes), y = N_Genes)) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
      geom_text(aes(label = N_Genes), vjust = -0.5, size = 4) +
      geom_hline(yintercept = 500, linetype = "dashed",
                 color = "red", alpha = 0.5) +
      annotate("text", x = 1, y = 550,
               label = "Target: ~500 genes",
               color = "red", size = 3) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)
      ) +
      labs(
        title = "Correlated Genes: Filtering Strategies",
        subtitle = paste0(CANCER, " - Original: ",
                          nrow(correlated_genes), " genes"),
        x = "Filtering Strategy",
        y = "Number of Genes"
      )

    # Save plot with explicit device settings
    pdf_file <- paste0(RESULTS_DIR, "/filtering_comparison.pdf")

    pdf(pdf_file, width = 10, height = 6)
    print(p)
    dev.off()

    cat("✓ Plot saved:", pdf_file, "\n")

    # Also save as PNG (more compatible)
    png_file <- paste0(RESULTS_DIR, "/filtering_comparison.png")
    ggsave(png_file, plot = p, width = 10, height = 6, dpi = 300)
    cat("✓ Plot also saved as PNG:", png_file, "\n")

  }, error = function(e) {
    cat("Warning: Could not create plot\n")
    cat("Error:", e$message, "\n")
    cat("This is OK - filtered CSV files are still saved!\n")
  })
}


# ============================================================================
# FINAL RECOMMENDATION
# ============================================================================

cat("\n=================================================================\n")
cat("RECOMMENDATION\n")
cat("=================================================================\n\n")

cat("Use: correlated_genes_filtered_HIGH_CONFIDENCE.csv\n\n")

cat("This file contains", nrow(filtered_4b), "genes that meet ALL criteria:\n")
cat("  ✓ Strong correlation with histone (|r| > 0.5)\n")
cat("  ✓ Differentially expressed (padj < 0.05)\n")
cat("  ✓ Well expressed (baseMean > 100)\n")
cat("  ✓ Substantial change (|log2FC| > 1)\n\n")

cat("Perfect for:\n")
cat("  - Pathway enrichment analysis\n")
cat("  - Network analysis\n")
cat("  - Clinical associations\n")
cat("  - Publication figures\n\n")

cat("=================================================================\n")
cat("FILTERING COMPLETE!\n")
cat("=================================================================\n")

cat("\nFiles saved in:", RESULTS_DIR, "/\n")
cat("\nNext steps:\n")
cat("  1. Review: correlated_genes_filtered_HIGH_CONFIDENCE.csv\n")
cat("  2. Run pathway analysis with these genes\n")
cat("  3. Create network visualizations\n")
cat("  4. Proceed to clinical integration\n\n")
