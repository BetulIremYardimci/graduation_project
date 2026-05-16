# ============================================================================
# Pan-Cancer Results Aggregation Script
# ============================================================================
# Purpose: Combine results from individual cancer analyses
# Run this AFTER you've analyzed all cancer types individually
# Memory: Light (~1-2 GB) - only reads result files
# ============================================================================

library(dplyr)
library(yaml)

cat("==================================================================\n")
cat("PAN-CANCER RESULTS AGGREGATION\n")
cat("==================================================================\n")

# Load configuration
config <- yaml::yaml.load_file("config.yaml")
all_cancer_types <- config$cancer_types

cat("Cancer types to aggregate:", length(all_cancer_types), "\n")
print(all_cancer_types)
cat("\n")

# ============================================================================
# Check which analyses are complete
# ============================================================================

completed_cancers <- c()
missing_cancers <- c()

for (cancer in all_cancer_types) {
  
  results_dir <- paste0("results/", cancer)
  consensus_file <- paste0(results_dir, "/consensus.csv")
  
  if (file.exists(consensus_file)) {
    completed_cancers <- c(completed_cancers, cancer)
  } else {
    missing_cancers <- c(missing_cancers, cancer)
  }
}

cat("Completed analyses:", length(completed_cancers), "\n")
if (length(completed_cancers) > 0) {
  cat("  ", paste(completed_cancers, collapse = ", "), "\n")
}

if (length(missing_cancers) > 0) {
  cat("\nMissing analyses:", length(missing_cancers), "\n")
  cat("  ", paste(missing_cancers, collapse = ", "), "\n")
  cat("\n⚠️  Run main_single_cancer.R for these cancer types first!\n")
}

if (length(completed_cancers) == 0) {
  stop("No completed analyses found. Run main_single_cancer.R first.")
}

cat("\nProceeding with", length(completed_cancers), "cancer types...\n")
cat("==================================================================\n")

# ============================================================================
# Aggregate Results
# ============================================================================

## --------------------------------------------------------------------
## 1. Pan-Cancer Summary Table
## --------------------------------------------------------------------
cat("\n>>> Creating Pan-Cancer Summary Table...\n")

summary_df <- data.frame()

for (cancer in completed_cancers) {
  
  cat("  Processing:", cancer, "\n")
  
  results_dir <- paste0("results/", cancer)
  
  # Load consensus
  consensus <- read.csv(paste0(results_dir, "/consensus.csv"))
  
  # Load metadata
  metadata <- readRDS(paste0(results_dir, "/analysis_metadata.rds"))
  
  # Create summary row
  summary_row <- data.frame(
    cancer_type = cancer,
    n_samples = metadata$n_samples,
    n_genes_analyzed = metadata$n_genes_analyzed,
    n_histone_genes = nrow(consensus),
    deseq2_sig = sum(consensus$deseq2_sig, na.rm = TRUE),
    edger_sig = sum(consensus$edger_sig, na.rm = TRUE),
    limma_sig = sum(consensus$limma_sig, na.rm = TRUE),
    consensus_sig = sum(consensus$consensus, na.rm = TRUE),
    consensus_upregulated = sum(consensus$consensus & 
                                consensus$deseq2_log2fc > 0, na.rm = TRUE),
    consensus_downregulated = sum(consensus$consensus & 
                                  consensus$deseq2_log2fc < 0, na.rm = TRUE),
    date_analyzed = as.character(metadata$date_analyzed)
  )
  
  summary_df <- rbind(summary_df, summary_row)
}

# Save summary
write.csv(summary_df, "results/pan_cancer_summary.csv", row.names = FALSE)
cat("\n✓ Summary table saved: results/pan_cancer_summary.csv\n")

print(summary_df)

## --------------------------------------------------------------------
## 2. Combined Histone Gene Results
## --------------------------------------------------------------------
cat("\n>>> Combining Histone Gene Results Across Cancers...\n")

combined_histone <- data.frame()

for (cancer in completed_cancers) {
  
  results_dir <- paste0("results/", cancer)
  consensus <- read.csv(paste0(results_dir, "/consensus.csv"))
  
  consensus$cancer_type <- cancer
  combined_histone <- rbind(combined_histone, consensus)
}

write.csv(combined_histone, "results/pan_cancer_histone_combined.csv", 
          row.names = FALSE)
cat("✓ Combined results saved: results/pan_cancer_histone_combined.csv\n")

## --------------------------------------------------------------------
## 3. Gene-Centric View (Each histone gene across cancers)
## --------------------------------------------------------------------
cat("\n>>> Creating Gene-Centric Summary...\n")

histone_genes <- unique(combined_histone$gene_name)

gene_centric <- data.frame()

for (gene in histone_genes) {
  
  gene_data <- combined_histone %>%
    filter(gene_name == gene)
  
  gene_row <- data.frame(
    gene_name = gene,
    ensembl_id = unique(gene_data$ensembl_id)[1],
    n_cancers_analyzed = nrow(gene_data),
    n_cancers_deseq2_sig = sum(gene_data$deseq2_sig, na.rm = TRUE),
    n_cancers_consensus = sum(gene_data$consensus, na.rm = TRUE),
    mean_log2fc_deseq2 = mean(gene_data$deseq2_log2fc, na.rm = TRUE),
    cancers_upregulated = paste(gene_data$cancer_type[
      gene_data$consensus & gene_data$deseq2_log2fc > 0
    ], collapse = ";"),
    cancers_downregulated = paste(gene_data$cancer_type[
      gene_data$consensus & gene_data$deseq2_log2fc < 0
    ], collapse = ";")
  )
  
  gene_centric <- rbind(gene_centric, gene_row)
}

write.csv(gene_centric, "results/pan_cancer_gene_centric_summary.csv", 
          row.names = FALSE)
cat("✓ Gene-centric summary saved: results/pan_cancer_gene_centric_summary.csv\n")

## --------------------------------------------------------------------
## 4. Correlation Summary
## --------------------------------------------------------------------
cat("\n>>> Aggregating Correlation Results...\n")

correlation_summary <- data.frame()

for (cancer in completed_cancers) {
  
  results_dir <- paste0("results/", cancer)
  corr_file <- paste0(results_dir, "/correlated_genes_summary.csv")
  
  if (file.exists(corr_file)) {
    
    corr_data <- read.csv(corr_file)
    
    corr_row <- data.frame(
      cancer_type = cancer,
      n_correlated_genes = nrow(corr_data),
      mean_n_histone = mean(corr_data$n_histone_correlated, na.rm = TRUE),
      max_n_histone = max(corr_data$n_histone_correlated, na.rm = TRUE),
      mean_correlation = mean(corr_data$mean_correlation, na.rm = TRUE)
    )
    
    correlation_summary <- rbind(correlation_summary, corr_row)
  }
}

if (nrow(correlation_summary) > 0) {
  write.csv(correlation_summary, 
            "results/pan_cancer_correlation_summary.csv", 
            row.names = FALSE)
  cat("✓ Correlation summary saved\n")
}

## --------------------------------------------------------------------
## 5. Create Comparison Matrices (for heatmaps)
## --------------------------------------------------------------------
cat("\n>>> Creating Comparison Matrices...\n")

# Matrix: Cancer x Histone Gene (Significance)
sig_matrix <- combined_histone %>%
  select(cancer_type, gene_name, consensus) %>%
  tidyr::pivot_wider(names_from = gene_name, values_from = consensus) %>%
  as.data.frame()

rownames(sig_matrix) <- sig_matrix$cancer_type
sig_matrix <- sig_matrix[, -1]

saveRDS(sig_matrix, "results/pan_cancer_significance_matrix.rds")

# Matrix: Cancer x Histone Gene (log2FC)
fc_matrix <- combined_histone %>%
  select(cancer_type, gene_name, deseq2_log2fc) %>%
  tidyr::pivot_wider(names_from = gene_name, values_from = deseq2_log2fc) %>%
  as.data.frame()

rownames(fc_matrix) <- fc_matrix$cancer_type
fc_matrix <- fc_matrix[, -1]

saveRDS(fc_matrix, "results/pan_cancer_log2fc_matrix.rds")

cat("✓ Comparison matrices saved\n")

## --------------------------------------------------------------------
## 6. Quick Heatmap Visualization
## --------------------------------------------------------------------
cat("\n>>> Creating Pan-Cancer Heatmap...\n")

library(pheatmap)

# Significance heatmap
pdf("results/pan_cancer_significance_heatmap.pdf", width = 8, height = 6)
pheatmap(as.matrix(sig_matrix),
         color = c("white", "red"),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Pan-Cancer Histone Gene Significance",
         display_numbers = TRUE,
         number_format = "%.0f",
         legend = FALSE)
dev.off()

# Log2FC heatmap
pdf("results/pan_cancer_log2fc_heatmap.pdf", width = 10, height = 8)
pheatmap(as.matrix(fc_matrix),
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-3, 3, length.out = 101),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Pan-Cancer Histone Gene Expression Changes",
         na_col = "grey90")
dev.off()

cat("✓ Heatmaps saved\n")

# ============================================================================
# Final Summary Report
# ============================================================================
cat("\n====================================================================\n")
cat("PAN-CANCER AGGREGATION COMPLETED\n")
cat("====================================================================\n")

summary_report <- sprintf("
Pan-Cancer Histone Variant Analysis Summary
===========================================

Date: %s

Cancer Types Analyzed: %d
  %s

Total Samples: %d
Average Genes Analyzed per Cancer: %.0f

Histone Genes Analyzed: %d

Consensus Significance Summary:
  Total significant findings: %d
  Average per cancer: %.1f
  Most significant cancers: %s

Files Generated:
  1. results/pan_cancer_summary.csv
  2. results/pan_cancer_histone_combined.csv
  3. results/pan_cancer_gene_centric_summary.csv
  4. results/pan_cancer_correlation_summary.csv
  5. results/pan_cancer_significance_matrix.rds
  6. results/pan_cancer_log2fc_matrix.rds
  7. results/pan_cancer_significance_heatmap.pdf
  8. results/pan_cancer_log2fc_heatmap.pdf

Next Steps:
  1. Review pan-cancer heatmaps
  2. Identify cancer-specific patterns
  3. Proceed to Phase 2: Visualization
  4. Proceed to Phase 3: Clinical Integration

===========================================
",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  length(completed_cancers),
  paste(completed_cancers, collapse = ", "),
  sum(summary_df$n_samples),
  mean(summary_df$n_genes_analyzed),
  unique(combined_histone$gene_name) %>% length(),
  sum(summary_df$consensus_sig),
  mean(summary_df$consensus_sig),
  paste(head(summary_df$cancer_type[order(-summary_df$consensus_sig)], 3), 
        collapse = ", ")
)

cat(summary_report)

writeLines(summary_report, "results/PAN_CANCER_AGGREGATION_SUMMARY.txt")

cat("\n✓ Aggregation complete!\n")
cat("✓ Report saved: results/PAN_CANCER_AGGREGATION_SUMMARY.txt\n")
cat("\n====================================================================\n")

