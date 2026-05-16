
# ============================================================================
# TCGA Histone Variant Analysis - Single Cancer Analysis
# Betül İrem Yardımcı
# ============================================================================
# Version: 2.1 - Single cancer run (memory efficient)
# Processing: One cancer at a time
# RAM Requirements: ~8-12 GB (with batch processing)
# ============================================================================

# Libraries
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(yaml)
  library(pheatmap)
})

# Source scripts
source("lib/tcga_data_extract.R")
source("lib/preprocessing.R")
source("lib/differential_expression_updated.R")
source("lib/correlation_analysis.R")
source("lib/expression_visualization.R")

# Load configuration
config <- yaml::yaml.load_file("config.yaml")

CURRENT_CANCER <- "LUAD"
cat("Current cancer:", CURRENT_CANCER, "\n")

data_category <- config$data_categories
histone_genes <- unlist(config$genes$tcga_gene_names)
ensembl_ids <- unlist(config$genes$ensembl_ids)
padj_threshold <- config$analysis$padj_threshold
log2fc_threshold <- config$analysis$log2fc_threshold
use_sva = config$analysis$use_sva
n_sv <- config$analysis$n_surrogate_variables

cat("Parameters:\n")
cat("  Histone genes:", length(histone_genes), "\n")
cat("  padj threshold:", padj_threshold, "\n")
cat("  log2FC threshold:", log2fc_threshold, "\n")
cat("  Use SVA:", use_sva, "\n\n")

cancer <- CURRENT_CANCER
project_name <- paste0("TCGA-", cancer)
category_name <- gsub(" ", "_", data_category$transcriptome)

results_dir <- paste0("results/", cancer)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------
# File paths
# --------------------------------------------------------------------
query_file <- paste0("data/queries/", project_name, "_",
                     category_name, "_query.RData")
expr_file <- paste0("data/processed/", cancer, "_expression.RData")
se_file <- paste0("data/processed/", cancer, "_se.RData")

# --------------------------------------------------------------------
# Step 1: Query + Download

cat("\n Step 1: Data Extraction\n")

if (file.exists(query_file)) {
  message("Loading existing query for ", cancer)
  load(query_file)
} else {
  message("Creating query and downloading data for ", cancer)
  query <- download_tcga_data(
    project_name = project_name,
    data_category = data_category$transcriptome
  )
  dir.create("data/queries", recursive = TRUE, showWarnings = FALSE)
  save(query, file = query_file)
}

# --------------------------------------------------------------------
# Step 2: Expression Matrix + SE (WITH BATCH PROCESSING!)

cat("\n Step 2: Expression Matrix Preparation\n")

if (file.exists(expr_file) && file.exists(se_file)) {
  message("Loading existing expression data for ", cancer)
  load(expr_file)   # loads: expr_data
  load(se_file)     # loads: se

  cat("Expression matrix:\n")
  cat("  Genes:", nrow(expr_data), "\n")
  cat("  Samples:", ncol(expr_data), "\n")

} else if (file.exists(expr_file) && !file.exists(se_file)) {

  message("Creating SE from existing expression matrix")
  load(expr_file)
  se <- convert_to_se(expr_data)
  save(se, file = se_file)

} else {

  message("Processing expression data from raw files")
  message("Using BATCH PROCESSING to manage memory...")

  # BATCH PROCESSING - Your existing function!
  expr_data <- process_in_batches_with_query(
    query = query,
    batch_size = 100  # Adjust based on your RAM (50-150)
  )

  save(expr_data, file = expr_file)

  # Filter samples
  message("Filtering samples...")
  filter_results <- sample_code_check(
    counts_matrix = expr_data,
    min_samples = 10,
    keep_codes = c("01", "11")  # Primary Tumor + Normal
  )

  expr_filtered <- filter_results$data_filtered

  # Convert to SummarizedExperiment
  se <- convert_to_se(expr_filtered)
  save(se, file = se_file)

  # Clean up
  rm(expr_data, expr_filtered)
  gc()
}

cat("\nSummarizedExperiment summary:\n")
cat("  Genes:", nrow(se), "\n")
cat("  Samples:", ncol(se), "\n")
cat("  Sample types:\n")
print(table(colData(se)$sample_type))

#-----------------------------------------
#ComBat-Seq Implementation

se_combat_file <- paste0("data/processed/", cancer, "_se_combat.RData")

if (file.exists(se_combat_file)) {
  cat("Loading ComBat-corrected SE...\n")
  load(se_combat_file)
} else {
  cat("Applying ComBat-seq...\n")
  se <- apply_combat_seq(se, min_samples = 3)
  save(se, file = se_combat_file)
}

#------------------------------------------
# Step 3: FULL TRANSCRIPTOME Differential Expression

gene_map <- setNames(histone_genes, ensembl_ids)

# DESeq2
deseq2_results <- run_deseq2_analysis(
  data = se,
  padj_threshold = padj_threshold,
  log2fc_threshold = log2fc_threshold,
  use_sva = FALSE,
  n_sv = NULL,
  condition = "sample_type",
  ref_group = "Normal",
  case_group = "Tumor"
)

deseq2_histone <- deseq2_results$results %>%
  filter(ensembl_id %in% ensembl_ids)

deseq2_histone$gene_name <- gene_map[deseq2_histone$ensembl_id]

#save full transcriptome
write.csv(deseq2_results$results,
          paste0(results_dir, "/deseq2_full_results.csv"),
          row.names = FALSE)
#save only histone genes
write.csv(deseq2_histone,
          paste0(results_dir, "/deseq2_histone.csv"),
          row.names = FALSE)


# Clean up to free memory
rm(se)
gc()

# edgeR
# Recreate se from saved file (memory efficient)
load(se_combat_file)
edger_results <- run_edger_analysis(
  data = se,
  padj_threshold = padj_threshold,
  log2fc_threshold = log2fc_threshold
)

edger_histone <- edger_results$results %>%
  filter(ensembl_id %in% ensembl_ids)

edger_histone$gene_name <- gene_map[edger_histone$ensembl_id]

#save full transcriptome
write.csv(edger_results$results,
          paste0(results_dir, "/edger_full_results.csv"),
          row.names = FALSE)
#save only histone genes
write.csv(edger_histone,
          paste0(results_dir, "/edger_histone.csv"),
          row.names = FALSE)

rm(se)
gc()

# limma-voom
load(se_combat_file)
limma_results <- run_limma_voom_analysis(
  data = se,
  padj_threshold = padj_threshold,
  log2fc_threshold = log2fc_threshold
)

limma_histone <- limma_results$results %>%
  filter(ensembl_id %in% ensembl_ids)

limma_histone$gene_name <- gene_map[limma_histone$ensembl_id]

#save full transcriptome
write.csv(limma_results$results,
          paste0(results_dir, "/limma_full_results.csv"),
          row.names = FALSE)
#save only histone genes
write.csv(limma_histone,
          paste0(results_dir, "/limma_histone.csv"),
          row.names = FALSE)

rm(se)
gc()

#control
cat("Histone genes found:\n")
cat("  DESeq2:", nrow(deseq2_histone), "\n")
cat("  edgeR:", nrow(edger_histone), "\n")
cat("  limma-voom:", nrow(limma_histone), "\n")

#--------------------------------------------
#Consensus Analysis
cat("\n Creating Method Consensus\n")

consensus <- compare_de_methods(
  deseq2_res = deseq2_results$results,
  edger_res = edger_results$results,
  limma_res = limma_results$results,
  histone_ensembl = ensembl_ids,
  tcga_gene_names = histone_genes,
  padj_threshold = padj_threshold,
  min_methods = 2
)

# save
write.csv(consensus,
          paste0(results_dir, "/consensus.csv"),
          row.names = FALSE)

# Normalized counts (for visualization)
norm_counts <- deseq2_results$normalized_counts
saveRDS(norm_counts,
        paste0(results_dir, "/normalized_counts.rds"))




