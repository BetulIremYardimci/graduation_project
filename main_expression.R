#updated: 20 May 2026
# Betül İrem YARDIMCI

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

#lib/
source("lib/tcga_data_extract.R")
source("lib/preprocessing.R")
source("lib/differential_expression_updated.R")
source("lib/nonparametric_tests.R")


config <- yaml::yaml.load_file("config.yaml")

data_category <- config$data_categories
histone_genes <- unlist(config$genes$tcga_gene_names)
ensembl_ids <- unlist(config$genes$ensembl_ids)
padj_threshold <- config$analysis$padj_threshold
log2fc_threshold <- config$analysis$log2fc_threshold
use_sva = config$analysis$use_sva
n_sv <- config$analysis$n_surrogate_variables

cancer <- "OV"
project_name <- paste0("TCGA-", cancer)
category_name <- data_category$transcriptome
results_dir <- paste0("results/", cancer)

#directory check for cancer type
if (dir.exists(results_dir)){
  message("Results directory exists for ", cancer)
  results_dir <- paste0("results/", cancer)
} else{
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
}

# data files check
query_file <- paste0("data/queries/", project_name, "_",
                     category_name, "_query.RData")
expr_file <- paste0("data/processed/", cancer, "_expression.RData")
se_file <- paste0("data/processed/", cancer, "_se.RData")

#----DATA PROCESSING

#1. Data Query
if (file.exists(query_file)) {
  message("Loading existing query for ", cancer)
  load(query_file)
} else {
  message("Creating query and downloading data for ", cancer)
  query <- download_tcga_data(
    project_name = project_name,
    data_category = category_name
  )
  dir.create("data/queries", recursive = TRUE, showWarnings = FALSE)
  save(query, file = query_file)
}

#2. Create Expression Matrix (SE) and Batch Processing
if (file.exists(expr_file) && file.exists(se_file)){
  message("Existing expression data is loading for ", cancer)
  load(expr_file)    # expr_data objesini yükler
  load(se_file)      # se objesini yükler

} else if (file.exists(expr_file) && !file.exists(se_file)) {
  message("Only Expression file exists for ", cancer)
  load(expr_file)    # expr_data yükler
  se <- convert_to_se(expr_data)
  save(se, file = se_file)
} else {
  message("Processing expression data from raw files by implementing Batch
          Processing")
  #batch processing:
  expr_data <- process_in_batches_with_query(
    query = query,
    batch_size = 100  #can be adjusted based on your RAM
  )
  save(expr_data, file = expr_file) #loaded file will be imported as expr_data

  # Filter samples
  message("Filtering samples...")
  filter_results <- sample_code_check(
    counts_matrix = expr_data,
    min_samples = 10,
    keep_codes = c("01", "11")  # Primary Tumor + Normal
  )

  expr_filtered <- filter_results$data_filtered
  se <- convert_to_se(expr_filtered)
  save(se, file = se_file)

  rm(expr_data, expr_filtered)
  gc()
}

# 3. Combat-Seq Implementation
se_combat_file <- paste0("data/processed/", cancer, "_se_combat.RData")

if (file.exists(se_combat_file)) {
  cat("Loading ComBat-corrected SE...\n")
  load(se_combat_file)   # se objesini yükler
} else {
  cat("Applying ComBat-seq...\n")
  se <- apply_combat_seq(se, min_samples = 3)
  save(se, file = se_combat_file)
}

#---- EXPRESSION ANALYSIS
gene_map <- setNames(histone_genes, ensembl_ids) #assing Ensembl ID to histone genes

#DeSeq2
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

write.csv(deseq2_results$results,
          paste0(results_dir, "/deseq2_full_results.csv"),
          row.names = FALSE)    #Full transcriptome

write.csv(deseq2_histone,
          paste0(results_dir, "/deseq2_histone.csv"),
          row.names = FALSE). #Only histone transcriptome

# Clean up to free memory
rm(se)
gc()

#Limma-voom
load(se_combat_file)
limma_results <- run_limma_voom_analysis(
  data = se,
  padj_threshold = padj_threshold,
  log2fc_threshold = log2fc_threshold
)

limma_histone <- limma_results$results %>%
  filter(ensembl_id %in% ensembl_ids)

limma_histone$gene_name <- gene_map[limma_histone$ensembl_id]

write.csv(limma_results$results,
          paste0(results_dir, "/limma_full_results.csv"),
          row.names = FALSE)
write.csv(limma_histone,
          paste0(results_dir, "/limma_histone.csv"),
          row.names = FALSE)

#edgeR
load(se_combat_file)
edger_results <- run_edger_analysis(
  data = se,
  padj_threshold = padj_threshold,
  log2fc_threshold = log2fc_threshold
)

edger_histone <- edger_results$results %>%
  filter(ensembl_id %in% ensembl_ids)

edger_histone$gene_name <- gene_map[edger_histone$ensembl_id]

write.csv(edger_results$results,
          paste0(results_dir, "/edger_full_results.csv"),
          row.names = FALSE)

write.csv(edger_histone,
          paste0(results_dir, "/edger_histone.csv"),
          row.names = FALSE)

#Control
cat("Histone genes found:\n")
cat("  DESeq2:", nrow(deseq2_histone), "\n")
cat("  edgeR:", nrow(edger_histone), "\n")
cat("  limma-voom:", nrow(limma_histone), "\n")

#Non-parametric Tests Application
vst_counts     <- readRDS(paste0("results/", cancer, "/vst_counts.rds"))
deseq2_histone <- read.csv(paste0("results/", cancer, "/deseq2_histone.csv"))

wilcox_results <- run_nonparametric_tests(
  expr_matrix = vst_counts,
  histone_df  = deseq2_histone,
  cancer      = cancer,
  padj_method = "BH"
)

#optional
norm_counts <- readRDS(paste0("results/", cancer, "/normalized_counts.rds"))
wilcox_norm <- run_nonparametric_tests(
  expr_matrix = norm_counts,
  histone_df  = deseq2_histone,
  cancer      = cancer,
  padj_method = "BH"
)


#----- Data Visualization
plot_dir <- paste0("results/", cancer, "/plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

p_enhanced <- plot_histone_boxplot_vst_enhanced(
  cancer        = cancer,
  vst_counts    = vst_counts,
  deseq2_histone = deseq2_histone,
  wilcox_results = wilcox_results
)

ggsave(paste0(plot_dir, "/histone_boxplot_deseq2_wilcoxon.jpeg"),
       plot  = p_enhanced,
       width = 18, height = 8, dpi = 300,
       device = "jpeg", bg = "white")

