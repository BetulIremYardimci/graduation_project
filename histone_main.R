
# Load required libraries
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(yaml)
})

source("lib/tcga_dataset_extraction.R")
source("lib/differential_expression.R")
source("lib/survival_analysis.R")
source("lib/visualization.R")

config <- yaml::yaml.load_file("config.yaml")

cancer_types <- config$cancer_types
cat("Cancer types to analyze:", length(cancer_types), "\n")

histone_genes <- unlist(config$genes$tcga_gene_names)
cat("Genes to analyze:", paste(histone_genes, collapse = ", "), "\n")

params <- config$analysis
cat("Analysis mode:", ifelse(params$use_full_dataset, "FULL DATASET", "TEST MODE"), "\n")
if (!params$use_full_dataset) {
  cat("Test sample size:", params$test_sample_size, "\n")
}
cat("Statistical thresholds: padj <", params$padj_threshold, 
    ", |log2FC| >", params$log2fc_threshold, "\n")

output_dirs <- config$output
dir.create(output_dirs$base_dir, showWarnings = FALSE)
dir.create(output_dirs$differential_expression, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dirs$survival, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dirs$figures, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dirs$tables, showWarnings = FALSE, recursive = TRUE)

cat("  -", output_dirs$differential_expression, "\n")
cat("  -", output_dirs$survival, "\n")
cat("  -", output_dirs$figures, "\n")
cat("  -", output_dirs$tables, "\n")

cat("Cancer types:", length(cancer_types), "\n")
cat("Genes analyzed:", paste(histone_genes, collapse = ", "), "\n")
cat("Statistical method:", params$differential_expression_method, "\n")
cat("Multiple testing correction:", params$multiple_testing_correction, "\n")


all_results <- list()
pan_cancer_summary <- data.frame()

for (cancer in cancer_types) {
  
  cat("\n")
  cat(" ANALYZING:", cancer, "\n")

  tryCatch({
    
    cat("\n[1/6] Downloading TCGA-", cancer, " data...\n", sep = "")
    
    sample_size <- if (params$use_full_dataset) NULL else params$test_sample_size
    
    query <- download_tcga_data(
      cancer_type = cancer,
      sample_size = sample_size
    )
    
    cat("[1/6] Preparing data...\n")
    data <- prepare_tcga_data(query)
    
    cat("      Samples:", ncol(data), "\n")
    cat("      Genes:", nrow(data), "\n")

    cat("\n[2/6] Quality control and filtering...\n")
    
    gene_indices <- which(rowData(data)$gene_name %in% histone_genes)
    
    if (length(gene_indices) == 0) {
      cat("      WARNING: No histone genes found in", cancer, "\n")
      next
    }
    
    cat("      Found", length(gene_indices), "histone genes\n")
    
    histone_data <- data[gene_indices, ]
    
    sample_types <- table(colData(histone_data)$sample_type)
    cat("      Sample types:\n")
    print(sample_types)
    
    if (!("Primary Tumor" %in% names(sample_types) && 
          "Solid Tissue Normal" %in% names(sample_types))) {
      cat("      WARNING: Missing tumor or normal samples in", cancer, "\n")
      next
    }
    
    cat("\n[3/6] Differential expression analysis (DESeq2)...\n")
    
    de_analysis <- run_deseq2_analysis(
      data = histone_data,
      histone_genes = histone_genes,
      padj_threshold = params$padj_threshold,
      log2fc_threshold = params$log2fc_threshold,
      correction_method = params$multiple_testing_correction
    )
    
    if (is.null(de_analysis)) {
      cat("      WARNING: Differential expression analysis failed\n")
      next
    }
    
    de_results_df <- de_analysis$results
    dds <- de_analysis$dds
    normalized_counts <- de_analysis$normalized_counts
    
    # Save results
    if (params$save_intermediate_results) {
      write.csv(
        de_results_df,
        file = file.path(output_dirs$differential_expression, 
                         paste0(cancer, "_DE_results.csv")),
        row.names = FALSE
      )
    }

    cat("\n[4/6] Survival analysis...\n")
    
    survival_analysis <- run_complete_survival_analysis(
      data = histone_data,
      normalized_counts = normalized_counts,
      genes = histone_genes,
      cancer_type = cancer,
      output_dir = output_dirs$survival
    )
    
    survival_results <- if (!is.null(survival_analysis)) {
      survival_analysis$univariate
    } else {
      NULL
    }

    cat("\n[5/6] Subgroup analysis...\n")
    
    subgroup_results <- NULL
    
    if (cancer == "BRCA") {
      # ER status stratification for breast cancer
      if ("er_status_by_ihc" %in% colnames(clinical)) {
        cat("      Analyzing ER+ vs ER- subgroups...\n")
        # Implement ER-stratified analysis here
        # This is where web tools fail - cannot do subtype-specific analysis
      }
    } else if (cancer == "LUAD") {
      # Stage stratification for lung cancer
      if ("tumor_stage" %in% colnames(clinical)) {
        cat("      Analyzing early vs advanced stage...\n")
        # Implement stage-stratified analysis
      }
    }
    
    if (params$create_plots) {
      cat("\n[6/6] Creating visualizations...\n")
      
      tryCatch({
        create_expression_boxplot(
          data = histone_data,
          genes = histone_genes,
          cancer_type = cancer,
          output_file = file.path(output_dirs$figures, 
                                  paste0(cancer, "_expression_boxplot.pdf"))
        )
        cat("      Saved expression boxplot\n")
      }, error = function(e) {
        cat("      Warning: Could not create boxplot -", e$message, "\n")
      })
      
      # Create volcano plot
      tryCatch({
        create_volcano_plot(
          de_results = de_results_df,
          padj_threshold = params$padj_threshold,
          log2fc_threshold = params$log2fc_threshold,
          cancer_type = cancer,
          output_file = file.path(output_dirs$figures, 
                                  paste0(cancer, "_volcano_plot.pdf"))
        )
        cat("      Saved volcano plot\n")
      }, error = function(e) {
        cat("      Warning: Could not create volcano plot -", e$message, "\n")
      })
    }
    
    all_results[[cancer]] <- list(
      differential_expression = de_results_df,
      survival = survival_results,
      subgroup = subgroup_results,
      n_samples = ncol(data),
      n_significant = sum(de_results_df$significant, na.rm = TRUE)
    )
    
    pan_cancer_summary <- rbind(
      pan_cancer_summary,
      data.frame(
        cancer_type = cancer,
        n_samples = ncol(data),
        n_tumor = sum(colData(data)$sample_type == "Primary Tumor"),
        n_normal = sum(colData(data)$sample_type == "Solid Tissue Normal"),
        n_genes_detected = length(gene_indices),
        n_significant_de = sum(de_results_df$significant, na.rm = TRUE),
        n_significant_survival = if (!is.null(survival_results)) {
          sum(survival_results$pvalue < 0.05, na.rm = TRUE)
        } else {
          0
        }
      )
    )
    
    cat("\n[COMPLETE]", cancer, "analysis finished\n")
    
  }, error = function(e) {
    cat("\n[ERROR] Failed to analyze", cancer, ":", e$message, "\n")
  })
  
  Sys.sleep(2)
}


cat("\n")
cat("  PAN-CANCER SUMMARY\n")

if (params$save_intermediate_results) {
  write.csv(
    pan_cancer_summary,
    file = file.path(output_dirs$tables, "pan_cancer_summary.csv"),
    row.names = FALSE
  )
  cat("Saved pan-cancer summary table\n")
}

cat("\nAnalysis Summary:\n")
print(pan_cancer_summary)

# Multi-cancer comparison in single figure

create_pancancer_heatmap(
  all_results = all_results,
  genes = histone_genes,
  cancer_types = names(all_results),
  output_file = file.path(output_dirs$figures, "pan_cancer_expression_heatmap.pdf")
)


cat("\nCreating pan-cancer survival forest plots...\n")

survival_compiled <- do.call(rbind, lapply(names(all_results), function(cancer) {
  all_results[[cancer]]$survival
}))

if (!is.null(survival_compiled) && nrow(survival_compiled) > 0) {
  
  for (gene in unique(survival_compiled$gene)) {
    create_forest_plot(
      survival_results = survival_compiled,
      gene = gene,
      output_file = file.path(output_dirs$figures, 
                              paste0("pan_cancer_", gene, "_forest_plot.pdf"))
    )
  }
  
  cat("Forest plots created for all genes\n")
}



cat("  ANALYSIS COMPLETE\n")
cat("\nTotal cancers analyzed:", nrow(pan_cancer_summary), "\n")
cat("Total samples:", sum(pan_cancer_summary$n_samples), "\n")
cat("Cancers with significant DE:", 
    sum(pan_cancer_summary$n_significant_de > 0), "\n")
cat("Cancers with prognostic genes:", 
    sum(pan_cancer_summary$n_significant_survival > 0), "\n")
cat("\nResults saved in:", output_dirs$base_dir, "\n")
cat("  - Differential expression:", output_dirs$differential_expression, "\n")
cat("  - Survival analysis:", output_dirs$survival, "\n")
cat("  - Figures:", output_dirs$figures, "\n")
cat("  - Tables:", output_dirs$tables, "\n")
cat("\n")

saveRDS(all_results, file = file.path(output_dirs$base_dir, "all_results.rds"))
cat("Complete results saved:", file.path(output_dirs$base_dir, "all_results.rds"), "\n")
