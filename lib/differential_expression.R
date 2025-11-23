
run_deseq2_analysis <- function(data, 
                                histone_genes,
                                padj_threshold = 0.05,
                                log2fc_threshold = 1.0,
                                correction_method = "BH") {
  
  library(DESeq2)
  library(SummarizedExperiment)
  library(dplyr)
  
  cat("\n=== DESeq2 Differential Expression Analysis ===\n")
  
  # Filter for histone genes
  gene_indices <- which(rowData(data)$gene_name %in% histone_genes)
  
  if (length(gene_indices) == 0) {
    cat("WARNING: No histone genes found in dataset\n")
    return(NULL)
  }
  
  cat("Genes detected:", length(gene_indices), "/", length(histone_genes), "\n")
  
  histone_data <- data[gene_indices, ]
  
  sample_types <- table(colData(histone_data)$sample_type)
  cat("Sample distribution:\n")
  print(sample_types)
  
  if (!all(c("Primary Tumor", "Solid Tissue Normal") %in% names(sample_types))) {
    cat("WARNING: Missing tumor or normal samples\n")
    return(NULL)
  }
  
  
  dds <- DESeqDataSet(histone_data, design = ~ sample_type)
  

  dds <- DESeq(dds)
  
  res <- results(
    dds,
    contrast = c("sample_type", "Primary Tumor", "Solid Tissue Normal"),
    alpha = padj_threshold,
    pAdjustMethod = correction_method
  )
  
  res_df <- as.data.frame(res)
  res_df$gene_name <- rowData(histone_data)$gene_name
  res_df$ensembl_id <- rownames(res_df)
  
  res_df <- res_df[, c("gene_name", "ensembl_id", "baseMean", 
                       "log2FoldChange", "lfcSE", "stat", 
                       "pvalue", "padj")]
  
  res_df$significant <- !is.na(res_df$padj) & 
    res_df$padj < padj_threshold &
    abs(res_df$log2FoldChange) > log2fc_threshold
  
  res_df$regulation <- ifelse(
    is.na(res_df$padj), "NS",
    ifelse(res_df$padj >= padj_threshold, "NS",
           ifelse(res_df$log2FoldChange > 0, "Upregulated", "Downregulated"))
  )
  
  res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
  
  if (sum(res_df$significant, na.rm = TRUE) > 0) {
    cat("  - Upregulated:", sum(res_df$regulation == "Upregulated"), "\n")
    cat("  - Downregulated:", sum(res_df$regulation == "Downregulated"), "\n")
    
    cat("\nTop significant genes:\n")
    print(res_df[res_df$significant, c("gene_name", "log2FoldChange", "padj")], 
          row.names = FALSE)
  } else {
    cat("No significant genes found\n")
  }
  
  return(list(
    results = res_df,
    dds = dds,
    normalized_counts = counts(dds, normalized = TRUE)
  ))
}

#------------------------- Edger Analysis ---------------
run_edger_analysis <- function(data,
                               histone_genes,
                               padj_threshold = 0.05,
                               log2fc_threshold = 1.0) {
  
  library(edgeR)
  library(SummarizedExperiment)
  
  gene_indices <- which(rowData(data)$gene_name %in% histone_genes)
  
  if (length(gene_indices) == 0) {
    return(NULL)
  }
  
  histone_data <- data[gene_indices, ]
  
  counts_matrix <- assay(histone_data, "unstranded")
  groups <- factor(colData(histone_data)$sample_type)
  
  dge <- DGEList(counts = counts_matrix, group = groups)
  
  keep <- filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  dge <- calcNormFactors(dge, method = "TMM")
  
  design <- model.matrix(~ groups)
  
  dge <- estimateDisp(dge, design)
  
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = 2)
  
  res <- topTags(qlf, n = Inf)$table
  res$gene_name <- rowData(histone_data[rownames(res), ])$gene_name
  
  res$significant <- res$FDR < padj_threshold & 
    abs(res$logFC) > log2fc_threshold
  
  cat("\nSignificant genes:", sum(res$significant), "\n")
  
  return(list(
    results = res,
    dge = dge
  ))
}


#-------- Differential Expression Analysis via Limma-Voom
run_limma_voom_analysis <- function(data,
                                    histone_genes,
                                    padj_threshold = 0.05,
                                    log2fc_threshold = 1.0) {
  
  library(limma)
  library(edgeR)
  library(SummarizedExperiment)

  gene_indices <- which(rowData(data)$gene_name %in% histone_genes)
  
  if (length(gene_indices) == 0) {
    return(NULL)
  }
  
  histone_data <- data[gene_indices, ]
  
  counts_matrix <- assay(histone_data, "unstranded")
  groups <- factor(colData(histone_data)$sample_type)
  
  dge <- DGEList(counts = counts_matrix)
  dge <- calcNormFactors(dge, method = "TMM")
  
  design <- model.matrix(~ groups)
  
  v <- voom(dge, design)
  
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  res <- topTable(fit, coef = 2, n = Inf)
  res$gene_name <- rowData(histone_data[rownames(res), ])$gene_name
  
  res$significant <- res$adj.P.Val < padj_threshold & 
    abs(res$logFC) > log2fc_threshold
  
  cat("\nSignificant genes:", sum(res$significant), "\n")
  
  return(list(
    results = res,
    voom_data = v
  ))
}


#--------------Comparasion of multiple DE methods
compare_de_methods <- function(data, histone_genes) {
  
  deseq2_res <- run_deseq2_analysis(data, histone_genes)
  edger_res <- run_edger_analysis(data, histone_genes)
  limma_res <- run_limma_voom_analysis(data, histone_genes)
  
  comparison <- data.frame(
    gene = histone_genes,
    deseq2_sig = NA,
    edger_sig = NA,
    limma_sig = NA,
    consensus = NA
  )
  
  for (gene in histone_genes) {
    if (!is.null(deseq2_res)) {
      comparison[comparison$gene == gene, "deseq2_sig"] <- 
        any(deseq2_res$results$gene_name == gene & 
              deseq2_res$results$significant, na.rm = TRUE)
    }
    
    if (!is.null(edger_res)) {
      comparison[comparison$gene == gene, "edger_sig"] <- 
        any(edger_res$results$gene_name == gene & 
              edger_res$results$significant, na.rm = TRUE)
    }
    
    if (!is.null(limma_res)) {
      comparison[comparison$gene == gene, "limma_sig"] <- 
        any(limma_res$results$gene_name == gene & 
              limma_res$results$significant, na.rm = TRUE)
    }
  }
  
  # Consensus: significant in at least 2 methods
  comparison$consensus <- rowSums(comparison[, c("deseq2_sig", "edger_sig", "limma_sig")], 
                                  na.rm = TRUE) >= 2
  
  cat("\n=== Method Comparison ===\n")
  print(comparison)
  
  return(list(
    deseq2 = deseq2_res,
    edger = edger_res,
    limma = limma_res,
    comparison = comparison
  ))
}

#------- Extraction of Normalized Counts
get_normalized_counts <- function(dds, genes = NULL) {
  
  library(DESeq2)
  
  norm_counts <- counts(dds, normalized = TRUE)
  
  if (!is.null(genes)) {
    gene_idx <- which(rowData(dds)$gene_name %in% genes)
    norm_counts <- norm_counts[gene_idx, ]
  }
  
  norm_counts_df <- as.data.frame(t(norm_counts))
  norm_counts_df$sample <- rownames(norm_counts_df)
  norm_counts_df$sample_type <- colData(dds)$sample_type
  
  return(norm_counts_df)
}


#--------------Batch effect correction using ComBat-seq
correct_batch_effects <- function(data, batch_variable = "plate") {
  
  library(sva)
  library(SummarizedExperiment)

  if (!batch_variable %in% colnames(colData(data))) {
    cat("WARNING: Batch variable", batch_variable, "not found\n")
    cat("Available variables:", colnames(colData(data)), "\n")
    return(data)
  }
  
  counts_matrix <- assay(data, "unstranded")
  batch <- colData(data)[[batch_variable]]
  
  complete_idx <- !is.na(batch)
  counts_matrix <- counts_matrix[, complete_idx]
  batch <- batch[complete_idx]
  
  cat("Samples before correction:", ncol(assay(data, "unstranded")), "\n")
  cat("Samples after removing NAs:", ncol(counts_matrix), "\n")
  cat("Number of batches:", length(unique(batch)), "\n")
  
  corrected_counts <- ComBat_seq(counts_matrix, batch = batch)
  
  assay(data, "unstranded") <- corrected_counts
  
  return(data)
}


#' Quality control plots for DE analysis
create_qc_plots <- function(dds, output_dir = "results/figures/qc") {
  
  library(DESeq2)
  library(ggplot2)
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("\nCreating quality control plots...\n")
  
  # 1. PCA plot
  vsd <- vst(dds, blind = FALSE)
  
  pdf(file.path(output_dir, "pca_plot.pdf"), width = 8, height = 6)
  plotPCA(vsd, intgroup = "sample_type") +
    theme_bw() +
    ggtitle("PCA: Tumor vs Normal Samples")
  dev.off()
  
  cat("PCA plot saved\n")
  
  # 2. Dispersion plot
  pdf(file.path(output_dir, "dispersion_plot.pdf"), width = 8, height = 6)
  plotDispEsts(dds)
  dev.off()
  
  cat("Dispersion plot saved\n")
  
  # 3. MA plot
  pdf(file.path(output_dir, "ma_plot.pdf"), width = 8, height = 6)
  plotMA(dds, ylim = c(-5, 5))
  dev.off()
  
  cat("MA plot saved\n")
  
  cat("Quality control plots saved in:", output_dir, "\n")
}