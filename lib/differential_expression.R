
run_deseq2_analysis <- function(data,
                                padj_threshold = 0.05,
                                log2fc_threshold = 1.0,
                                correction_method = "BH") {

  library(DESeq2)
  library(SummarizedExperiment)

  cat("\n=== DESeq2 Differential Expression Analysis ===\n")

  # Data summary
  cat("Genes:", nrow(data), "\n")
  cat("Samples:", ncol(data), "\n")
  cat("Sample distribution:\n")
  print(table(colData(data)$sample_type))

  # Create DESeqDataSet
  dds <- DESeqDataSet(data, design = ~ sample_type)

  # Run DESeq2
  dds <- DESeq(dds)

  # Get results (Tumor vs Normal)
  res <- results(
    dds,
    contrast = c("sample_type", "Tumor", "Normal"),  # İsimler değişti
    alpha = padj_threshold,
    pAdjustMethod = correction_method
  )

  # Convert to dataframe
  res_df <- as.data.frame(res)
  res_df$ensembl_id <- rownames(res_df)

  res_df <- res_df[, c("ensembl_id",
                       "baseMean",
                       "log2FoldChange",
                       "lfcSE",
                       "stat",
                       "pvalue",
                       "padj")]


  # Add significance flags
  res_df$significant <- !is.na(res_df$padj) &
    res_df$padj < padj_threshold &
    abs(res_df$log2FoldChange) > log2fc_threshold

  res_df$regulation <- ifelse(
    is.na(res_df$padj), "NS",
    ifelse(res_df$padj >= padj_threshold, "NS",
           ifelse(res_df$log2FoldChange > 0, "Upregulated", "Downregulated"))
  )

  # Sort by p-value
  res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

  # Print results
  cat("\n=== Results Summary ===\n")
  cat("Upregulated genes:", sum(res_df$regulation == "Upregulated"), "\n")
  cat("Downregulated genes:", sum(res_df$regulation == "Downregulated"), "\n")
  cat("Not significant:", sum(res_df$regulation == "NS"), "\n")

  if (sum(res_df$significant, na.rm = TRUE) > 0) {
    cat("\nSignificant genes:\n")
    print(res_df[res_df$significant, c("ensembl_id", "log2FoldChange", "padj")],
          row.names = FALSE)
  } else {
    cat("\nNo significant genes found\n")
  }

  # Get normalized counts (düzeltme!)
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)

  # Return results
  return(list(
    results = res_df,
    dds = dds,
    normalized_counts = norm_counts
  ))
}

#------------------------- Edger Analysis ---------------
run_edger_analysis <- function(data,
                               padj_threshold = 0.05,
                               log2fc_threshold = 1) {

  library(edgeR)

  message("\n=== edgeR Differential Expression Analysis (full transcriptome) ===")

  counts <- assay(data, "counts")
  group  <- colData(data)$sample_type

  ## DGEList
  dge <- DGEList(counts = counts, group = group)

  ## Filtering lowly expressed genes (standard)
  keep <- filterByExpr(dge)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]

  ## Normalization
  dge <- calcNormFactors(dge)

  ## Design matrix
  design <- model.matrix(~ group)

  ## Dispersion estimation
  dge <- estimateDisp(dge, design)

  ## Fit model
  fit <- glmQLFit(dge, design)

  ## Test Tumor vs Normal
  qlf <- glmQLFTest(fit, coef = 2)

  ## Results (ALL genes)
  res <- topTags(qlf, n = Inf)$table

  ## Add ensembl_id explicitly
  res$ensembl_id <- rownames(res)

  ## Significance flag (optional, convenient)
  res$significant <-
    res$FDR < padj_threshold &
    abs(res$logFC) >= log2fc_threshold

  ## Reorder columns
  res <- res[, c("ensembl_id",
                 "logFC",
                 "logCPM",
                 "F",
                 "PValue",
                 "FDR",
                 "significant")]

  rownames(res) <- NULL

  return(list(results = res))
}



#-------- Differential Expression Analysis via Limma-Voom
run_limma_voom_analysis <- function(data,
                                    padj_threshold = 0.05,
                                    log2fc_threshold = 1) {

  library(edgeR)
  library(limma)

  message("\n=== limma-voom Differential Expression Analysis (full transcriptome) ===")

  counts <- assay(data, "counts")
  group  <- colData(data)$sample_type

  ## DGEList
  dge <- DGEList(counts = counts)

  ## Filtering low expression
  keep <- filterByExpr(dge, group = group)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]

  ## Normalization
  dge <- calcNormFactors(dge)

  ## Design matrix
  design <- model.matrix(~ group)

  ## voom transformation
  v <- voom(dge, design, plot = FALSE)

  ## Linear model
  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  ## Results (ALL genes)
  res <- topTable(fit,
                  coef = 2,
                  number = Inf,
                  sort.by = "none")

  ## Add ensembl_id explicitly
  res$ensembl_id <- rownames(res)

  ## Significance flag
  res$significant <-
    res$adj.P.Val < padj_threshold &
    abs(res$logFC) >= log2fc_threshold

  ## Reorder columns
  res <- res[, c("ensembl_id",
                 "logFC",
                 "AveExpr",
                 "t",
                 "P.Value",
                 "adj.P.Val",
                 "significant")]

  rownames(res) <- NULL

  return(list(results = res))
}


#--------------Comparasion of multiple DE methods
compare_de_methods <- function(deseq2_res,
                               edger_res,
                               limma_res,
                               histone_ensembl,
                               tcga_gene_names,
                               padj_threshold = 0.05,
                               min_methods = 2) {

  ## güvenlik: isimler eşleşiyor mu
  stopifnot(
    length(histone_ensembl) == length(tcga_gene_names),
    !is.null(names(tcga_gene_names)) ||
      is.null(names(tcga_gene_names))
  )

  ## ensembl -> gene_name map
  gene_map <- setNames(tcga_gene_names, histone_ensembl)

  ## initialize comparison table (gene_name görünür olacak)
  comparison <- data.frame(
    gene_name = tcga_gene_names,
    ensembl_id = histone_ensembl,
    deseq2_sig = FALSE,
    edger_sig  = FALSE,
    limma_sig  = FALSE,
    stringsAsFactors = FALSE
  )

  ## DESeq2
  sig_deseq2 <- deseq2_res$ensembl_id[
    deseq2_res$padj < padj_threshold
  ]

  comparison$deseq2_sig <-
    comparison$ensembl_id %in% sig_deseq2

  ## edgeR
  sig_edger <- edger_res$ensembl_id[
    edger_res$FDR < padj_threshold
  ]

  comparison$edger_sig <-
    comparison$ensembl_id %in% sig_edger

  ## limma
  sig_limma <- limma_res$ensembl_id[
    limma_res$adj.P.Val < padj_threshold
  ]

  comparison$limma_sig <-
    comparison$ensembl_id %in% sig_limma

  ## consensus (>= min_methods)
  comparison$consensus <-
    rowSums(comparison[, c("deseq2_sig",
                           "edger_sig",
                           "limma_sig")]) >= min_methods

  ## final table: gene_name önde
  comparison <- comparison[, c("gene_name",
                               "ensembl_id",
                               "deseq2_sig",
                               "edger_sig",
                               "limma_sig",
                               "consensus")]

  return(comparison)
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

  # 1. PCA plot - Use rlog instead of vst for small gene sets
  tryCatch({
    if (nrow(dds) < 1000) {
      # For small gene sets, use rlog
      cat("Using rlog transformation (small gene set)...\n")
      rld <- rlog(dds, blind = FALSE)
      plot_data <- rld
    } else {
      # For large gene sets, use vst
      cat("Using VST transformation...\n")
      vsd <- vst(dds, blind = FALSE)
      plot_data <- vsd
    }

    # PCA plot
    pdf(file.path(output_dir, "pca_plot.pdf"), width = 8, height = 6)
    p <- plotPCA(plot_data, intgroup = "sample_type") +
      theme_bw() +
      ggtitle("PCA: Tumor vs Normal Samples") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    print(p)
    dev.off()

    cat("✓ PCA plot saved\n")

  }, error = function(e) {
    cat("Warning: Could not create PCA plot:", e$message, "\n")
  })

  # 2. Dispersion plot
  tryCatch({
    pdf(file.path(output_dir, "dispersion_plot.pdf"), width = 8, height = 6)
    plotDispEsts(dds, main = "Dispersion Estimates")
    dev.off()
    cat("✓ Dispersion plot saved\n")
  }, error = function(e) {
    cat("Warning: Could not create dispersion plot:", e$message, "\n")
  })

  # 3. MA plot
  tryCatch({
    pdf(file.path(output_dir, "ma_plot.pdf"), width = 8, height = 6)
    plotMA(dds, ylim = c(-5, 5), main = "MA Plot")
    dev.off()
    cat("✓ MA plot saved\n")
  }, error = function(e) {
    cat("Warning: Could not create MA plot:", e$message, "\n")
  })

  # 4. Sample distance heatmap (works better for small gene sets)
  tryCatch({
    if (nrow(dds) < 1000) {
      rld <- rlog(dds, blind = FALSE)
    } else {
      rld <- vst(dds, blind = FALSE)
    }

    sampleDists <- dist(t(assay(rld)))
    sampleDistMatrix <- as.matrix(sampleDists)

    library(pheatmap)
    pdf(file.path(output_dir, "sample_distance_heatmap.pdf"),
        width = 10, height = 10)
    pheatmap(sampleDistMatrix,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists,
             main = "Sample Distance Heatmap",
             show_rownames = FALSE,
             show_colnames = FALSE)
    dev.off()
    cat("✓ Sample distance heatmap saved\n")

  }, error = function(e) {
    cat("Warning: Could not create distance heatmap:", e$message, "\n")
  })

  cat("\nQuality control plots saved in:", output_dir, "\n")
}
