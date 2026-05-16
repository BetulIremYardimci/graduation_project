
# ============================================================================
# Differential Expression Analysis Functions
# ============================================================================
# Updated version with SVA support and enhanced features
# ============================================================================

run_deseq2_analysis <- function(data,
                                padj_threshold = 0.05,
                                log2fc_threshold = 1.0,
                                correction_method = "BH",
                                use_sva,
                                n_sv,
                                condition = "sample_type",
                                ref_group = "Normal",
                                case_group = "Tumor") {

  cat("\n=== DESeq2 Differential Expression Analysis ===\n")

  # Data summary
  cat("Genes:", nrow(data), "\n")
  cat("Samples:", ncol(data), "\n")
  cat("Sample distribution:\n")
  #print(table(colData(data)$sample_type))

  if (!is.factor(colData(data)$sample_type)) {
    colData(data)$sample_type <- factor(colData(data)$sample_type)
  }

  # Set reference level
  colData(data)$sample_type <- relevel(
    colData(data)$sample_type,
    ref = "Normal"
  )

  # Create DESeqDataSet with optional SVA
  if (use_sva) {
    cat("\n>>> Using SVA for batch effect correction...\n")
    library(sva)

    # Initial model without SVs
    dds_temp <- DESeqDataSet(data, design = ~ sample_type)
    dds_temp <- estimateSizeFactors(dds_temp)

    # Detect surrogate variables
    dat <- counts(dds_temp, normalized = TRUE)
    mod <- model.matrix(~ sample_type, colData(dds_temp))
    mod0 <- model.matrix(~ 1, colData(dds_temp))

    # Filter genes with very low counts
    idx <- rowMeans(dat) > 1
    cat("   Genes passing filter:", sum(idx), "/", nrow(dat), "\n")

    # Estimate number of SVs if not specified
    if (is.null(n_sv)) {
      n_sv <- num.sv(dat[idx,], mod, method = "leek")
      cat("   Estimated number of SVs:", n_sv, "\n")
    }

    # Run svaseq
    svseq <- svaseq(dat[idx,], mod, mod0, n.sv = n_sv)

    # Add SVs to colData
    for (i in 1:svseq$n.sv) {
      colData(data)[[paste0("SV", i)]] <- svseq$sv[,i]
    }

    # Create design formula with SVs
    sv_terms <- paste0("SV", 1:svseq$n.sv, collapse = " + ")
    design_formula <- as.formula(paste("~", sv_terms, "+ sample_type"))

    cat("   Design formula:", deparse(design_formula), "\n")

    # Create DESeqDataSet with SVs
    dds <- DESeqDataSet(data, design = design_formula)
    cat("✓ Added", svseq$n.sv, "surrogate variables\n")

  } else {
    dds <- DESeqDataSet(data, design = ~ sample_type)
  }

  #check values
  if (!is.integer(assay(data))) {
    cat(">>> Warning: Converting counts to integer mode...\n")
    assay(data) <- round(assay(data))
    mode(assay(data)) <- "integer"
  }

  # Run DESeq2
  cat("\n>>> Running DESeq2 pipeline...\n")
  dds <- DESeq(dds)

  # Get results (Tumor vs Normal)
  res <- results(
    dds,
    contrast = c(condition, case_group, ref_group),
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
  cat("Total genes analyzed:", nrow(res_df), "\n")
  cat("Upregulated genes:", sum(res_df$regulation == "Upregulated"), "\n")
  cat("Downregulated genes:", sum(res_df$regulation == "Downregulated"), "\n")
  cat("Not significant:", sum(res_df$regulation == "NS"), "\n")

  # Top significant genes
  if (sum(res_df$significant, na.rm = TRUE) > 0) {
    top_sig <- head(res_df[res_df$significant, ], 10)
    cat("\nTop 10 significant genes:\n")
    print(top_sig[, c("ensembl_id", "log2FoldChange", "padj")])
  } else {
    cat("\n! No significant genes found\n")
  }

  # Get normalized counts
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)

  # Return results
  return(list(
    results = res_df,
    dds = dds,
    normalized_counts = norm_counts,
    sva_used = use_sva
  ))
}


#------------------------- edgeR Analysis ---------------
run_edger_analysis <- function(data,
                               padj_threshold = 0.05,
                               log2fc_threshold = 1) {

  library(edgeR)

  message("\n=== edgeR Differential Expression Analysis ===")

  counts <- assay(data, "counts")
  group  <- colData(data)$sample_type

  # DGEList
  dge <- DGEList(counts = counts, group = group)

  # Filtering lowly expressed genes
  keep <- filterByExpr(dge)
  cat("Genes before filtering:", nrow(dge), "\n")
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  cat("Genes after filtering:", nrow(dge), "\n")

  # Normalization
  dge <- calcNormFactors(dge)

  # Design matrix
  design <- model.matrix(~ group)

  # Dispersion estimation
  dge <- estimateDisp(dge, design)

  # Fit model
  fit <- glmQLFit(dge, design)

  # Test Tumor vs Normal
  qlf <- glmQLFTest(fit, coef = 2)

  # Results (ALL genes)
  res <- topTags(qlf, n = Inf)$table

  # Add ensembl_id explicitly
  res$ensembl_id <- rownames(res)

  # Significance flag
  res$significant <-
    res$FDR < padj_threshold &
    abs(res$logFC) >= log2fc_threshold

  res$regulation <- ifelse(
    res$FDR >= padj_threshold, "NS",
    ifelse(res$logFC > 0, "Upregulated", "Downregulated")
  )

  # Reorder columns
  res <- res[, c("ensembl_id",
                 "logFC",
                 "logCPM",
                 "F",
                 "PValue",
                 "FDR",
                 "significant",
                 "regulation")]

  rownames(res) <- NULL

  # Summary
  cat("\n=== Results Summary ===\n")
  cat("Upregulated:", sum(res$regulation == "Upregulated"), "\n")
  cat("Downregulated:", sum(res$regulation == "Downregulated"), "\n")
  cat("Not significant:", sum(res$regulation == "NS"), "\n")

  return(list(results = res))
}


#-------- Limma-Voom Analysis
run_limma_voom_analysis <- function(data,
                                    padj_threshold = 0.05,
                                    log2fc_threshold = 1) {

  library(edgeR)
  library(limma)

  message("\n=== limma-voom Differential Expression Analysis ===")

  counts <- assay(data, "counts")
  group  <- colData(data)$sample_type

  # DGEList
  dge <- DGEList(counts = counts)

  # Filtering
  keep <- filterByExpr(dge, group = group)
  cat("Genes before filtering:", nrow(dge), "\n")
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  cat("Genes after filtering:", nrow(dge), "\n")

  # Normalization
  dge <- calcNormFactors(dge)

  # Design matrix
  design <- model.matrix(~ group)

  # voom transformation
  v <- voom(dge, design, plot = FALSE)

  # Linear model
  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  # Results (ALL genes)
  res <- topTable(fit,
                  coef = 2,
                  number = Inf,
                  sort.by = "none")

  # Add ensembl_id
  res$ensembl_id <- rownames(res)

  # Significance flag
  res$significant <-
    res$adj.P.Val < padj_threshold &
    abs(res$logFC) >= log2fc_threshold

  res$regulation <- ifelse(
    res$adj.P.Val >= padj_threshold, "NS",
    ifelse(res$logFC > 0, "Upregulated", "Downregulated")
  )

  # Reorder columns
  res <- res[, c("ensembl_id",
                 "logFC",
                 "AveExpr",
                 "t",
                 "P.Value",
                 "adj.P.Val",
                 "significant",
                 "regulation")]

  rownames(res) <- NULL

  # Summary
  cat("\n=== Results Summary ===\n")
  cat("Upregulated:", sum(res$regulation == "Upregulated"), "\n")
  cat("Downregulated:", sum(res$regulation == "Downregulated"), "\n")
  cat("Not significant:", sum(res$regulation == "NS"), "\n")

  return(list(results = res))
}


#--------------Comparison of DE Methods
compare_de_methods <- function(deseq2_res,
                                          edger_res,
                                          limma_res,
                                          histone_ensembl,
                                          tcga_gene_names,
                                          padj_threshold = 0.05,
                                          min_methods = 2) {

  cat("\n=== Method Comparison (Using 'significant' columns) ===\n")

  # Gene map
  gene_map <- setNames(tcga_gene_names, histone_ensembl)

  # Initialize
  comparison <- data.frame(
    gene_name = tcga_gene_names,
    ensembl_id = histone_ensembl,
    n_methods_sig = 0,
    consensus = FALSE,
    deseq2_sig = FALSE,
    edger_sig = FALSE,
    limma_sig = FALSE,
    deseq2_log2fc = NA,
    edger_logfc = NA,
    limma_logfc = NA,
    deseq2_padj = NA,
    edger_fdr = NA,
    limma_adjpval = NA,
    stringsAsFactors = FALSE
  )

  # For each histone gene
  for (i in 1:nrow(comparison)) {

    ens_id <- comparison$ensembl_id[i]
    gene_name <- comparison$gene_name[i]

    ## ----------------------------------------------------------------
    ## DESeq2 - USE 'significant' COLUMN!
    ## ----------------------------------------------------------------
    deseq2_match <- deseq2_res[deseq2_res$ensembl_id == ens_id, ]

    if (nrow(deseq2_match) > 0) {
      comparison$deseq2_sig[i] <- deseq2_match$significant  # ✅ USE THIS!
      comparison$deseq2_log2fc[i] <- deseq2_match$log2FoldChange
      comparison$deseq2_padj[i] <- deseq2_match$padj
    }

    ## ----------------------------------------------------------------
    ## edgeR - USE 'significant' COLUMN!
    ## ----------------------------------------------------------------
    edger_match <- edger_res[edger_res$ensembl_id == ens_id, ]

    if (nrow(edger_match) > 0) {
      comparison$edger_sig[i] <- edger_match$significant  # ✅ USE THIS!
      comparison$edger_logfc[i] <- edger_match$logFC
      comparison$edger_fdr[i] <- edger_match$FDR
    }

    ## ----------------------------------------------------------------
    ## limma - USE 'significant' COLUMN!
    ## ----------------------------------------------------------------
    limma_match <- limma_res[limma_res$ensembl_id == ens_id, ]

    if (nrow(limma_match) > 0) {
      comparison$limma_sig[i] <- limma_match$significant  # ✅ USE THIS!
      comparison$limma_logfc[i] <- limma_match$logFC
      comparison$limma_adjpval[i] <- limma_match$adj.P.Val
    }

    ## ----------------------------------------------------------------
    ## Calculate consensus
    ## ----------------------------------------------------------------
    comparison$n_methods_sig[i] <-
      sum(comparison$deseq2_sig[i],
          comparison$edger_sig[i],
          comparison$limma_sig[i],
          na.rm = TRUE)

    comparison$consensus[i] <-
      comparison$n_methods_sig[i] >= min_methods

    # Debug output
    cat(sprintf("%-12s | DESeq2: %5s | edgeR: %5s | limma: %5s | Consensus: %5s (%d/3)\n",
                gene_name,
                comparison$deseq2_sig[i],
                comparison$edger_sig[i],
                comparison$limma_sig[i],
                comparison$consensus[i],
                comparison$n_methods_sig[i]))
  }

  cat("\n=== Summary ===\n")
  cat("Consensus significant (>=", min_methods, " methods):",
      sum(comparison$consensus), "/", nrow(comparison), "\n\n")

  return(comparison)
}

#------- Quality Control Plots
create_qc_plots <- function(dds, output_dir = "results/figures/qc") {

  library(DESeq2)
  library(ggplot2)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n>>> Creating quality control plots...\n")

  # Transformation for visualization
  tryCatch({
    if (nrow(dds) < 1000) {
      cat("   Using rlog transformation (small gene set)...\n")
      rld <- rlog(dds, blind = FALSE)
      plot_data <- rld
    } else {
      cat("   Using VST transformation...\n")
      vsd <- vst(dds, blind = FALSE)
      plot_data <- vsd
    }

    # 1. PCA plot
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

  # 4. Sample distance heatmap
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

  cat("\n✓ Quality control plots completed\n")
  cat("   Location:", output_dir, "\n")
}
