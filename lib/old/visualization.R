
#expression boxplots (tumor vs normal)

create_expression_boxplot <- function(data, 
                                      genes, 
                                      cancer_type = "",
                                      output_file = NULL) {
  
  library(ggplot2)
  library(ggpubr)
  library(SummarizedExperiment)
  library(dplyr)
  
  cat("\n=== Creating Expression Boxplot ===\n")
  
  # Filter for histone genes
  gene_indices <- which(rowData(data)$gene_name %in% genes)
  
  if (length(gene_indices) == 0) {
    cat("WARNING: No genes found\n")
    return(NULL)
  }
  
  histone_data <- data[gene_indices, ]
  
  # Prepare plot data
  plot_data <- data.frame(
    expression = as.vector(t(assay(histone_data, "unstranded"))),
    gene = rep(rowData(histone_data)$gene_name, each = ncol(histone_data)),
    sample_type = rep(colData(histone_data)$sample_type, 
                      times = nrow(histone_data))
  )
  
  # Calculate statistics for each gene
  stats <- plot_data %>%
    group_by(gene, sample_type) %>%
    summarize(
      median = median(expression),
      n = n(),
      .groups = "drop"
    )
  
  cat("Sample sizes:\n")
  print(stats)
  
  # Create plot
  p <- ggplot(plot_data, aes(x = sample_type, y = log2(expression + 1), 
                             fill = sample_type)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
    facet_wrap(~gene, scales = "free_y", ncol = 2) +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.signif",
                       label.y.npc = 0.95) +
    scale_fill_manual(values = c("Primary Tumor" = "#E41A1C", 
                                 "Solid Tissue Normal" = "#377EB8")) +
    theme_bw(base_size = 12) +
    labs(
      title = paste("Histone Variant Expression in", cancer_type),
      x = "Sample Type",
      y = "log2(Count + 1)",
      fill = "Sample Type"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold")
    )
  
  # Save plot
  if (!is.null(output_file)) {
    ggsave(
      filename = output_file,
      plot = p,
      width = 8,
      height = 6,
      dpi = 300
    )
    cat("Saved boxplot:", output_file, "\n")
  }
  
  return(p)
}


#' Create volcano plot for differential expression
#' 
#' @param de_results Data frame from DESeq2 analysis
#' @param padj_threshold Adjusted p-value threshold
#' @param log2fc_threshold Log2 fold change threshold
#' @param cancer_type Cancer type label
#' @param output_file Path to save plot
create_volcano_plot <- function(de_results,
                                padj_threshold = 0.05,
                                log2fc_threshold = 1.0,
                                cancer_type = "",
                                output_file = NULL) {
  
  library(ggplot2)
  library(ggrepel)
  
  cat("\n=== Creating Volcano Plot ===\n")
  
  # Prepare data
  volcano_data <- de_results
  volcano_data$neg_log10_padj <- -log10(volcano_data$padj)
  
  # Color by significance
  volcano_data$significance <- "NS"
  volcano_data$significance[
    !is.na(volcano_data$padj) & 
      volcano_data$padj < padj_threshold & 
      volcano_data$log2FoldChange > log2fc_threshold
  ] <- "Upregulated"
  
  volcano_data$significance[
    !is.na(volcano_data$padj) & 
      volcano_data$padj < padj_threshold & 
      volcano_data$log2FoldChange < -log2fc_threshold
  ] <- "Downregulated"
  
  # Label significant genes
  volcano_data$label <- ifelse(
    volcano_data$significance != "NS",
    volcano_data$gene_name,
    ""
  )
  
  cat("Upregulated:", sum(volcano_data$significance == "Upregulated"), "\n")
  cat("Downregulated:", sum(volcano_data$significance == "Downregulated"), "\n")
  
  # Create plot
  p <- ggplot(volcano_data, aes(x = log2FoldChange, y = neg_log10_padj)) +
    geom_point(aes(color = significance), alpha = 0.6, size = 2) +
    geom_text_repel(aes(label = label), 
                    size = 4, 
                    max.overlaps = 20,
                    box.padding = 0.5) +
    geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), 
               linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(padj_threshold), 
               linetype = "dashed", color = "grey50") +
    scale_color_manual(
      values = c("Upregulated" = "#E41A1C", 
                 "Downregulated" = "#377EB8", 
                 "NS" = "grey70")
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = paste("Differential Expression Volcano Plot:", cancer_type),
      subtitle = paste0("Tumor vs Normal (padj < ", padj_threshold, 
                        ", |log2FC| > ", log2fc_threshold, ")"),
      x = "log2 Fold Change",
      y = "-log10(adjusted p-value)",
      color = "Significance"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  
  # Save plot
  if (!is.null(output_file)) {
    ggsave(
      filename = output_file,
      plot = p,
      width = 8,
      height = 6,
      dpi = 300
    )
    cat("Saved volcano plot:", output_file, "\n")
  }
  
  return(p)
}


#' Create pan-cancer expression heatmap
#' 
#' @param all_results List of differential expression results per cancer
#' @param genes Vector of gene names
#' @param cancer_types Vector of cancer type codes
#' @param output_file Path to save plot
create_pancancer_heatmap <- function(all_results,
                                     genes,
                                     cancer_types,
                                     output_file = NULL) {
  
  library(pheatmap)
  library(RColorBrewer)
  
  cat("\n=== Creating Pan-Cancer Heatmap ===\n")
  
  # Create matrix for log2 fold changes
  fc_matrix <- matrix(NA, 
                      nrow = length(genes), 
                      ncol = length(cancer_types))
  rownames(fc_matrix) <- genes
  colnames(fc_matrix) <- cancer_types
  
  # Fill matrix
  for (cancer in cancer_types) {
    if (cancer %in% names(all_results) && 
        !is.null(all_results[[cancer]]$differential_expression)) {
      
      de_data <- all_results[[cancer]]$differential_expression
      
      for (gene in genes) {
        gene_data <- de_data[de_data$gene_name == gene, ]
        if (nrow(gene_data) > 0 && !is.na(gene_data$log2FoldChange[1])) {
          fc_matrix[gene, cancer] <- gene_data$log2FoldChange[1]
        }
      }
    }
  }
  
  cat("Matrix dimensions:", dim(fc_matrix), "\n")
  cat("Non-NA values:", sum(!is.na(fc_matrix)), "\n")
  
  # Create color palette
  color_palette <- colorRampPalette(
    rev(brewer.pal(11, "RdBu"))
  )(100)
  
  # Create heatmap
  if (!is.null(output_file)) {
    pdf(output_file, width = 10, height = 6)
  }
  
  pheatmap(
    fc_matrix,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = color_palette,
    breaks = seq(-3, 3, length.out = 101),
    main = "Pan-Cancer Histone Variant Expression\n(log2 Fold Change: Tumor vs Normal)",
    fontsize = 10,
    na_col = "grey90",
    border_color = "grey60",
    cellwidth = 25,
    cellheight = 25,
    display_numbers = TRUE,
    number_format = "%.2f",
    number_color = "black",
    fontsize_number = 8
  )
  
  if (!is.null(output_file)) {
    dev.off()
    cat("Saved pan-cancer heatmap:", output_file, "\n")
  }
  
  return(fc_matrix)
}


#' Create forest plot for survival results
#' 
#' @param survival_results Data frame with Cox regression results
#' @param gene Gene name to plot
#' @param output_file Path to save plot
create_forest_plot <- function(survival_results,
                               gene = NULL,
                               output_file = NULL) {
  
  library(ggplot2)
  library(dplyr)
  
  cat("\n=== Creating Forest Plot ===\n")
  
  # Filter for specific gene if provided
  if (!is.null(gene)) {
    plot_data <- survival_results %>% filter(gene == !!gene)
    title_text <- paste("Pan-Cancer Prognostic Value of", gene)
  } else {
    plot_data <- survival_results
    title_text <- "Pan-Cancer Survival Analysis"
  }
  
  # Filter significant results
  plot_data <- plot_data %>%
    filter(pvalue < 0.05) %>%
    arrange(HR)
  
  if (nrow(plot_data) == 0) {
    cat("No significant results to plot\n")
    return(NULL)
  }
  
  cat("Plotting", nrow(plot_data), "significant associations\n")
  
  # Create plot
  p <- ggplot(plot_data, aes(x = HR, y = reorder(cancer_type, HR))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 1) +
    geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), 
                   height = 0.2, size = 0.8) +
    geom_point(size = 3, color = "#E41A1C") +
    scale_x_continuous(trans = "log10",
                       breaks = c(0.5, 1, 2, 4, 8)) +
    theme_bw(base_size = 12) +
    labs(
      title = title_text,
      subtitle = "Hazard Ratios (95% CI) for Overall Survival (p < 0.05)",
      x = "Hazard Ratio (log scale)",
      y = "Cancer Type"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid.major.y = element_line(color = "grey90")
    ) +
    geom_text(aes(label = sprintf("HR=%.2f\np=%.3f", HR, pvalue)),
              hjust = -0.2, size = 3)
  
  # Save plot
  if (!is.null(output_file)) {
    ggsave(
      filename = output_file,
      plot = p,
      width = 8,
      height = max(4, nrow(plot_data) * 0.4),
      dpi = 300
    )
    cat("Saved forest plot:", output_file, "\n")
  }
  
  return(p)
}


#' Create multi-panel figure combining multiple plots
#' 
#' @param plots List of ggplot objects
#' @param labels Vector of panel labels (A, B, C, etc.)
#' @param ncol Number of columns
#' @param output_file Path to save plot
create_multipanel_figure <- function(plots,
                                     labels = LETTERS[1:length(plots)],
                                     ncol = 2,
                                     output_file = NULL) {
  
  library(cowplot)
  
  cat("\n=== Creating Multi-Panel Figure ===\n")
  cat("Number of panels:", length(plots), "\n")
  
  # Combine plots
  combined <- plot_grid(
    plotlist = plots,
    labels = labels,
    ncol = ncol,
    align = "hv",
    axis = "tb"
  )
  
  # Save
  if (!is.null(output_file)) {
    ggsave(
      filename = output_file,
      plot = combined,
      width = 8 * ncol,
      height = 6 * ceiling(length(plots) / ncol),
      dpi = 300,
      limitsize = FALSE
    )
    cat("Saved multi-panel figure:", output_file, "\n")
  }
  
  return(combined)
}


#' Create gene correlation heatmap
#' 
#' @param normalized_counts Normalized expression matrix
#' @param genes Vector of gene names
#' @param cancer_type Cancer type label
#' @param output_file Path to save plot
create_correlation_heatmap <- function(normalized_counts,
                                       genes,
                                       cancer_type = "",
                                       output_file = NULL) {
  
  library(pheatmap)
  library(RColorBrewer)
  
  cat("\n=== Creating Correlation Heatmap ===\n")
  
  # Filter for genes
  gene_expr <- normalized_counts[rownames(normalized_counts) %in% genes, ]
  
  if (nrow(gene_expr) < 2) {
    cat("WARNING: Need at least 2 genes for correlation\n")
    return(NULL)
  }
  
  # Calculate correlation
  cor_matrix <- cor(t(gene_expr), method = "spearman")
  
  cat("Correlation matrix:\n")
  print(round(cor_matrix, 3))
  
  # Create heatmap
  if (!is.null(output_file)) {
    pdf(output_file, width = 6, height = 6)
  }
  
  pheatmap(
    cor_matrix,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
    breaks = seq(-1, 1, length.out = 101),
    main = paste("Gene Correlation Heatmap:", cancer_type),
    display_numbers = TRUE,
    number_format = "%.2f",
    fontsize = 12,
    cellwidth = 40,
    cellheight = 40
  )
  
  if (!is.null(output_file)) {
    dev.off()
    cat("Saved correlation heatmap:", output_file, "\n")
  }
  
  return(cor_matrix)
}


#' Create summary table plot
#' 
#' @param summary_table Data frame with analysis summary
#' @param output_file Path to save plot
create_summary_table_plot <- function(summary_table,
                                      output_file = NULL) {
  
  library(ggplot2)
  library(gridExtra)
  
  cat("\n=== Creating Summary Table ===\n")
  
  # Format table
  formatted_table <- summary_table
  
  # Create table plot
  table_plot <- tableGrob(
    formatted_table,
    rows = NULL,
    theme = ttheme_default(
      base_size = 10,
      padding = unit(c(4, 4), "mm")
    )
  )
  
  # Save
  if (!is.null(output_file)) {
    ggsave(
      filename = output_file,
      plot = table_plot,
      width = 12,
      height = max(6, nrow(formatted_table) * 0.3),
      dpi = 300
    )
    cat("Saved summary table:", output_file, "\n")
  }
  
  return(table_plot)
}


#' Create comprehensive figure set for a cancer type
#' 
#' @param data SummarizedExperiment object
#' @param de_results Differential expression results
#' @param survival_results Survival analysis results
#' @param genes Vector of gene names
#' @param cancer_type Cancer type label
#' @param output_dir Directory to save figures
create_cancer_figure_set <- function(data,
                                     de_results,
                                     survival_results,
                                     genes,
                                     cancer_type,
                                     output_dir) {
  
  cat("\n")
  cat("========================================\n")
  cat("Creating Figure Set for", cancer_type, "\n")
  cat("========================================\n")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 1. Expression boxplot
  create_expression_boxplot(
    data, 
    genes, 
    cancer_type,
    file.path(output_dir, paste0(cancer_type, "_boxplot.pdf"))
  )
  
  # 2. Volcano plot
  create_volcano_plot(
    de_results,
    cancer_type = cancer_type,
    output_file = file.path(output_dir, paste0(cancer_type, "_volcano.pdf"))
  )
  
  # 3. Correlation heatmap (if normalized counts available)
  # This would require normalized_counts as input
  
  cat("\nFigure set complete for", cancer_type, "\n")
}