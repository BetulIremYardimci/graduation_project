library(ggplot2)
library(ggpubr)
library(pheatmap)
library(patchwork)

# Expression Comparasion : Tumor vs Normal
# Function 1: Boxplot
plot_boxplot <- function(data_long, title = "Histone Variant Expression") {

  p <- ggplot(data_long, aes(x = gene_name, y = log2(expression + 1),  # gene_id -> gene_name
                             fill = sample_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(alpha = 0.1, size = 0.3,
                position = position_jitterdodge(jitter.width = 0.2)) +
    scale_fill_manual(values = c("Tumor" = "#E74C3C", "Normal" = "#3498DB")) +
    theme_minimal(base_size = 12) +
    labs(title = title,
         x = "Gene",
         y = "log2(Expression + 1)",
         fill = "Sample Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"))

  return(p)
}

# Function 2: Violin Plot
plot_violin <- function(data_long, title = "Expression Distribution") {

  p <- ggplot(data_long, aes(x = gene_name, y = log2(expression + 1),
                             fill = sample_type)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.15, alpha = 0.8,
                 position = position_dodge(0.9)) +
    scale_fill_manual(values = c("Tumor" = "#E74C3C", "Normal" = "#3498DB")) +
    theme_minimal(base_size = 12) +
    labs(title = title,
         x = "Gene",
         y = "log2(Expression + 1)",
         fill = "Sample Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"))

  return(p)
}

# Function 3: Faceted Plot
plot_faceted <- function(data_long, title = "Expression by Gene") {

  p <- ggplot(data_long, aes(x = sample_type, y = log2(expression + 1),
                             fill = sample_type)) +
    geom_boxplot(alpha = 0.8) +
    facet_wrap(~ gene_name, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("Tumor" = "#E74C3C", "Normal" = "#3498DB")) +
    theme_bw(base_size = 12) +
    labs(title = title,
         x = "",
         y = "log2(Expression + 1)",
         fill = "Sample Type") +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          strip.background = element_rect(fill = "lightgray"))

  return(p)
}


# Function 4: Statistical Comparison
plot_with_stats <- function(data_long, title = "Expression with Statistics") {

  p <- ggplot(data_long, aes(x = gene_name, y = log2(expression + 1),
                             fill = sample_type)) +
    geom_boxplot(alpha = 0.8) +
    stat_compare_means(aes(group = sample_type),
                       label = "p.signif",
                       method = "wilcox.test",
                       label.y.npc = 0.95) +
    scale_fill_manual(values = c("Tumor" = "#E74C3C", "Normal" = "#3498DB")) +
    theme_minimal(base_size = 12) +
    labs(title = title,
         x = "Gene",
         y = "log2(Expression + 1)",
         fill = "Sample Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"))

  return(p)
}

# Function 5: Heatmap
plot_heatmap <- function(data_matrix, annotation_df,
                         title = "Expression Heatmap") {

  # Log transform
  data_log <- log2(data_matrix + 1)

  # Plot
  pheatmap(data_log,
           annotation_col = annotation_df,
           show_colnames = FALSE,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = title,
           fontsize = 10,
           annotation_colors = list(
             sample_type = c("Tumor" = "#E74C3C", "Normal" = "#3498DB")
           ))
}


# Function 6: Combined Plot (Publication Quality)
plot_combined <- function(data_long, output_file = NULL) {

  # Create individual plots
  p1 <- plot_boxplot(data_long, "A. Expression Levels")
  p2 <- plot_violin(data_long, "B. Expression Distribution")
  p3 <- plot_with_stats(data_long, "C. Statistical Comparison")
  p4 <- plot_faceted(data_long, "D. Individual Gene Expression")

  # Combine
  combined <- (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = "Histone Variant Expression: Tumor vs Normal",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )

  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, combined, width = 16, height = 12, dpi = 300)
    message("✓ Plot saved: ", output_file)
  }

  return(combined)
}


# Function 7: Summary Statistics
get_summary_stats <- function(data_long) {

  library(dplyr)

  stats <- data_long %>%
    group_by(gene_name, sample_type) %>%
    summarise(
      n = n(),
      mean = mean(expression, na.rm = TRUE),
      median = median(expression, na.rm = TRUE),
      sd = sd(expression, na.rm = TRUE),
      min = min(expression, na.rm = TRUE),
      max = max(expression, na.rm = TRUE),
      .groups = "drop"
    )

  return(stats)
}


# Function 8: Export All Plots
export_all_plots <- function(data_long, data_matrix, annotation_df,
                             output_dir = "figures") {

  # Create directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Individual plots
  ggsave(paste0(output_dir, "/boxplot.png"),
         plot_boxplot(data_long), width = 10, height = 6, dpi = 300)

  ggsave(paste0(output_dir, "/violin.png"),
         plot_violin(data_long), width = 10, height = 6, dpi = 300)

  ggsave(paste0(output_dir, "/faceted.png"),
         plot_faceted(data_long), width = 10, height = 8, dpi = 300)

  ggsave(paste0(output_dir, "/with_stats.png"),
         plot_with_stats(data_long), width = 10, height = 6, dpi = 300)

  # Heatmap (save as PDF)
  pdf(paste0(output_dir, "/heatmap.pdf"), width = 12, height = 8)
  plot_heatmap(data_matrix, annotation_df)
  dev.off()

  # Combined plot
  plot_combined(data_long, paste0(output_dir, "/combined_plot.png"))

  message("✓ All plots exported to: ", output_dir)
}

#####################
#Differential Expression

plot_expression_heatmap <- function(norm_counts, sample_metadata, gene_names) {

  library(pheatmap)

  # Log transform
  norm_log <- log2(norm_counts + 1)

  # Gene names as rownames
  rownames(norm_log) <- gene_names

  # Sample annotation
  annotation_col <- data.frame(
    Type = sample_metadata$sample_type,
    row.names = colnames(norm_log)
  )

  # Plot
  pheatmap(norm_log,
           annotation_col = annotation_col,
           show_colnames = FALSE,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "Normalized Expression Heatmap",
           fontsize_row = 12,
           annotation_colors = list(
             Type = c("Normal" = "#3498DB", "Tumor" = "#E74C3C")
           ))
}

#-------------

plot_volcano_full <- function(de_results,
                              title = "Volcano Plot",
                              padj_cutoff = 0.05,
                              log2fc_cutoff = 1) {

  library(ggplot2)

  res_df <- de_results$results

  ## classify genes
  res_df$regulation <- "NS"
  res_df$regulation[res_df$padj < padj_cutoff &
                      res_df$log2FoldChange >= log2fc_cutoff] <- "Up"
  res_df$regulation[res_df$padj < padj_cutoff &
                      res_df$log2FoldChange <= -log2fc_cutoff] <- "Down"

  res_df$regulation <- factor(
    res_df$regulation,
    levels = c("Up", "Down", "NS")
  )

  ggplot(res_df,
         aes(x = log2FoldChange,
             y = -log10(padj),
             color = regulation)) +

    ## background points
    geom_point(data = subset(res_df, regulation == "NS"),
               size = 1.2, alpha = 0.4) +

    ## significant points
    geom_point(data = subset(res_df, regulation != "NS"),
               size = 1.8, alpha = 0.9) +

    ## thresholds
    geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff),
               linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(padj_cutoff),
               linetype = "dashed", color = "gray50") +

    ## colors
    scale_color_manual(
      values = c(
        "Up" = "#D73027",
        "Down" = "#4575B4",
        "NS" = "gray70"
      )
    ) +

    ## limits (prevents extreme tails ruining scale)
    coord_cartesian(
      xlim = quantile(res_df$log2FoldChange, c(0.01, 0.99), na.rm = TRUE),
      ylim = c(0, quantile(-log10(res_df$padj), 0.99, na.rm = TRUE))
    ) +

    theme_classic(base_size = 14) +
    labs(
      title = title,
      x = "log2 Fold Change (Tumor vs Normal)",
      y = expression(-log[10]("adjusted p-value")),
      color = NULL
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "top",
      axis.line = element_line(color = "black"),
      axis.text = element_text(color = "black")
    )
}


#--------------



