library(ggplot2)
library(ggpubr)
library(pheatmap)
library(patchwork)
library(dplyr)

# ============================================================
# Expression Comparison: Tumor vs Normal
# Updated: JPEG output + Advanced boxplot design
# ============================================================

# Function 1: Advanced Boxplot (side-by-side, all genes in one panel)
plot_boxplot <- function(data_long, title = "Histone Variant Expression") {

  p <- ggplot(data_long, aes(x = gene_name, y = log2(expression + 1),
                             fill = sample_type)) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.85,
      width = 0.65,
      linewidth = 0.5,
      color = "grey30",
      position = position_dodge(0.75),
      notch = TRUE,
      notchwidth = 0.6
    ) +
    geom_jitter(
      aes(color = sample_type),
      alpha = 0.15,
      size = 0.5,
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
      show.legend = FALSE
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 23,
      size = 2.5,
      fill = "white",
      color = "grey20",
      position = position_dodge(0.75),
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = c("Tumor" = "#C0392B", "Normal" = "#2980B9")
    ) +
    scale_color_manual(
      values = c("Tumor" = "#922B21", "Normal" = "#1F618D")
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(linewidth = 0.8, color = "grey40"),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 11),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
      plot.margin = margin(10, 15, 10, 10)
    ) +
    labs(
      title = title,
      x = "Gene",
      y = expression(log[2](Expression + 1)),
      fill = "Sample Type"
    )

  return(p)
}

# Function 2: Faceted Boxplot with log2FC & padj annotations
# (matches the ALS-style figure)
plot_faceted <- function(data_long,
                         de_stats = NULL,
                         title = "Expression by Gene",
                         subtitle = NULL,
                         condition_col = "sample_type",
                         condition_levels = NULL,
                         condition_colors = NULL,
                         condition_label = "Condition",
                         y_lab = expression(log[2](Expression + 1))) {

  # --- set up condition factor and colors ---
  if (!is.null(condition_levels)) {
    data_long[[condition_col]] <- factor(data_long[[condition_col]],
                                         levels = condition_levels)
  }

  if (is.null(condition_colors)) {
    lvls <- levels(factor(data_long[[condition_col]]))
    condition_colors <- setNames(c("#2980B9", "#C0392B"), lvls)
  }

  # --- compute DE annotation labels per gene if provided ---
  if (!is.null(de_stats)) {
    # de_stats should be a data.frame with columns: gene_name, log2FoldChange, padj
    de_stats$label <- paste0(
      "log2FC = ", formatC(de_stats$log2FoldChange, format = "f", digits = 2), "\n",
      "padj = ", formatC(de_stats$padj, format = "e", digits = 2)
    )
  }

  # --- base plot ---
  p <- ggplot(data_long, aes(x = .data[[condition_col]],
                             y = log2(expression + 1),
                             fill = .data[[condition_col]])) +
    geom_boxplot(
      alpha = 0.85,
      width = 0.55,
      linewidth = 0.5,
      color = "grey30",
      outlier.size = 1.2,
      outlier.alpha = 0.5
    ) +
    facet_wrap(~ gene_name, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = condition_colors) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90", color = "grey60"),
      strip.text = element_text(face = "bold", size = 12),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey30"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(10, 15, 10, 10)
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "",
      y = y_lab,
      fill = condition_label
    )

  # --- add log2FC / padj annotation ---
  if (!is.null(de_stats)) {
    # compute y position for annotation (top of each facet)
    y_positions <- data_long %>%
      group_by(gene_name) %>%
      summarise(y_pos = max(log2(expression + 1), na.rm = TRUE) * 1.02,
                .groups = "drop")

    annotation_df <- merge(de_stats, y_positions, by = "gene_name")

    # x position: centered between the two groups
    annotation_df$x_pos <- 1.5

    p <- p +
      geom_text(data = annotation_df,
                aes(x = x_pos, y = y_pos, label = label),
                inherit.aes = FALSE,
                size = 3.5, vjust = 1, hjust = 0.5, color = "grey20")
  }

  return(p)
}

# Function 3: Violin Plot
plot_violin <- function(data_long, title = "Expression Distribution") {

  p <- ggplot(data_long, aes(x = gene_name, y = log2(expression + 1),
                             fill = sample_type)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.15, alpha = 0.8,
                 position = position_dodge(0.9)) +
    scale_fill_manual(values = c("Tumor" = "#C0392B", "Normal" = "#2980B9")) +
    theme_minimal(base_size = 12) +
    labs(title = title,
         x = "Gene",
         y = expression(log[2](Expression + 1)),
         fill = "Sample Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"))

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
    scale_fill_manual(values = c("Tumor" = "#C0392B", "Normal" = "#2980B9")) +
    theme_minimal(base_size = 12) +
    labs(title = title,
         x = "Gene",
         y = expression(log[2](Expression + 1)),
         fill = "Sample Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"))

  return(p)
}

# Function 5: Heatmap
plot_heatmap <- function(data_matrix, annotation_df,
                         title = "Expression Heatmap") {

  data_log <- log2(data_matrix + 1)

  if (is.null(rownames(data_log))) {
    rownames(data_log) <- c("H2AFX", "H2AFY", "H2AFY2", "H2AFZ")
  }

  pheatmap(data_log,
           annotation_col = annotation_df,
           show_colnames = FALSE,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = title,
           fontsize = 10,
           annotation_colors = list(
             sample_type = c("Tumor" = "#C0392B", "Normal" = "#2980B9")
           ))
}

# Function 6: Combined Plot (Publication Quality)
plot_combined <- function(data_long, output_file = NULL) {

  p1 <- plot_boxplot(data_long, "A. Expression Levels")
  p2 <- plot_violin(data_long, "B. Expression Distribution")
  p3 <- plot_with_stats(data_long, "C. Statistical Comparison")
  p4 <- plot_faceted(data_long, title = "D. Individual Gene Expression")

  combined <- (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = "Histone Variant Expression: Tumor vs Normal",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )

  if (!is.null(output_file)) {
    ggsave(output_file, combined, width = 16, height = 12, dpi = 300,
           device = "jpeg", bg = "white")
    message("\u2713 Plot saved: ", output_file)
  }

  return(combined)
}

# Function 7: Summary Statistics
get_summary_stats <- function(data_long) {

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

# Function 8: Export All Plots (JPEG format)
export_all_plots <- function(data_long, data_matrix, annotation_df,
                             de_stats = NULL,
                             output_dir = "figures") {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  ggsave(paste0(output_dir, "/boxplot.jpeg"),
         plot_boxplot(data_long), width = 10, height = 6, dpi = 300,
         device = "jpeg", bg = "white")

  ggsave(paste0(output_dir, "/violin.jpeg"),
         plot_violin(data_long), width = 10, height = 6, dpi = 300,
         device = "jpeg", bg = "white")

  ggsave(paste0(output_dir, "/faceted.jpeg"),
         plot_faceted(data_long, de_stats = de_stats),
         width = 4 * length(unique(data_long$gene_name)), height = 5, dpi = 300,
         device = "jpeg", bg = "white")

  ggsave(paste0(output_dir, "/with_stats.jpeg"),
         plot_with_stats(data_long), width = 10, height = 6, dpi = 300,
         device = "jpeg", bg = "white")

  # Heatmap as JPEG
  jpeg(paste0(output_dir, "/heatmap.jpeg"), width = 12, height = 8,
       units = "in", res = 300, quality = 95)
  plot_heatmap(data_matrix, annotation_df)
  dev.off()

  # Combined plot as JPEG
  plot_combined(data_long, paste0(output_dir, "/combined_plot.jpeg"))

  message("\u2713 All plots exported as JPEG to: ", output_dir)
}

# ============================================================
# Differential Expression Plots
# ============================================================

plot_expression_heatmap <- function(norm_counts, sample_metadata, gene_names) {

  library(pheatmap)

  norm_log <- log2(norm_counts + 1)
  rownames(norm_log) <- gene_names

  annotation_col <- data.frame(
    Type = sample_metadata$sample_type,
    row.names = colnames(norm_log)
  )

  pheatmap(norm_log,
           annotation_col = annotation_col,
           show_colnames = FALSE,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "Normalized Expression Heatmap",
           fontsize_row = 12,
           annotation_colors = list(
             Type = c("Normal" = "#2980B9", "Tumor" = "#C0392B")
           ))
}

# Volcano Plot
plot_volcano_full <- function(de_results,
                              title = "Volcano Plot",
                              padj_cutoff = 0.05,
                              log2fc_cutoff = 1) {

  res_df <- de_results$results

  res_df$regulation <- "NS"
  res_df$regulation[res_df$padj < padj_cutoff &
                      res_df$log2FoldChange >= log2fc_cutoff] <- "Up"
  res_df$regulation[res_df$padj < padj_cutoff &
                      res_df$log2FoldChange <= -log2fc_cutoff] <- "Down"

  res_df$regulation <- factor(res_df$regulation,
                              levels = c("Up", "Down", "NS"))

  ggplot(res_df,
         aes(x = log2FoldChange,
             y = -log10(padj),
             color = regulation)) +
    geom_point(data = subset(res_df, regulation == "NS"),
               size = 1.2, alpha = 0.4) +
    geom_point(data = subset(res_df, regulation != "NS"),
               size = 1.8, alpha = 0.9) +
    geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff),
               linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(padj_cutoff),
               linetype = "dashed", color = "gray50") +
    scale_color_manual(
      values = c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "gray70")
    ) +
    coord_cartesian(
      xlim = quantile(res_df$log2FoldChange, c(0.01, 0.99), na.rm = TRUE),
      ylim = c(0, quantile(-log10(res_df$padj), 0.99, na.rm = TRUE))
    ) +
    theme_classic(base_size = 14) +
    labs(
      title = title,
      x = expression(log[2]~"Fold Change (Tumor vs Normal)"),
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

