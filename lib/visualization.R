library(ggplot2)
library(ggpubr)
library(pheatmap)
library(patchwork)

# ============================================================
# Expression Comparison: Tumor vs Normal
# Updated: JPEG output + Advanced boxplot design
# ============================================================

# Function 1: Advanced Boxplot
plot_boxplot <- function(data_long, title = "Histone Variant Expression") {

  p <- ggplot(data_long, aes(x = gene_name, y = log2(expression + 1),
                             fill = sample_type)) +
    # Boxplot with refined aesthetics
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.85,
      width = 0.65,
      linewidth = 0.5,
      color = "grey30",
      position = position_dodge(0.75),
      notch = TRUE,           # notch for median CI
      notchwidth = 0.6
    ) +
    # Jittered points with better visibility
    geom_jitter(
      aes(color = sample_type),
      alpha = 0.15,
      size = 0.5,
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
      show.legend = FALSE
    ) +
    # Median diamond marker
    stat_summary(
      fun = median,
      geom = "point",
      shape = 23,            # diamond
      size = 2.5,
      fill = "white",
      color = "grey20",
      position = position_dodge(0.75),
      show.legend = FALSE
    ) +
    # Color palette
    scale_fill_manual(
      values = c("Tumor" = "#C0392B", "Normal" = "#2980B9"),
      labels = c("Tumor" = "Tumor", "Normal" = "Normal")
    ) +
    scale_color_manual(
      values = c("Tumor" = "#922B21", "Normal" = "#1F618D")
    ) +
    # Theme
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

# Function 2: Violin Plot
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

# Function 3: Faceted Plot
plot_faceted <- function(data_long, title = "Expression by Gene") {

  p <- ggplot(data_long, aes(x = sample_type, y = log2(expression + 1),
                             fill = sample_type)) +
    geom_boxplot(alpha = 0.8) +
    facet_wrap(~ gene_name, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("Tumor" = "#C0392B", "Normal" = "#2980B9")) +
    theme_bw(base_size = 12) +
    labs(title = title,
         x = "",
         y = expression(log[2](Expression + 1)),
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
  p4 <- plot_faceted(data_long, "D. Individual Gene Expression")

  combined <- (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = "Histone Variant Expression: Tumor vs Normal",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )

  if (!is.null(output_file)) {
    ggsave(output_file, combined, width = 16, height = 12, dpi = 300)
    message("\u2713 Plot saved: ", output_file)
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

# Function 8: Export All Plots (JPEG format)
export_all_plots <- function(data_long, data_matrix, annotation_df,
                             output_dir = "figures") {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # All plots as JPEG (quality = 95 for publication)
  ggsave(paste0(output_dir, "/boxplot.jpeg"),
         plot_boxplot(data_long), width = 10, height = 6, dpi = 300,
         device = "jpeg", bg = "white")

  ggsave(paste0(output_dir, "/violin.jpeg"),
         plot_violin(data_long), width = 10, height = 6, dpi = 300,
         device = "jpeg", bg = "white")

  ggsave(paste0(output_dir, "/faceted.jpeg"),
         plot_faceted(data_long), width = 10, height = 8, dpi = 300,
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
  ggsave(paste0(output_dir, "/combined_plot.jpeg"),
         plot_combined(data_long), width = 16, height = 12, dpi = 300,
         device = "jpeg", bg = "white")

  message("\u2713 All plots exported as JPEG to: ", output_dir)
}
