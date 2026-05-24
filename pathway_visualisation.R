suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

output_dir <- "results/BRCA/gsea_pathway_analysis"
all_gsea <- readRDS(paste0(output_dir, "/all_gsea_results.rds"))

histone_genes <- c("H2AX", "H2AZ1", "MACROH2A1", "MACROH2A2")

#------TOP PATHWAYS
extract_top_pathways <- function(gsea_result, n = 20, type = "GO") {

  if (nrow(gsea_result@result) == 0) return(NULL)

  gsea_result@result %>%
    arrange(p.adjust) %>%
    head(n) %>%
    mutate(
      GeneCount = sapply(strsplit(core_enrichment, "/"), length),
      neg_log10_padj = -log10(p.adjust),
      Type = type
    ) %>%
    select(ID, Description, NES, p.adjust, neg_log10_padj, GeneCount, Type)
}

all_go_pathways <- list()

for (gene in histone_genes) {
  go_result <- all_gsea[[gene]]$GO@result

  if (nrow(go_result) > 0) {
    top_pathways <- go_result %>%
      filter(p.adjust < 0.05) %>%
      arrange(p.adjust) %>%
      head(50)  # Top 50 per gene

    all_go_pathways[[gene]] <- top_pathways$Description
  }
}

# Find overlaps
pathway_counts <- table(unlist(all_go_pathways))
common_pathways <- names(pathway_counts[pathway_counts >= 2])  # In at least 2 genes

cat("Total unique pathways:", length(pathway_counts), "\n")
cat("Pathways in ≥2 genes:", length(common_pathways), "\n\n")

# If many common pathways, create combined plot
if (length(common_pathways) >= 10) {

  cat("Creating COMBINED dotplot for", length(common_pathways), "common pathways\n\n")

  # Extract data for common pathways
  combined_data <- data.frame()

  for (gene in histone_genes) {
    go_result <- all_gsea[[gene]]$GO@result

    common_data <- go_result %>%
      filter(Description %in% common_pathways) %>%
      arrange(p.adjust) %>%
      head(20) %>%
      mutate(
        Gene = gene,
        GeneCount = as.numeric(gsub("/.*", "", setSize)),
        neg_log10_padj = -log10(p.adjust)
      ) %>%
      select(Gene, Description, NES, p.adjust, neg_log10_padj, GeneCount)

    combined_data <- rbind(combined_data, common_data)
  }

  # Plot: Combined dotplot
  p_combined <- ggplot(combined_data,
                       aes(x = NES, y = reorder(Description, neg_log10_padj),
                           color = Gene, size = neg_log10_padj)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c(
      "H2AX" = "#E74C3C",
      "H2AZ1" = "#3498DB",
      "MACROH2A1" = "#2ECC71",
      "MACROH2A2" = "#F39C12"
    )) +
    scale_size_continuous(range = c(3, 10)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = "Common GSEA Pathways Across 4 Histone Variants",
      subtitle = paste("BRCA -", length(unique(combined_data$Description)), "pathways enriched in ≥2 genes"),
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      color = "Histone Variant",
      size = "-log10(padj)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      axis.text.y = element_text(size = 9),
      legend.position = "right"
    )

  ggsave(paste0(output_dir, "/combined_pathways_dotplot.jpeg"),
         plot = p_combined, width = 14, height = 12)

  cat("✓ Saved: combined_pathways_dotplot.jpeg\n\n")

} else {
  cat("Too few common pathways - creating individual plots only\n\n")
}


#------ TOP 20 PATHWAYS PER GENE
for (gene in histone_genes) {

  cat("Plotting:", gene, "\n")

  # GO-BP Top 20
  go_top20 <- extract_top_pathways(all_gsea[[gene]]$GO, n = 20, type = "GO-BP")

  if (!is.null(go_top20) && nrow(go_top20) > 0) {

    p_go <- ggplot(go_top20,
                   aes(x = GeneCount, y = reorder(Description, neg_log10_padj),
                       color = NES, size = neg_log10_padj)) +
      geom_point(alpha = 0.8) +
      scale_color_gradient2(
        low = "#2E86AB", mid = "white", high = "#A23B72",
        midpoint = 0,
        limits = c(-max(abs(go_top20$NES)), max(abs(go_top20$NES)))
      ) +
      scale_size_continuous(range = c(3, 10)) +
      labs(
        title = paste(gene, "- Top 20 GO Biological Process Pathways"),
        subtitle = "GSEA with full ranked gene list (no threshold)",
        x = "Gene Set Size",
        y = NULL,
        color = "NES",
        size = "-log10(padj)"
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", size = 13),
        axis.text.y = element_text(size = 9)
      )

    ggsave(paste0(output_dir, "/", gene, "_GO_top20.jpeg"),
           plot = p_go, width = 12, height = 10)
  }

  # KEGG Top 20
  kegg_top20 <- extract_top_pathways(all_gsea[[gene]]$KEGG, n = 20, type = "KEGG")

  if (!is.null(kegg_top20) && nrow(kegg_top20) > 0) {

    p_kegg <- ggplot(kegg_top20,
                     aes(x = GeneCount, y = reorder(Description, neg_log10_padj),
                         color = NES, size = neg_log10_padj)) +
      geom_point(alpha = 0.8) +
      scale_color_gradient2(
        low = "#2E86AB", mid = "white", high = "#A23B72",
        midpoint = 0,
        limits = c(-max(abs(kegg_top20$NES)), max(abs(kegg_top20$NES)))
      ) +
      scale_size_continuous(range = c(3, 10)) +
      labs(
        title = paste(gene, "- Top 20 KEGG Pathways"),
        subtitle = "GSEA with full ranked gene list (no threshold)",
        x = "Gene Set Size",
        y = NULL,
        color = "NES",
        size = "-log10(padj)"
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", size = 13),
        axis.text.y = element_text(size = 9)
      )

    ggsave(paste0(output_dir, "/", gene, "_KEGG_top20.jpeg"),
           plot = p_kegg, width = 12, height = 10)
  }

  cat("  ✓", gene, "plots saved\n")
}

