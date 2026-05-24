# ============================================================================
# Strong Correlations (r > 0.5) for Each Gene
# ============================================================================

library(dplyr)

output_dir <- "results/BRCA/h2ax_network"

# Load correlation matrix (az önce hesapladık)
cor_matrix <- read.csv(paste0(output_dir, "/full_correlation_matrix.csv"),
                       row.names = 1, check.names = FALSE)

cat("============================================\n")
cat("STRONG CORRELATIONS (|r| > 0.5)\n")
cat("============================================\n\n")

# -------------------------
# FUNCTION: Extract Strong Correlations
# -------------------------

get_strong_correlations <- function(cor_matrix, threshold = 0.5) {

  all_results <- list()

  for (gene in rownames(cor_matrix)) {

    # Get correlations for this gene
    gene_cors <- cor_matrix[gene, ]

    # Filter: |r| > threshold, exclude self
    strong_cors <- gene_cors[abs(gene_cors) > threshold & names(gene_cors) != gene]

    if (length(strong_cors) > 0) {
      strong_cors <- sort(strong_cors, decreasing = TRUE)

      all_results[[gene]] <- data.frame(
        Target_Gene = gene,
        Correlated_Gene = names(strong_cors),
        Correlation = strong_cors,
        Direction = ifelse(strong_cors > 0, "Positive", "Negative"),
        row.names = NULL
      )
    }
  }

  return(all_results)
}

# -------------------------
# EXTRACT r > 0.5
# -------------------------

strong_cors <- get_strong_correlations(cor_matrix, threshold = 0.5)

# Summary
cat("=== Summary ===\n")
for (gene in names(strong_cors)) {
  n_cors <- nrow(strong_cors[[gene]])
  cat(sprintf("%-10s: %d strong correlations\n", gene, n_cors))
}

# -------------------------
# DETAILED RESULTS
# -------------------------

cat("\n============================================\n")
cat("DETAILED RESULTS (r > 0.5)\n")
cat("============================================\n\n")

for (gene in names(strong_cors)) {

  cat("---", gene, "---\n")
  print(strong_cors[[gene]])
  cat("\n")
}

# -------------------------
# COMBINE ALL RESULTS
# -------------------------

all_strong_cors <- bind_rows(strong_cors)

cat("=== Total Strong Correlations ===\n")
cat("Total pairs:", nrow(all_strong_cors), "\n")
cat("Positive:", sum(all_strong_cors$Direction == "Positive"), "\n")
cat("Negative:", sum(all_strong_cors$Direction == "Negative"), "\n\n")

# Save
write.csv(all_strong_cors,
          paste0(output_dir, "/strong_correlations_r05.csv"),
          row.names = FALSE)

cat("✓ Saved: strong_correlations_r05.csv\n\n")

# -------------------------
# NETWORK EDGE LIST
# -------------------------

# Create edge list (undirected: remove duplicates)
edge_list <- all_strong_cors %>%
  rowwise() %>%
  mutate(
    Gene1 = min(Target_Gene, Correlated_Gene),
    Gene2 = max(Target_Gene, Correlated_Gene)
  ) %>%
  ungroup() %>%
  select(Gene1, Gene2, Correlation, Direction) %>%
  distinct()

cat("=== Network Edge List ===\n")
cat("Unique edges:", nrow(edge_list), "\n\n")

write.csv(edge_list,
          paste0(output_dir, "/network_edges_r05.csv"),
          row.names = FALSE)

cat("✓ Saved: network_edges_r05.csv\n")

# -------------------------
# VISUALIZATION: Network Plot
# -------------------------

library(igraph)

# Create graph
g <- graph_from_data_frame(
  d = edge_list[, c("Gene1", "Gene2")],
  directed = FALSE
)

# Add edge weights
E(g)$weight <- edge_list$Correlation
E(g)$color <- ifelse(edge_list$Direction == "Positive", "#E74C3C", "#3498DB")
E(g)$width <- abs(E(g)$weight) * 5

# Node colors
node_colors <- ifelse(V(g)$name == "H2AX", "#F39C12", "#95A5A6")

pdf(paste0(output_dir, "/h2ax_network_graph.pdf"), width = 10, height = 10)

plot(g,
     vertex.color = node_colors,
     vertex.size = 25,
     vertex.label.color = "black",
     vertex.label.cex = 1.2,
     vertex.label.font = 2,
     edge.color = E(g)$color,
     edge.width = E(g)$width,
     layout = layout_with_fr(g),
     main = "H2AX Co-expression Network\n(|r| > 0.5, BRCA Tumors)")

# Add legend
legend("bottomright",
       legend = c("Positive (r>0.5)", "Negative (r<-0.5)", "H2AX"),
       col = c("#E74C3C", "#3498DB", "#F39C12"),
       pch = c(15, 15, 19),
       pt.cex = 2,
       bty = "n")

dev.off()

cat("\n✓ Saved: h2ax_network_graph.pdf\n")

cat("\n============================================\n")
cat("Analysis complete!\n")
cat("============================================\n")
