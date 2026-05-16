
# Analysis metadata
metadata <- list(
  cancer_type = cancer,
  date_analyzed = Sys.time(),
  n_samples = ncol(deseq2_results$dds),
  n_genes_analyzed = nrow(deseq2_results$results),
  n_histone_genes = length(ensembl_ids),
  parameters = list(
    padj_threshold = padj_threshold,
    log2fc_threshold = log2fc_threshold,
    use_sva = use_sva,
    correlation_threshold = 0.4
  ),
  session_info = sessionInfo()
)

saveRDS(metadata, paste0(results_dir, "/analysis_metadata.rds"))

cat("✓ Results saved in:", results_dir, "\n")

# --------------------------------------------------------------------
# Step 8: QC Plots
cat("\n>>> Step 8: Creating QC Plots\n")

create_qc_plots(
  dds = deseq2_results$dds,
  output_dir = paste0(results_dir, "/qc")
)

# --------------------------------------------------------------------
#Visualization

cat("\n Creating Visualizations\n")

# Correlation network plot
if (nrow(correlated_genes) > 0) {
  plot_correlation_network(
    correlation_results,
    top_n = 20,
    output_file = paste0(results_dir, "/correlation_network.pdf")
  )
}

# Volcano plot (full transcriptome)
if (exists("plot_volcano_full")) {
  p <- plot_volcano_full(
    de_results = deseq2_results,
    title = paste0(cancer, " – Differential Expression")
  )

  ggsave(
    paste0(results_dir, "/volcano_plot.pdf"),
    plot = p,
    width = 8,
    height = 7
  )
}


# Clean up
rm(list = setdiff(ls(), c("cancer", "results_dir")))
gc()
