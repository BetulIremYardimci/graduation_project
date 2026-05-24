# ============================================================
# TCGA-UCEC Mutation Analysis - Histone Variants
# Phase 2: Somatic Mutation Profiling
# ============================================================

# --- Package Installation (run once) ---
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("maftools", quietly = TRUE)) {
  cat("Installing maftools from Bioconductor...\n")
  BiocManager::install("maftools", update = FALSE, ask = FALSE)
  cat("\u2713 maftools installed successfully.\n")
}

# --- Load Libraries ---
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(maftools)
  library(dplyr)
  library(ggplot2)
  library(yaml)
})

# Load config
config <- yaml::yaml.load_file("config.yaml")
histone_genes <- unlist(config$genes$tcga_gene_names)

cancer <- "UCEC"
project_name <- paste0("TCGA-", cancer)

results_dir <- paste0("results/", cancer, "/mutation")
plot_dir    <- paste0(results_dir, "/plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

cat("============================================\n")
cat("Mutation Analysis:", cancer, "\n")
cat("Target genes:", paste(histone_genes, collapse = ", "), "\n")
cat("============================================\n\n")

# ============================================================
# Step 1: Download MAF Data
# ============================================================

maf_file <- paste0("data/processed/", cancer, "_maf.RData")

if (file.exists(maf_file)) {
  cat("Loading existing MAF data...\n")
  load(maf_file)
} else {
  cat("Step 1: Downloading mutation data from GDC...\n")

  maf_query <- GDCquery(
    project  = project_name,
    data.category = "Simple Nucleotide Variation",
    data.type     = "Masked Somatic Mutation",
    access        = "open"
  )

  GDCdownload(maf_query)
  maf_data <- GDCprepare(maf_query)

  dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
  save(maf_data, file = maf_file)
  cat("\u2713 MAF data saved:", maf_file, "\n")
}

# ============================================================
# Step 2: Create maftools Object
# ============================================================

cat("\nStep 2: Creating maftools object...\n")

# TCGAbiolinks output -> maftools
maf <- read.maf(maf = maf_data, isTCGA = TRUE)

cat("  Total mutations:", nrow(maf@data), "\n")
cat("  Total samples:", length(unique(maf@data$Tumor_Sample_Barcode)), "\n")
cat("  Genes mutated:", length(unique(maf@data$Hugo_Symbol)), "\n\n")

# ============================================================
# Step 3: General Cohort Summary
# ============================================================

cat("Step 3: Cohort-level mutation summary...\n")

# MAF summary dashboard
jpeg(paste0(plot_dir, "/maf_summary.jpeg"),
     width = 14, height = 10, units = "in", res = 300, quality = 95)
plotmafSummary(maf = maf, rmOutlier = TRUE,
               addStat = "median", dashboard = TRUE)
dev.off()
cat("\u2713 Saved: maf_summary.jpeg\n")

# Oncoplot - top 20 mutated genes
jpeg(paste0(plot_dir, "/oncoplot_top20.jpeg"),
     width = 12, height = 8, units = "in", res = 300, quality = 95)
oncoplot(maf = maf, top = 20, fontSize = 0.7)
dev.off()
cat("\u2713 Saved: oncoplot_top20.jpeg\n")

# TMB comparison across TCGA
jpeg(paste0(plot_dir, "/tcga_compare.jpeg"),
     width = 12, height = 8, units = "in", res = 300, quality = 95)
tcgaCompare(maf = maf, cohortName = cancer)
dev.off()
cat("\u2713 Saved: tcga_compare.jpeg\n")

# ============================================================
# Step 4: Histone Variant Mutation Analysis
# ============================================================

cat("\nStep 4: Histone variant-specific mutation analysis...\n")

# Filter mutations for histone genes
histone_muts <- maf@data %>%
  filter(Hugo_Symbol %in% histone_genes)

cat("  Histone variant mutations found:", nrow(histone_muts), "\n")

if (nrow(histone_muts) > 0) {

  cat("  Mutated histone genes:",
      paste(unique(histone_muts$Hugo_Symbol), collapse = ", "), "\n")

  # --- 4a. Mutation summary per gene ---
  histone_mut_summary <- histone_muts %>%
    group_by(Hugo_Symbol) %>%
    summarise(
      n_mutations      = n(),
      n_samples        = n_distinct(Tumor_Sample_Barcode),
      mutation_types   = paste(unique(Variant_Classification), collapse = "; "),
      most_common_type = names(sort(table(Variant_Classification),
                                    decreasing = TRUE))[1],
      .groups = "drop"
    ) %>%
    arrange(desc(n_mutations))

  # Total sample count for frequency calculation
  total_samples <- length(unique(maf@data$Tumor_Sample_Barcode))
  histone_mut_summary$mutation_freq_pct <- round(
    histone_mut_summary$n_samples / total_samples * 100, 2
  )

  write.csv(histone_mut_summary,
            paste0(results_dir, "/histone_mutation_summary.csv"),
            row.names = FALSE)
  cat("\u2713 Saved: histone_mutation_summary.csv\n")
  print(histone_mut_summary)

  # --- 4b. Oncoplot for histone genes ---
  jpeg(paste0(plot_dir, "/oncoplot_histone.jpeg"),
       width = 12, height = 6, units = "in", res = 300, quality = 95)
  tryCatch({
    oncoplot(maf = maf, genes = histone_genes, fontSize = 0.8)
  }, error = function(e) {
    cat("  Note: Oncoplot skipped -", e$message, "\n")
    plot.new()
    text(0.5, 0.5, "Insufficient mutations for oncoplot", cex = 1.5)
  })
  dev.off()
  cat("\u2713 Saved: oncoplot_histone.jpeg\n")

  # --- 4c. Lollipop plots per gene ---
  for (gene in unique(histone_muts$Hugo_Symbol)) {
    gene_muts <- histone_muts %>% filter(Hugo_Symbol == gene)

    if (nrow(gene_muts) >= 2) {
      jpeg(paste0(plot_dir, "/lollipop_", gene, ".jpeg"),
           width = 10, height = 5, units = "in", res = 300, quality = 95)
      tryCatch({
        lollipopPlot(maf = maf, gene = gene,
                     showMutationRate = TRUE,
                     labelPos = "all")
      }, error = function(e) {
        cat("  Lollipop skipped for", gene, "-", e$message, "\n")
        plot.new()
        text(0.5, 0.5, paste("Lollipop not available for", gene), cex = 1.2)
      })
      dev.off()
      cat("\u2713 Saved: lollipop_", gene, ".jpeg\n")
    }
  }

  # --- 4d. Mutation type distribution bar plot ---
  mut_type_df <- histone_muts %>%
    group_by(Hugo_Symbol, Variant_Classification) %>%
    summarise(count = n(), .groups = "drop")

  p_mut_type <- ggplot(mut_type_df,
                       aes(x = reorder(Hugo_Symbol, -count, sum),
                           y = count,
                           fill = Variant_Classification)) +
    geom_bar(stat = "identity", alpha = 0.9, color = "white", linewidth = 0.3) +
    scale_fill_brewer(palette = "Set2") +
    theme_classic(base_size = 13) +
    labs(
      title = paste0("Histone Variant Mutations \u2014 ", cancer),
      subtitle = paste0("Total: ", nrow(histone_muts), " mutations across ",
                        n_distinct(histone_muts$Tumor_Sample_Barcode), " samples"),
      x = "Gene",
      y = "Number of Mutations",
      fill = "Mutation Type"
    ) +
    theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "grey40"),
      axis.text.x   = element_text(face = "italic", size = 11),
      legend.position = "right"
    )

  ggsave(paste0(plot_dir, "/histone_mutation_types.jpeg"),
         plot = p_mut_type, width = 10, height = 6, dpi = 300,
         device = "jpeg", bg = "white")
  cat("\u2713 Saved: histone_mutation_types.jpeg\n")

  # --- 4e. Mutation frequency bar plot ---
  p_freq <- ggplot(histone_mut_summary,
                   aes(x = reorder(Hugo_Symbol, -mutation_freq_pct),
                       y = mutation_freq_pct)) +
    geom_bar(stat = "identity", fill = "#E63946", alpha = 0.85,
             width = 0.6, color = "white") +
    geom_text(aes(label = paste0(mutation_freq_pct, "%")),
              vjust = -0.5, size = 4, fontface = "bold") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_classic(base_size = 13) +
    labs(
      title = paste0("Histone Variant Mutation Frequency \u2014 ", cancer),
      subtitle = paste0("n = ", total_samples, " samples"),
      x = "Gene",
      y = "Mutation Frequency (%)"
    ) +
    theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "grey40"),
      axis.text.x   = element_text(face = "italic", size = 11)
    )

  ggsave(paste0(plot_dir, "/histone_mutation_frequency.jpeg"),
         plot = p_freq, width = 8, height = 6, dpi = 300,
         device = "jpeg", bg = "white")
  cat("\u2713 Saved: histone_mutation_frequency.jpeg\n")

} else {
  cat("  No mutations found in histone variant genes for", cancer, "\n")
  cat("  This is itself a finding - low mutation rate may indicate\n")
  cat("  dysregulation occurs at expression/epigenetic level rather than genomic.\n")
}

# ============================================================
# Step 5: Co-occurrence & Mutual Exclusivity
# ============================================================

cat("\nStep 5: Co-occurrence analysis...\n")

# Histone genes vs top mutated genes
top_genes <- getGeneSummary(maf)$Hugo_Symbol[1:20]
test_genes <- unique(c(histone_genes, top_genes))

tryCatch({
  somaticInteractions_result <- somaticInteractions(
    maf = maf,
    genes = test_genes,
    pvalue = c(0.05, 0.01)
  )

  write.csv(somaticInteractions_result,
            paste0(results_dir, "/somatic_interactions.csv"),
            row.names = FALSE)
  cat("\u2713 Saved: somatic_interactions.csv\n")

  # Filter for histone gene interactions
  histone_interactions <- somaticInteractions_result %>%
    filter(gene1 %in% histone_genes | gene2 %in% histone_genes) %>%
    filter(pValue < 0.05) %>%
    arrange(pValue)

  if (nrow(histone_interactions) > 0) {
    write.csv(histone_interactions,
              paste0(results_dir, "/histone_interactions_significant.csv"),
              row.names = FALSE)
    cat("\u2713 Significant histone interactions found:", nrow(histone_interactions), "\n")
    print(histone_interactions)
  } else {
    cat("  No significant co-occurrence/exclusivity with histone genes.\n")
  }

}, error = function(e) {
  cat("  Co-occurrence analysis skipped:", e$message, "\n")
})

# ============================================================
# Step 6: Detailed Mutation Table
# ============================================================

cat("\nStep 6: Exporting detailed mutation data...\n")

# Full mutation details for histone genes
if (nrow(histone_muts) > 0) {
  histone_muts_export <- histone_muts %>%
    select(Hugo_Symbol, Chromosome, Start_Position, End_Position,
           Variant_Classification, Variant_Type,
           Reference_Allele, Tumor_Seq_Allele2,
           HGVSc, HGVSp_Short,
           Tumor_Sample_Barcode, t_depth, t_ref_count, t_alt_count,
           IMPACT, SIFT, PolyPhen) %>%
    arrange(Hugo_Symbol, Start_Position)

  write.csv(histone_muts_export,
            paste0(results_dir, "/histone_mutations_detailed.csv"),
            row.names = FALSE)
  cat("\u2713 Saved: histone_mutations_detailed.csv\n")
}

# ============================================================
# Step 7: Summary Report
# ============================================================

cat("\n============================================\n")
cat("MUTATION ANALYSIS COMPLETE:", cancer, "\n")
cat("============================================\n")
cat("Total cohort mutations:", nrow(maf@data), "\n")
cat("Total samples:", length(unique(maf@data$Tumor_Sample_Barcode)), "\n")
cat("Histone mutations:", nrow(histone_muts), "\n")

if (nrow(histone_muts) > 0) {
  cat("\nPer-gene breakdown:\n")
  for (i in seq_len(nrow(histone_mut_summary))) {
    cat(sprintf("  %-10s : %3d mutations in %3d samples (%.2f%%)\n",
                histone_mut_summary$Hugo_Symbol[i],
                histone_mut_summary$n_mutations[i],
                histone_mut_summary$n_samples[i],
                histone_mut_summary$mutation_freq_pct[i]))
  }
} else {
  cat("\nNo histone variant mutations detected.\n")
  cat("Interpretation: Dysregulation likely occurs at transcriptional/\n")
  cat("epigenetic level rather than through somatic mutations.\n")
}

cat("\nResults saved to:", normalizePath(results_dir, mustWork = FALSE), "\n")
cat("Plots saved to:", normalizePath(plot_dir, mustWork = FALSE), "\n")

# Clean up
rm(maf_data)
gc()
