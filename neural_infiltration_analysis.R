
# ============================================================
# BRCA — Histone Variant × Neural Gene Sets (ssGSEA + Correlation)
# ============================================================
library(GSVA)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

dir.create("results/BRCA/neural", recursive = TRUE, showWarnings = FALSE)

# ---- 1. Expression matrisini yükle --------------------------------
# Mevcut BRCA expression RDS dosyasının yolunu güncelle
brca_rds <- file.path("data", "processed", "BRCA_se_combat.RData")
expr_rds <- "results/BRCA/expression_matrix.rds"   # <-- kendi yoluna göre düzenle
cat("Expression matrisi yükleniyor:", expr_rds, "\n")
expr_mat <- readRDS(expr_rds)

# log2(counts+1) normalize edilmiş olduğunu varsayıyoruz
# Eğer ham count ise: expr_mat <- log2(expr_mat + 1)
cat("Matris boyutu:", dim(expr_mat), "\n")

# ---- 2. Neural Gene Sets (MSigDB) ---------------------------------
cat("Neural gen setleri MSigDB'den çekiliyor...\n")

# KEGG — Neuroactive ligand-receptor + Neurotrophin
neural_kegg <- msigdbr(species = "Homo sapiens",
                       category = "C2",
                       subcategory = "CP:KEGG_LEGACY") %>%
  filter(gs_name %in% c(
    "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION",
    "KEGG_NEUROTROPHIN_SIGNALING_PATHWAY"
  ))

# REACTOME — Neuronal system, synaptic transmission
neural_reactome <- msigdbr(species = "Homo sapiens",
                           category = "C2",
                           subcategory = "CP:REACTOME") %>%
  filter(gs_name %in% c(
    "REACTOME_NEURONAL_SYSTEM",
    "REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES",
    "REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION",
    "REACTOME_NEUROTRANSMITTER_RELEASE_CYCLE"
  ))

# GO:BP — Axon guidance, synapse assembly, action potential
neural_go <- msigdbr(species = "Homo sapiens",
                     category = "C5",
                     subcategory = "GO:BP") %>%
  filter(gs_name %in% c(
    "GOBP_AXON_GUIDANCE",
    "GOBP_SYNAPSE_ASSEMBLY",
    "GOBP_ACTION_POTENTIAL",
    "GOBP_REGULATION_OF_SYNAPTIC_TRANSMISSION",
    "GOBP_GLUTAMATE_RECEPTOR_SIGNALING_PATHWAY"
  ))

# Tüm setleri birleştir ve named list formatına çevir
neural_combined <- bind_rows(neural_kegg, neural_reactome, neural_go)

gene_sets <- split(neural_combined$gene_symbol, neural_combined$gs_name)
cat("Toplam", length(gene_sets), "neural gen seti hazır\n")
cat("Set isimleri:\n")
print(names(gene_sets))

# ---- 3. ssGSEA Skorları -------------------------------------------
cat("\nssGSEA hesaplanıyor (biraz zaman alabilir)...\n")

# GSVA 1.x için
gsva_param <- gsvaParam(expr_mat, gene_sets, method = "ssgsea",
                        kcdf = "Gaussian")
ssgsea_scores <- gsva(gsva_param, verbose = TRUE)

# ssGSEA sonuç matrisi: gen_seti × hasta
cat("ssGSEA tamamlandı. Boyut:", dim(ssgsea_scores), "\n")
saveRDS(ssgsea_scores, "results/BRCA/neural/ssgsea_neural_scores.rds")

# ---- 4. Histone Variant Ekspresyonunu Çıkar ----------------------
histone_genes <- c("H2AFZ",    # H2AZ1
                   "MACROH2A1", # H2AFY
                   "MACROH2A2", # H2AFY2
                   "H2AFX")    # H2AX

# Genlerin matristeki isimlerini kontrol et
available <- intersect(histone_genes, rownames(expr_mat))
cat("\nMevcut histone genleri:", paste(available, collapse = ", "), "\n")

histone_expr <- vst_mat[available, , drop = FALSE]

# ---- 5. Korelasyon Matrisi ----------------------------------------
cat("\nKorelasyon matrisi hesaplanıyor...\n")

# Ortak hastalar
histone_expr <- vst_mat[histone_ids$ENSEMBL, ]
rownames(histone_expr) <- histone_ids$SYMBOL

common_samples <- intersect(colnames(histone_expr), colnames(ssgsea_scores))
cat("Ortak hasta sayısı:", length(common_samples), "\n")

histone_sub  <- histone_expr[, common_samples]
ssgsea_sub   <- ssgsea_scores[, common_samples]

# Spearman korelasyon: her histone × her neural set
cor_mat <- matrix(NA,
                  nrow = nrow(histone_sub),
                  ncol = nrow(ssgsea_sub),
                  dimnames = list(rownames(histone_sub),
                                  rownames(ssgsea_sub)))

pval_mat <- cor_mat

for (h in rownames(histone_sub)) {
  for (s in rownames(ssgsea_sub)) {
    test <- cor.test(as.numeric(histone_sub[h, ]),
                     as.numeric(ssgsea_sub[s, ]),
                     method = "spearman", exact = FALSE)
    cor_mat[h, s]  <- test$estimate
    pval_mat[h, s] <- test$p.value
  }
}

saveRDS(cor_mat,  "results/BRCA/neural/correlation_matrix.rds")
saveRDS(pval_mat, "results/BRCA/neural/pvalue_matrix.rds")

write.csv(cor_mat,  "results/BRCA/neural/correlation_matrix.csv")
write.csv(pval_mat, "results/BRCA/neural/pvalue_matrix.csv")
cat("Korelasyon matrisleri kaydedildi.\n")

# ---- 6. Görselleştirme: Heatmap -----------------------------------
cat("\nHeatmap çiziliyor...\n")

# Anlamlılık yıldızları
sig_mat <- ifelse(pval_mat < 0.001, "***",
                  ifelse(pval_mat < 0.01,  "**",
                         ifelse(pval_mat < 0.05,  "*", "")))

# Set isimlerini kısalt (heatmap için)
short_names <- gsub("KEGG_|REACTOME_|GOBP_", "", rownames(cor_mat %>% t() %>% as.data.frame()))
# ya da:
col_labels <- gsub("KEGG_|REACTOME_|GOBP_", "", colnames(cor_mat))
col_labels <- gsub("_", " ", col_labels)
col_labels <- stringr::str_wrap(col_labels, 30)

jpeg("results/BRCA/neural/histone_neural_correlation_heatmap.jpeg",
     width = 2800, height = 1400, res = 200, quality = 95)

pheatmap(cor_mat,
         display_numbers = sig_mat,
         number_color    = "black",
         fontsize_number = 14,
         color           = colorRampPalette(
           rev(brewer.pal(11, "RdBu")))(100),
         breaks          = seq(-0.5, 0.5, length.out = 101),
         border_color    = "white",
         cellwidth       = 80,
         cellheight      = 50,
         cluster_rows    = FALSE,
         cluster_cols    = TRUE,
         labels_col      = col_labels,
         fontsize_row    = 13,
         fontsize_col    = 10,
         main            = "Histone Variants × Neural Gene Sets\n(Spearman Correlation, TCGA-BRCA)",
         angle_col       = 45)

dev.off()
cat("Heatmap kaydedildi.\n")

# ---- 7. Scatter plots: MACROH2A2 × her neural set ----------------
cat("\nMAROCH2A2 scatter plotları çiziliyor...\n")

dir.create("results/BRCA/neural/scatter", showWarnings = FALSE)

macroh2a2_expr <- as.numeric(histone_sub["MACROH2A2", ])

for (set_name in rownames(ssgsea_sub)) {
  neural_score <- as.numeric(ssgsea_sub[set_name, ])

  df_plot <- data.frame(
    MACROH2A2    = macroh2a2_expr,
    Neural_Score = neural_score
  )

  test <- cor.test(macroh2a2_expr, neural_score,
                   method = "spearman", exact = FALSE)

  short <- gsub("KEGG_|REACTOME_|GOBP_", "", set_name)
  short <- gsub("_", " ", short)

  p <- ggplot(df_plot, aes(x = MACROH2A2, y = Neural_Score)) +
    geom_point(alpha = 0.3, size = 1.2, color = "#D6604D") +
    geom_smooth(method = "lm", color = "black", linewidth = 1) +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("r = %.3f\np = %.2e",
                             test$estimate, test$p.value),
             hjust = -0.1, vjust = 1.3, size = 4.5) +
    labs(title = short,
         x = "MACROH2A2 Expression (log2)",
         y = "ssGSEA Score") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(size = 10, face = "bold"))

  fname <- paste0("results/BRCA/neural/scatter/MACROH2A2_vs_",
                  gsub("[^A-Za-z0-9]", "_", set_name), ".jpeg")
  ggsave(fname, p, width = 5, height = 4, dpi = 200)
}

cat("\n=== Tüm analizler tamamlandı ===\n")
cat("Çıktılar: results/BRCA/neural/\n")

