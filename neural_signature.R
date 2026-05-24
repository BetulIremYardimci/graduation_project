
library(msigdbr)
library(dplyr)

# KEGG Legacy'den nöral setler
neural_kegg <- msigdbr(species = "Homo sapiens",
                       category = "C2",
                       subcategory = "CP:KEGG_LEGACY") %>%
  filter(grepl("NEURO|SYNAP|GLUTAMATE|GABA|SEROTONIN|DOPAMINE", gs_name))

# Reactome'dan nöral setler
neural_reactome <- msigdbr(species = "Homo sapiens",
                           category = "C2",
                           subcategory = "CP:REACTOME") %>%
  filter(grepl("NEURON|SYNAP|TRANSMIT|NEURO|GABA|GLUTAMATE", gs_name))

# GO Biological Process'ten nöral setler
neural_go <- msigdbr(species = "Homo sapiens",
                     category = "C5",
                     subcategory = "GO:BP") %>%
  filter(grepl("NEURON|SYNAP|NEURO|AXON|DENDRIT|NEUROTRANSMIT", gs_name))

# Hangi set isimleri bulundu?
cat("=== KEGG ===\n")
neural_kegg$gs_name %>% unique() %>% print()

cat("=== REACTOME ===\n")
neural_reactome$gs_name %>% unique() %>% print()

cat("=== GO:BP ===\n")
neural_go$gs_name %>% unique() %>% length() %>%
  paste("adet GO:BP seti bulundu") %>% cat("\n")
neural_go$gs_name %>% unique() %>% head(20) %>% print()




# === NEURAL SIGNATURE GEN SETİ OLUŞTURMA ===

# Çakışan fonksiyonları dplyr'a bağla
select    <- dplyr::select
filter    <- dplyr::filter
rename    <- dplyr::rename
mutate    <- dplyr::mutate
summarise <- dplyr::summarise

# Katman 1: KEGG - Genel nöroaktif ligand-reseptör (geniş kapsam)
sets_kegg <- c(
  "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION",
  "KEGG_NEUROTROPHIN_SIGNALING_PATHWAY"
)

# Katman 2: REACTOME - Sinaptik/nörotransmitter (mekanistik)
sets_reactome <- c(
  "REACTOME_NEURONAL_SYSTEM",                              # Ana üst set
  "REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES",       # Sinaptik iletim
  "REACTOME_NEUROTRANSMITTER_RELEASE_CYCLE",              # NT salınımı
  "REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION",
  "REACTOME_PROTEIN_PROTEIN_INTERACTIONS_AT_SYNAPSES",    # Sinaps yapısı
  "REACTOME_NEUREXINS_AND_NEUROLIGINS"                    # Sinaptik adhezyon
)

# Katman 3: GO:BP - Fonksiyonel nöral süreçler
sets_go <- c(
  "GOBP_AXON_DEVELOPMENT",
  "GOBP_AXONAL_TRANSPORT",
  "GOBP_ANTEROGRADE_AXONAL_TRANSPORT"
)

# Tüm setleri birleştir - Neural Signature gen listesi
neural_genes_combined <- bind_rows(
  msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_LEGACY") %>%
    filter(gs_name %in% sets_kegg),
  msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
    filter(gs_name %in% sets_reactome),
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
    filter(gs_name %in% sets_go)
) %>%
  select(gs_name, gene_symbol, entrez_gene) %>%
  distinct()

# Kaç unique gen var?
cat("Toplam unique neural gen:", neural_genes_combined$gene_symbol %>% unique() %>% length(), "\n")
cat("Hangi setlerden kaç gen:\n")

library(tibble)

# Ya da doğrudan base R ile:
neural_geneset_list <- split(
  neural_genes_combined$gene_symbol,
  neural_genes_combined$gs_name
)
neural_geneset_list <- lapply(neural_geneset_list, unique)

cat("Gen seti sayısı:", length(neural_geneset_list), "\n")
cat("Toplam unique gen:", neural_genes_combined$gene_symbol %>% unique() %>% length(), "\n")

neural_genes_combined %>%
  group_by(gs_name) %>%
  summarise(n_genes = n_distinct(gene_symbol)) %>%
  arrange(desc(n_genes)) %>%
  print()


# VST matrisini yükle
vst_mat <- readRDS("results/BRCA/vst_counts.rds")

# Yapısını kontrol et
cat("Boyut:", dim(vst_mat), "\n")  # gen x örnek olmalı
cat("Sınıf:", class(vst_mat), "\n")
cat("İlk satır/sütun isimleri:\n")
rownames(vst_mat)[1:5]
colnames(vst_mat)[1:5]


library(org.Hs.eg.db)
library(AnnotationDbi)

# ENSEMBL → SYMBOL dönüşümü
gene_map <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = rownames(vst_mat),
                                  columns = c("SYMBOL"),
                                  keytype = "ENSEMBL")

cat("Eşleşen gen sayısı:", sum(!is.na(gene_map$SYMBOL)), "\n")

# Duplicate ENSEMBL'leri kaldır, NA'ları at
gene_map_clean <- gene_map %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

# VST matrisini filtrele ve satır isimlerini değiştir
vst_symbol <- vst_mat[gene_map_clean$ENSEMBL, ]
rownames(vst_symbol) <- gene_map_clean$SYMBOL

# Duplicate SYMBOL varsa ortalama al
vst_symbol <- vst_symbol[!duplicated(rownames(vst_symbol)), ]

cat("Final matris boyutu:", dim(vst_symbol), "\n")
cat("İlk gen isimleri:", rownames(vst_symbol)[1:5], "\n")


library(GSVA)

# Neural gen setlerinde kaç gen matrisimizde var?
neural_genes_in_matrix <- lapply(neural_geneset_list, function(genes) {
  intersect(genes, rownames(vst_symbol))
})

cat("Her setteki örtüşen gen sayısı:\n")
sapply(neural_genes_in_matrix, length) %>% sort(decreasing = TRUE) %>% print()


# GSVA parametreleri
gsva_param <- gsvaParam(
  exprData = vst_symbol,          # gene x sample matris
  geneSets = neural_geneset_list,
  kcdf = "Gaussian"            # VST verisi için Gaussian
)

# GSVA skorlarını hesapla → gene_set x sample matris
neural_scores <- gsva(gsva_param, verbose = TRUE)

cat("neural_scores boyutu:", dim(neural_scores), "\n")
# Beklenen: 11 set × N örnek


# === MASTER NEURAL SCORE ===
# 11 setin ortalaması → her hasta için tek skor
master_neural_score <- colMeans(neural_scores)

# === H2AFY2 EKSPRESYONu ===
h2ax_expr <- vst_symbol["H2AX", ]

# İkisini data frame'e al
df_corr <- data.frame(
  sample_id      = names(master_neural_score),
  neural_score   = master_neural_score,
  h2ax_expr    = as.numeric(h2ax_expr)
)

cat("Veri hazır:\n")
cat("Satır sayısı:", nrow(df_corr), "\n")
head(df_corr)

# === SPEARMAN KORELASYON ===
cor_test <- cor.test(df_corr$h2ax_expr,
                     df_corr$neural_score,
                     method = "spearman")

cat("\n=== SPEARMAN KORELASYON SONUCU ===\n")
cat("rho:", round(cor_test$estimate, 3), "\n")
cat("p-value:", cor_test$p.value, "\n")



library(ggplot2)
library(ggpubr)

# Tüm histonlar için df hazırla
plot_df <- data.frame(
  neural_score = master_neural_score
)

for (gene in histone_genes) {
  plot_df[[gene]] <- as.numeric(vst_symbol[gene, ])
}

# Her histone için scatter plot
plot_list <- list()

for (gene in histone_genes) {
  rho_val <- cor_results$rho[cor_results$gene == gene]
  pval    <- cor_results$pvalue[cor_results$gene == gene]
  sig_val <- cor_results$sig[cor_results$gene == gene]

  p <- ggplot(plot_df, aes_string(x = gene, y = "neural_score")) +
    geom_point(alpha = 0.3, size = 1.2, color = ifelse(sig_val == "***", "#C0392B", "#7F8C8D")) +
    geom_smooth(method = "lm", se = TRUE, color = "navy", linewidth = 0.8) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("rho = ", rho_val, "\np = ", formatC(pval, format = "e", digits = 2)),
             hjust = 1.1, vjust = 1.5, size = 3.5, fontface = "italic") +
    labs(
      title = gene,
      x = paste(gene, "Expression (VST)"),
      y = "Master Neural Score"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  plot_list[[gene]] <- p
}

# Grid olarak kaydet
library(cowplot)
final_plot <- plot_grid(plotlist = plot_list, ncol = 3, nrow = 2)

ggsave("results/BRCA/neural_score_histone_correlation.jpeg",
       final_plot, width = 15, height = 10, dpi = 300)

cat("Plot kaydedildi!\n")



# Master score yerine tek tek setlere bak
# Belki MACROH2A2 spesifik bir nöral süreçle ilişkilidir

for (i in 1:nrow(neural_scores)) {
  set_name <- rownames(neural_scores)[i]
  ct <- cor.test(as.numeric(vst_symbol["MACROH2A2", ]),
                 neural_scores[i, ],
                 method = "spearman")
  cat(sprintf("%-70s rho=%6.3f p=%.4f\n", set_name, ct$estimate, ct$p.value))
}

