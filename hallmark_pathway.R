suppressPackageStartupMessages({
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)       # MSigDB Hallmark gene setleri için
  library(GSVA)          # GSVA skor hesaplama için
  library(dplyr)
  library(ggplot2)
  library(matrixStats)   # rowMedians için
})

histone_genes <- c("H2AX", "H2AZ1", "MACROH2A1", "MACROH2A2")
ensembl_ids <- c(
  "ENSG00000188486",  # H2AX
  "ENSG00000164032",  # H2AZ1
  "ENSG00000134986",  # MACROH2A1
  "ENSG00000172264"   # MACROH2A2
)
gene_map <- setNames(histone_genes, ensembl_ids)

cancer     <- "UCEC"
output_dir <- paste0("results/", cancer, "/gsea_pathway_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────
# 1. Veri yükleme ve VST normalizasyonu
# ─────────────────────────────────────────────
load("data/processed/UCEC_se_combat.RData")
se_tumor  <- se[, colData(se)$sample_type == "Tumor"]

dds       <- DESeqDataSet(se_tumor, design = ~ 1)
dds       <- estimateSizeFactors(dds)
vst_obj   <- vst(dds, blind = FALSE)
vst_matrix <- assay(vst_obj)

# Ensemble ID → Entrez yardımcı fonksiyonu
ensembl_to_entrez <- function(ensembl_vec) {
  mapIds(
    org.Hs.eg.db,
    keys      = gsub("\\..*", "", ensembl_vec),
    column    = "ENTREZID",
    keytype   = "ENSEMBL",
    multiVals = "first"
  )
}

# Ensemble ID → Gene Symbol yardımcı fonksiyonu
ensembl_to_symbol <- function(ensembl_vec) {
  mapIds(
    org.Hs.eg.db,
    keys      = gsub("\\..*", "", ensembl_vec),
    column    = "SYMBOL",
    keytype   = "ENSEMBL",
    multiVals = "first"
  )
}

cat("VST matrix:", nrow(vst_matrix), "×", ncol(vst_matrix), "\n\n")

# ─────────────────────────────────────────────
# 2. MSigDB Hallmark gene set'lerini hazırla
#    (GSEA için: list of Entrez ID vectors)
#    (GSVA için: ayrıca Symbol versiyonu)
# ─────────────────────────────────────────────
cat("MSigDB Hallmark gene setleri yükleniyor...\n")

hallmark_entrez <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene) %>%
  dplyr::mutate(entrez_gene = as.character(entrez_gene))

hallmark_list_entrez <- split(hallmark_entrez$entrez_gene,
                              hallmark_entrez$gs_name)

hallmark_symbol <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

hallmark_list_symbol <- split(hallmark_symbol$gene_symbol,
                              hallmark_symbol$gs_name)

cat("  Toplam Hallmark set sayısı:", length(hallmark_list_entrez), "\n\n")

# ─────────────────────────────────────────────
# 3. GSEA analizi (GO-BP + KEGG + Hallmark)
# ─────────────────────────────────────────────
all_gsea_results <- list()

for (i in seq_along(ensembl_ids)) {

  ens_id    <- ensembl_ids[i]
  gene_name <- gene_map[ens_id]

  cat("===================================\n")
  cat("GSEA:", gene_name, "\n")
  cat("===================================\n")

  # Hedef gen ekspresyonu ile tüm genler arasında korelasyon
  target_expr  <- vst_matrix[ens_id, ]

  correlations <- apply(vst_matrix, 1, function(gene_expr) {
    tryCatch(
      cor(target_expr, gene_expr, method = "spearman", use = "complete.obs"),
      error = function(e) NA
    )
  })

  gene_list <- correlations[!is.na(correlations)]
  gene_list <- gene_list[abs(gene_list) < 1]   # kendisiyle korelasyonu çıkar
  gene_list <- sort(gene_list, decreasing = TRUE)

  # Entrez ID'ye çevir, duplikat temizle
  gene_list_entrez <- ensembl_to_entrez(names(gene_list))

  gene_df <- data.frame(
    ensembl     = names(gene_list),
    entrez      = gene_list_entrez,
    correlation = gene_list,
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(entrez)) %>%
    group_by(entrez) %>%
    slice_max(abs(correlation), n = 1, with_ties = FALSE) %>%
    ungroup()

  gene_list_final <- setNames(gene_df$correlation, gene_df$entrez)
  gene_list_final <- sort(gene_list_final, decreasing = TRUE)

  cat("  Korelasyon listesi:", length(gene_list_final), "gen\n")

  # 3a. GO Biological Process
  cat("  GSEA (GO-BP) çalıştırılıyor...\n")
  gsea_go <- gseGO(
    geneList      = gene_list_final,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    verbose       = FALSE
  )
  cat("    GO-BP sonuç:", nrow(gsea_go@result), "pathway\n")

  # 3b. KEGG
  cat("  GSEA (KEGG) çalıştırılıyor...\n")
  gsea_kegg <- gseKEGG(
    geneList      = gene_list_final,
    organism      = "hsa",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    verbose       = FALSE
  )
  cat("    KEGG sonuç:", nrow(gsea_kegg@result), "pathway\n")

  # 3c. MSigDB Hallmark — GSEA
  cat("  GSEA (MSigDB Hallmark) çalıştırılıyor...\n")
  gsea_hallmark <- GSEA(
    geneList      = gene_list_final,
    TERM2GENE     = hallmark_entrez %>%
      dplyr::rename(term = gs_name, gene = entrez_gene),
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    verbose       = FALSE
  )
  cat("    Hallmark sonuç:", nrow(gsea_hallmark@result), "pathway\n\n")

  # CSV olarak kaydet
  if (nrow(gsea_go@result) > 0)
    write.csv(gsea_go@result,
              file.path(output_dir, paste0(gene_name, "_GSEA_GO.csv")),
              row.names = FALSE)

  if (nrow(gsea_kegg@result) > 0)
    write.csv(gsea_kegg@result,
              file.path(output_dir, paste0(gene_name, "_GSEA_KEGG.csv")),
              row.names = FALSE)

  if (nrow(gsea_hallmark@result) > 0)
    write.csv(gsea_hallmark@result,
              file.path(output_dir, paste0(gene_name, "_GSEA_Hallmark.csv")),
              row.names = FALSE)

  all_gsea_results[[gene_name]] <- list(
    GO       = gsea_go,
    KEGG     = gsea_kegg,
    Hallmark = gsea_hallmark,
    gene_list = gene_list_final
  )
}

saveRDS(all_gsea_results,
        file.path(output_dir, "all_gsea_results.rds"))

cat("\n✓ GSEA tamamlandı (GO-BP + KEGG + Hallmark)\n\n")

# ─────────────────────────────────────────────
# 4. GSVA analizi — Hallmark set aktivasyon skoru
#    Her örnek için pathway düzeyinde skor hesapla
#    (RiePath makalesindeki Tablo 2 ile karşılaştırılabilir)
# ─────────────────────────────────────────────
cat("===================================\n")
cat("GSVA: Hallmark pathway aktivasyon skorları\n")
cat("===================================\n")

# VST matrisini gene symbol'e çevir (GSVA symbol list ile uyum için)
row_symbols <- ensembl_to_symbol(rownames(vst_matrix))

# Sadece symbol'ü olan ve duplikat olmayan genleri al
symbol_df <- data.frame(
  ensembl = rownames(vst_matrix),
  symbol  = row_symbols,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(symbol)) %>%
  group_by(symbol) %>%
  slice(1) %>%   # duplikat symbol'lerde ilkini al
  ungroup()

vst_symbol <- vst_matrix[symbol_df$ensembl, ]
rownames(vst_symbol) <- symbol_df$symbol

cat("  GSVA için matris:", nrow(vst_symbol), "gen ×", ncol(vst_symbol), "örnek\n")

# GSVA nesnesini oluştur ve çalıştır
gsva_param  <- gsvaParam(vst_symbol, hallmark_list_symbol)
gsva_scores <- gsva(gsva_param, verbose = FALSE)

cat("  GSVA skor matrisi:", nrow(gsva_scores), "pathway ×", ncol(gsva_scores), "örnek\n\n")

# GSVA skoru CSV olarak kaydet
write.csv(gsva_scores,
          file.path(output_dir, "GSVA_Hallmark_scores.csv"))

saveRDS(gsva_scores,
        file.path(output_dir, "gsva_hallmark_scores.rds"))

# ─────────────────────────────────────────────
# 5. Histone genleriyle GSVA skor korelasyonu
#    Her histone gen → her Hallmark pathway korelasyonu
# ─────────────────────────────────────────────
cat("Histone gen × Hallmark pathway korelasyonları hesaplanıyor...\n")

# Histone genlerini symbol olarak al
histone_symbols <- c("H2AX", "H2AZ1", "MACROH2A1", "MACROH2A2")

cor_results_list <- list()

for (h_sym in histone_symbols) {
  if (!h_sym %in% rownames(vst_symbol)) {
    cat("  UYARI:", h_sym, "VST matrisinde bulunamadı, atlanıyor.\n")
    next
  }

  h_expr <- vst_symbol[h_sym, ]

  pathway_cors <- apply(gsva_scores, 1, function(pw_scores) {
    tryCatch(
      cor(h_expr, pw_scores, method = "spearman", use = "complete.obs"),
      error = function(e) NA
    )
  })

  cor_df <- data.frame(
    Pathway     = names(pathway_cors),
    Correlation = pathway_cors,
    Gene        = h_sym,
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(abs(Correlation)))

  cor_results_list[[h_sym]] <- cor_df

  write.csv(cor_df,
            file.path(output_dir,
                      paste0(h_sym, "_GSVA_Hallmark_correlation.csv")),
            row.names = FALSE)

  cat("  ", h_sym, "→ en güçlü korelasyon:",
      cor_df$Pathway[1], "(r =", round(cor_df$Correlation[1], 3), ")\n")
}

# ─────────────────────────────────────────────
# 6. Görselleştirme — Her varyant için dot plot + bar plot
# ─────────────────────────────────────────────
# Kaynak:
#   • Dot plot  → GSEA sonuçları (NES, p.adjust, gene set size)
#   • Bar plot  → GSVA skorları (median ± IQR, High vs Low ekspresyon grubu)
# ─────────────────────────────────────────────
cat("\nDot plot ve bar plot görselleri oluşturuluyor...\n")

# Yardımcı: HALLMARK_ önekini kaldır, _ → boşluk, Title Case
clean_pathway <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

# Renk paleti: pozitif NES → kırmızı, negatif → mavi (standart GSEA convention)
gsea_colors <- c("Activated" = "#B2182B", "Suppressed" = "#2166AC")

for (h_sym in histone_symbols) {

  if (!h_sym %in% rownames(vst_symbol)) {
    cat("  WARNING:", h_sym, "not found, skipping.\n")
    next
  }
  cat(" ", h_sym, "plots being prepared...\n")

  # ════════════════════════════════════════════
  # PLOT 1 — GSEA Hallmark Dot Plot
  # X: NES  |  dot size: gene set size  |  color: direction  |  y: pathway
  # ════════════════════════════════════════════
  hallmark_obj <- all_gsea_results[[h_sym]]$Hallmark
  gsea_res <- if (!is.null(hallmark_obj) && isS4(hallmark_obj)) {
    hallmark_obj@result
  } else {
    data.frame()
  }

  if (nrow(gsea_res) > 0) {

    # Anlamlı sonuçları hazırla
    dot_df <- gsea_res %>%
      filter(p.adjust < 0.05) %>%
      mutate(
        Direction = ifelse(NES > 0, "Activated", "Suppressed"),
        Pathway   = clean_pathway(ID),
        # Pathway ismini en fazla 40 karakterde kırp
        Pathway   = ifelse(nchar(Pathway) > 40,
                           paste0(substr(Pathway, 1, 38), ".."),
                           Pathway),
        Pathway   = factor(Pathway, levels = Pathway[order(NES)])
      ) %>%
      arrange(NES)

    if (nrow(dot_df) > 0) {
      p_dot <- ggplot(dot_df,
                      aes(x = NES, y = Pathway,
                          size = setSize, color = Direction)) +
        geom_point(alpha = 0.85) +
        geom_vline(xintercept = 0, linetype = "dashed",
                   color = "grey50", linewidth = 0.4) +
        scale_color_manual(values = gsea_colors,
                           name   = "Direction") +
        scale_size_continuous(range = c(3, 9),
                              name  = "Gene Set Size") +
        scale_x_continuous(expand = expansion(mult = 0.15)) +
        labs(
          title    = paste0(cancer, " — ", h_sym,
                            ": Hallmark GSEA (Dot Plot)"),
          subtitle = paste0("Significant pathways (adjusted p < 0.05), n = ",
                            nrow(dot_df)),
          x        = "Normalized Enrichment Score (NES)",
          y        = NULL
        ) +
        theme_bw(base_size = 12) +
        theme(
          plot.title      = element_text(face = "bold", size = 13),
          plot.subtitle   = element_text(color = "grey40", size = 10),
          axis.text.y     = element_text(size = 10),
          axis.text.x     = element_text(size = 10),
          legend.position = "right",
          panel.grid.major.y = element_line(color = "grey92"),
          panel.grid.minor   = element_blank()
        )

      out_dot <- file.path(output_dir,
                           paste0(h_sym, "_GSEA_Hallmark_dotplot.jpeg"))
      ggsave(out_dot, p_dot,
             width = 9, height = max(4, nrow(dot_df) * 0.38 + 2),
             dpi = 220, units = "in")
      cat("    →", basename(out_dot), "saved.\n")
    } else {
      cat("    No significant Hallmark GSEA results for", h_sym, "\n")
    }
  }

  # ════════════════════════════════════════════
  # PLOT 2 — GSVA Hallmark Bar Plot
  # Median GSVA score per pathway, split by High vs Low expression group
  # Top 20 pathways by absolute median difference (High - Low)
  # ════════════════════════════════════════════

  # High / Low grupları (üst ve alt tertil)
  h_expr   <- vst_symbol[h_sym, ]
  q_lo     <- quantile(h_expr, 0.33)
  q_hi     <- quantile(h_expr, 0.66)
  high_idx <- which(h_expr > q_hi)
  low_idx  <- which(h_expr < q_lo)

  # Her pathway için High ve Low grubun median GSVA skoru
  med_high <- rowMedians(gsva_scores[, high_idx])
  med_low  <- rowMedians(gsva_scores[, low_idx])
  names(med_high) <- rownames(gsva_scores)
  names(med_low)  <- rownames(gsva_scores)

  delta <- med_high - med_low   # fark skoru

  # Top 20 pathway (en büyük mutlak fark)
  top20 <- names(sort(abs(delta), decreasing = TRUE)[1:20])

  bar_df <- data.frame(
    Pathway   = rep(clean_pathway(top20), 2),
    Group     = rep(c(paste0(h_sym, " High"), paste0(h_sym, " Low")),
                    each = length(top20)),
    MedianGSVA = c(med_high[top20], med_low[top20]),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      Pathway = factor(Pathway,
                       levels = clean_pathway(top20[order(delta[top20])])),
      Group   = factor(Group,
                       levels = c(paste0(h_sym, " High"),
                                  paste0(h_sym, " Low")))
    )

  p_bar <- ggplot(bar_df,
                  aes(x = MedianGSVA, y = Pathway, fill = Group)) +
    geom_col(position = position_dodge(width = 0.7),
             width = 0.65, alpha = 0.9) +
    geom_vline(xintercept = 0, color = "grey30",
               linewidth = 0.4, linetype = "solid") +
    scale_fill_manual(
      values = c(
        setNames("#B2182B", paste0(h_sym, " High")),
        setNames("#4393C3", paste0(h_sym, " Low"))
      ),
      name = paste0(h_sym, " Expression")
    ) +
    scale_x_continuous(expand = expansion(mult = 0.12)) +
    labs(
      title    = paste0(cancer, " — ", h_sym,
                        ": Hallmark GSVA Scores (Bar Plot)"),
      subtitle = paste0("Median GSVA score — High vs. Low expression groups",
                        " (top 20 pathways by differential score)"),
      x        = "Median GSVA Score",
      y        = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(color = "grey40", size = 10),
      axis.text.y        = element_text(size = 10),
      axis.text.x        = element_text(size = 10),
      legend.position    = "top",
      legend.title       = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey92"),
      panel.grid.minor   = element_blank()
    )

  out_bar <- file.path(output_dir,
                       paste0(h_sym, "_GSVA_Hallmark_barplot.jpeg"))
  ggsave(out_bar, p_bar,
         width = 9, height = 7, dpi = 220, units = "in")
  cat("    →", basename(out_bar), "saved.\n")
}

cat("\n  All plots completed.\n\n")

# ─────────────────────────────────────────────
# 7. Özet tablo
# ─────────────────────────────────────────────
summary_df <- data.frame(
  Gene           = character(),
  GO_pathways    = integer(),
  KEGG_pathways  = integer(),
  Hallmark_GSEA  = integer(),
  stringsAsFactors = FALSE
)

for (gene_name in histone_genes) {
  if (!is.null(all_gsea_results[[gene_name]])) {
    res <- all_gsea_results[[gene_name]]
    summary_df <- rbind(summary_df, data.frame(
      Gene          = gene_name,
      GO_pathways   = nrow(res$GO@result),
      KEGG_pathways = nrow(res$KEGG@result),
      Hallmark_GSEA = nrow(res$Hallmark@result)
    ))
  }
}

cat("═══════════════════════════════════════\n")
cat("ÖZET — GSEA Sonuç Sayıları\n")
cat("═══════════════════════════════════════\n")
print(summary_df)

cat("\n═══════════════════════════════════════\n")
cat("GSVA Hallmark: En Güçlü Histone Korelasyonları\n")
cat("═══════════════════════════════════════\n")
for (h_sym in names(cor_results_list)) {
  top3 <- head(cor_results_list[[h_sym]], 3)
  cat(h_sym, ":\n")
  for (j in 1:nrow(top3)) {
    cat("  ", top3$Pathway[j], "→ r =", round(top3$Correlation[j], 3), "\n")
  }
}

cat("\n✓ Tüm analizler tamamlandı.\n")
cat("Çıktı dizini:", output_dir, "\n")
