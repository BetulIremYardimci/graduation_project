# =============================================================================
# BRCA — Histone Varyant GSEA Pathway Analizi
# Genler: H2AX (H2AFX), H2AZ1 (H2AFZ), MACROH2A1 (H2AFY), MACROH2A2 (H2AFY2)
# Yöntemler: GO-BP + KEGG + MSigDB Hallmark (GSEA) + GSVA
# Görselleştirme: (1) Her varyant için ayrı Hallmark dot plot
#                 (2) 4 varyantı karşılaştıran ortak "bubble" dot plot
#                     (yüklenen görseldeki formatta)
# =============================================================================
#
# LİTERATÜR BAĞLAMLARI (kod içinde ilgili yerlerde referans verilmiştir):
#
#  [L1] H2AX / γH2AX:
#       Moyal et al. (2022) PMC8900000 — yüksek H2AX ekspresyonu BRCA'da
#       kötü prognozla ilişkili; DNA hasarı tamir yolakları aktive.
#
#  [L2] H2AZ1 (H2A.Z):
#       - Svotelis et al. (2010) PMID 20023423 — H2A.Z fazla ekspresyonu
#         meme kanseri hücrelerinde proliferasyonu artırır.
#       - Li et al. (2024) s12964-024-01823 — H2AZ1, RELA-HIF1A-EGFR
#         eksenini aktive ederek kanser progresyonunu artırır.
#       - Hsieh et al. (2021) s41419-021-03895 — H2A.Z, meme ve prostat
#         kanserinde bilinen onkojenik rol.
#
#  [L3] MACROH2A1 (mH2A1):
#       - Broggi et al. (2020) Frontiers Oncol — primer meme kanserinde
#         metastatik olgularda mH2A1 ekspresyonu değişir.
#       - Douet et al. (2022) PMC9016624 — TNBC'de mH2A1.1 upregüle,
#         RNA Pol II duraksatılmış genleri regüle eder; bağlamsal tümör
#         baskılayıcı/onkojenik rol.
#       - Gaspar-Maia et al. (2023) PMC9950461 — mH2A varyantları
#         enhancer aktivitesini baskılayarak onkojenik programları kısar.
#
#  [L4] MACROH2A2 (mH2A2):
#       - Filipescu et al. (2023) Nature Cell Biology — mH2A2 kaybı
#         tümör yükünü artırır; kanser ilişkili fibroblastlarda inflamatuar
#         gen ekspresyonunu kontrol eder.
#       - Gaspar-Maia et al. (2023) PMC9950461 — mH2A2, melanoma, akciğer,
#         mesane ve MEME kanserlerinde düşük eksprese; tümör baskılayıcı.
#       - Spallotta et al. (2021) biorxiv — mH2A2, GBM'de BRD4 okupansını
#         negatif düzenler; self-renewal programlarına karşı epigenetik bariyer.
#
# =============================================================================

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(sva)             # ComBat-seq
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(GSVA)
  library(dplyr)
  library(ggplot2)
  library(matrixStats)
  library(ggrepel)
  library(RColorBrewer)
  library(patchwork)
})

# ─────────────────────────────────────────────────────────────────────────────
# 0. KONFIGÜRASYON
# ─────────────────────────────────────────────────────────────────────────────
cancer      <- "BRCA"
output_dir  <- file.path("results", cancer, "gsea_pathway_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Hedef genler: sembol → Ensembl ID eşleştirmesi
# H2AFX yeni adı H2AX; H2AFZ yeni adı H2AZ1 — her ikisi de kullanılabilir
TARGET_GENES <- list(
  H2AX      = "ENSG00000188486",   # H2AFX / γH2AX
  H2AZ1     = "ENSG00000164032",   # H2AFZ / H2A.Z
  MACROH2A1 = "ENSG00000134986",   # H2AFY  / mH2A1
  MACROH2A2 = "ENSG00000172264"    # H2AFY2 / mH2A2
)

histone_symbols <- names(TARGET_GENES)
ensembl_ids     <- unlist(TARGET_GENES)

# Renkler (görseldeki paleti koru)
VARIANT_COLORS <- c(
  H2AX      = "#F4A582",   # somon/turuncu-pembe
  H2AZ1     = "#92C5DE",   # mavi
  MACROH2A1 = "#4DAF4A",   # yeşil
  MACROH2A2 = "#F5C242"    # sarı/altın
)

cat("=== BRCA Histone Varyant GSEA Analizi ===\n")
cat("Hedef genler:", paste(histone_symbols, collapse = ", "), "\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 1. VERİ İNDİRME ve HAZIRLIK
# ─────────────────────────────────────────────────────────────────────────────
brca_rds <- file.path("data", "processed", "BRCA_se_combat.RData")

if (!file.exists(brca_rds)) {

  cat("TCGA-BRCA verisi indiriliyor...\n")
  dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
  dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

  query <- GDCquery(
    project           = "TCGA-BRCA",
    data.category     = "Transcriptome Profiling",
    data.type         = "Gene Expression Quantification",
    workflow.type     = "STAR - Counts",
    experimental.strategy = "RNA-Seq"
  )
  GDCdownload(query, method = "api", files.per.chunk = 50,
              directory = "data/raw")
  se_brca <- GDCprepare(query, directory = "data/raw")

  # Tümör örneklerini ayır (sample type kodu 01 = Primary Solid Tumor)
  se_tumor <- se_brca[, se_brca$sample_type == "Primary Tumor"]
  se_normal <- se_brca[, se_brca$sample_type == "Solid Tissue Normal"]

  cat("Tümör örnek sayısı:", ncol(se_tumor), "\n")
  cat("Normal doku örnek sayısı:", ncol(se_normal), "\n\n")

  # ComBat-seq batch düzeltmesi (plate / TSS bazlı)
  # TCGA'da batch bilgisi genellikle 'plate' veya TSS (barkodun 6-7. karakterleri)
  batch_ids <- substr(colnames(se_tumor), 6, 7)
  counts_raw <- assay(se_tumor, "unstranded")

  cat("ComBat-seq batch düzeltmesi uygulanıyor...\n")
  counts_corrected <- ComBat_seq(
    counts  = counts_raw,
    batch   = batch_ids,
    group   = NULL
  )
  assay(se_tumor, "unstranded") <- counts_corrected

  se <- se_tumor
  save(se, file = brca_rds)
  cat("✓ Veri kaydedildi:", brca_rds, "\n\n")

} else {
  cat("Kayıtlı veri yükleniyor:", brca_rds, "\n\n")
  load(brca_rds)
}

# ─────────────────────────────────────────────────────────────────────────────
# 2. VST NORMALİZASYON
# ─────────────────────────────────────────────────────────────────────────────
cat("DESeq2 VST normalizasyonu...\n")

dds <- DESeqDataSet(se, design = ~ 1)
dds <- estimateSizeFactors(dds)

# Düşük eksprese genleri filtrele (en az 10 toplam okuma)
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

vst_obj    <- vst(dds, blind = FALSE)
vst_matrix <- assay(vst_obj)

cat("VST matrisi:", nrow(vst_matrix), "gen ×", ncol(vst_matrix), "örnek\n\n")

# Hedef genlerin varlığını kontrol et
present_genes <- intersect(ensembl_ids, rownames(vst_matrix))
missing_genes  <- setdiff(ensembl_ids, rownames(vst_matrix))
if (length(missing_genes) > 0) {
  cat("UYARI — Matrisde bulunamayan genler:",
      paste(names(TARGET_GENES)[ensembl_ids %in% missing_genes], collapse = ", "), "\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# 3. YARDIMCİ FONKSİYONLAR
# ─────────────────────────────────────────────────────────────────────────────

# Ensembl → Entrez
ensembl_to_entrez <- function(ensembl_vec) {
  mapIds(org.Hs.eg.db,
         keys      = gsub("\\..*", "", ensembl_vec),
         column    = "ENTREZID",
         keytype   = "ENSEMBL",
         multiVals = "first")
}

# Ensembl → Symbol
ensembl_to_symbol <- function(ensembl_vec) {
  mapIds(org.Hs.eg.db,
         keys      = gsub("\\..*", "", ensembl_vec),
         column    = "SYMBOL",
         keytype   = "ENSEMBL",
         multiVals = "first")
}

# Pathway adını temizle (HALLMARK_ öneki, alt çizgi → boşluk)
clean_pathway <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

# ─────────────────────────────────────────────────────────────────────────────
# 4. MSigDB GENE SET'LERİ HAZIRLA
# ─────────────────────────────────────────────────────────────────────────────
cat("MSigDB gene set'leri yükleniyor...\n")

# --- Hallmark ---
hallmark_entrez <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, entrez_gene) %>%
  mutate(entrez_gene = as.character(entrez_gene))


hallmark_list_entrez <- split(hallmark_entrez$entrez_gene,
                              hallmark_entrez$gs_name)

hallmark_symbol <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, gene_symbol)

hallmark_list_symbol <- split(hallmark_symbol$gene_symbol,
                              hallmark_symbol$gs_name)

cat("  Hallmark set sayısı:", length(hallmark_list_entrez), "\n")

# --- KEGG (C2) ---
# msigdbr v7.5+ sürümünde "CP:KEGG" → "CP:KEGG_LEGACY" veya "CP:KEGG_MEDICUS"
# olarak ikiye ayrıldı. Dinamik olarak mevcut string'i buluyoruz.
all_collections <- msigdbr_collections()

# msigdbr 2026: sütun adı "gs_subcollection"; eski versiyonlarda "gs_subcat"
# Her iki ismi de dene
sub_col <- if ("gs_subcollection" %in% names(all_collections)) {
  "gs_subcollection"
} else if ("gs_subcat" %in% names(all_collections)) {
  "gs_subcat"
} else {
  stop("msigdbr_collections() çıktısında subcategory sütunu bulunamadı.")
}
cat_col <- if ("gs_collection" %in% names(all_collections)) {
  "gs_collection"
} else "gs_cat"

kegg_sub <- all_collections[[sub_col]][
  grepl("KEGG", all_collections[[sub_col]], ignore.case = TRUE) &
    all_collections[[cat_col]] == "C2"
]

if (length(kegg_sub) == 0) {
  stop("KEGG subcategory bulunamadı. msigdbr_collections() çıktısını kontrol edin.")
}


# Tercih sırası: LEGACY > MEDICUS > ilk eşleşen
if ("CP:KEGG_LEGACY" %in% kegg_sub) {
  kegg_sub_use <- "CP:KEGG_LEGACY"
} else if ("CP:KEGG_MEDICUS" %in% kegg_sub) {
  kegg_sub_use <- "CP:KEGG_MEDICUS"
} else {
  kegg_sub_use <- kegg_sub[1]
}
cat("  Kullanılan KEGG subcategory:", kegg_sub_use, "\n")

kegg_sets <- msigdbr(species = "Homo sapiens",
                     category = "C2", subcategory = kegg_sub_use) %>%
  select(gs_name, entrez_gene) %>%
  mutate(entrez_gene = as.character(entrez_gene))

cat("  KEGG set sayısı:", length(unique(kegg_sets$gs_name)), "\n")

# --- Reactome (C2:CP:REACTOME) ---
# 1839 set — KEGG'den çok daha kapsamlı ve güncel; aktif olarak güncelleniyor.
# Mitokondriyal biyogenez, DNA onarımı, hücre döngüsü checkpoint gibi
# histone varyantlarıyla ilişkili yolaklar için zengin kaynak.
reactome_sets <- msigdbr(species = "Homo sapiens",
                         category = "C2", subcategory = "CP:REACTOME") %>%
  select(gs_name, entrez_gene) %>%
  mutate(entrez_gene = as.character(entrez_gene))

cat("  Reactome set sayısı:", length(unique(reactome_sets$gs_name)), "\n")

# --- C6: Oncogenic Signatures ---
# 189 set — MYC, RAS, E2F, EGFR, HER2, PI3K, KRAS gibi kanser driver
# imzaları. H2AZ1'in RELA-HIF1A-EGFR [L2] ve H2AX'in DNA damage
# ilişkisi için kritik; BRCA'ya özgü onkojenik bağlamı yakalar.
oncogenic_sets <- msigdbr(species = "Homo sapiens", category = "C6") %>%
  select(gs_name, entrez_gene) %>%
  mutate(entrez_gene = as.character(entrez_gene))

cat("  C6 Oncogenic set sayısı:", length(unique(oncogenic_sets$gs_name)), "\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 5. GSEA ANALİZİ — Her varyant için GO-BP + KEGG + Reactome + Hallmark + C6
# ─────────────────────────────────────────────────────────────────────────────
# Strateji: Korelasyon bazlı sıralı liste (Spearman r)
# Referans: UCEC analizinde de kullanılan yaklaşım; tek gen için GSEA
#           için en güçlü ve kabul gören strateji.
# ─────────────────────────────────────────────────────────────────────────────
all_gsea_results <- list()

for (h_name in histone_symbols) {

  ens_id <- TARGET_GENES[[h_name]]

  if (!ens_id %in% rownames(vst_matrix)) {
    cat("ATLANADI:", h_name, "(matrisde yok)\n")
    next
  }

  cat("═══════════════════════════════════\n")
  cat("GSEA:", h_name, "\n")
  cat("═══════════════════════════════════\n")

  # Korelasyon hesapla
  target_expr  <- vst_matrix[ens_id, ]
  correlations <- apply(vst_matrix, 1, function(g) {
    tryCatch(
      cor(target_expr, g, method = "spearman", use = "complete.obs"),
      error = function(e) NA_real_
    )
  })

  # Kendisiyle korelasyonu ve NA'ları çıkar
  gene_list <- correlations[!is.na(correlations)]
  gene_list <- gene_list[abs(gene_list) < 0.9999]
  gene_list <- sort(gene_list, decreasing = TRUE)

  # Entrez'e çevir, duplikat temizle
  entrez_map <- ensembl_to_entrez(names(gene_list))

  gene_df <- data.frame(
    ensembl     = names(gene_list),
    entrez      = entrez_map,
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

  # ── 5a. GO Biological Process ──
  cat("  GSEA (GO-BP)...\n")
  gsea_go <- tryCatch(
    gseGO(
      geneList      = gene_list_final,
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      minGSSize     = 15,
      maxGSSize     = 500,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      verbose       = FALSE,
      seed          = 42
    ),
    error = function(e) {
      cat("    HATA (GO-BP):", conditionMessage(e), "\n")
      NULL
    }
  )
  n_go <- if (!is.null(gsea_go)) nrow(gsea_go@result) else 0
  cat("    GO-BP:", n_go, "pathway\n")

  # ── 5b. KEGG (MSigDB CP:KEGG) ──
  cat("  GSEA (KEGG via MSigDB)...\n")
  gsea_kegg <- tryCatch(
    GSEA(
      geneList      = gene_list_final,
      TERM2GENE     = kegg_sets %>% rename(term = gs_name, gene = entrez_gene),
      minGSSize     = 15,
      maxGSSize     = 500,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      verbose       = FALSE,
      seed          = 42
    ),
    error = function(e) {
      cat("    HATA (KEGG):", conditionMessage(e), "\n")
      NULL
    }
  )
  n_kegg <- if (!is.null(gsea_kegg)) nrow(gsea_kegg@result) else 0
  cat("    KEGG:", n_kegg, "pathway\n")

  # ── 5c. Reactome ──
  cat("  GSEA (Reactome)...\n")
  gsea_reactome <- tryCatch(
    GSEA(
      geneList      = gene_list_final,
      TERM2GENE     = reactome_sets %>% rename(term = gs_name, gene = entrez_gene),
      minGSSize     = 15,
      maxGSSize     = 500,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      verbose       = FALSE,
      seed          = 42
    ),
    error = function(e) {
      cat("    HATA (Reactome):", conditionMessage(e), "\n")
      NULL
    }
  )
  n_reactome <- if (!is.null(gsea_reactome)) nrow(gsea_reactome@result) else 0
  cat("    Reactome:", n_reactome, "pathway\n")

  # ── 5d. Hallmark ──
  cat("  GSEA (Hallmark)...\n")
  gsea_hallmark <- tryCatch(
    GSEA(
      geneList      = gene_list_final,
      TERM2GENE     = hallmark_entrez %>% rename(term = gs_name, gene = entrez_gene),
      minGSSize     = 15,
      maxGSSize     = 500,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      verbose       = FALSE,
      seed          = 42
    ),
    error = function(e) {
      cat("    HATA (Hallmark):", conditionMessage(e), "\n")
      NULL
    }
  )
  n_hallmark <- if (!is.null(gsea_hallmark)) nrow(gsea_hallmark@result) else 0
  cat("    Hallmark:", n_hallmark, "pathway\n")

  # ── 5e. C6 Oncogenic Signatures ──
  # Özellikle H2AZ1 (EGFR/HER2/MYC) ve H2AX (DNA damage / E2F) için kritik
  cat("  GSEA (C6 Oncogenic)...\n")
  gsea_oncogenic <- tryCatch(
    GSEA(
      geneList      = gene_list_final,
      TERM2GENE     = oncogenic_sets %>% rename(term = gs_name, gene = entrez_gene),
      minGSSize     = 10,   # C6 setleri küçük olabilir, eşiği düşür
      maxGSSize     = 500,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      verbose       = FALSE,
      seed          = 42
    ),
    error = function(e) {
      cat("    HATA (C6 Oncogenic):", conditionMessage(e), "\n")
      NULL
    }
  )
  n_oncogenic <- if (!is.null(gsea_oncogenic)) nrow(gsea_oncogenic@result) else 0
  cat("    C6 Oncogenic:", n_oncogenic, "imza\n\n")

  # CSV kaydet
  if (!is.null(gsea_go) && n_go > 0)
    write.csv(gsea_go@result,
              file.path(output_dir, paste0(h_name, "_GSEA_GO-BP.csv")),
              row.names = FALSE)

  if (!is.null(gsea_kegg) && n_kegg > 0)
    write.csv(gsea_kegg@result,
              file.path(output_dir, paste0(h_name, "_GSEA_KEGG.csv")),
              row.names = FALSE)

  if (!is.null(gsea_reactome) && n_reactome > 0)
    write.csv(gsea_reactome@result,
              file.path(output_dir, paste0(h_name, "_GSEA_Reactome.csv")),
              row.names = FALSE)

  if (!is.null(gsea_hallmark) && n_hallmark > 0)
    write.csv(gsea_hallmark@result,
              file.path(output_dir, paste0(h_name, "_GSEA_Hallmark.csv")),
              row.names = FALSE)

  if (!is.null(gsea_oncogenic) && n_oncogenic > 0)
    write.csv(gsea_oncogenic@result,
              file.path(output_dir, paste0(h_name, "_GSEA_C6_Oncogenic.csv")),
              row.names = FALSE)

  all_gsea_results[[h_name]] <- list(
    GO         = gsea_go,
    KEGG       = gsea_kegg,
    Reactome   = gsea_reactome,
    Hallmark   = gsea_hallmark,
    Oncogenic  = gsea_oncogenic,
    gene_list  = gene_list_final
  )
}

saveRDS(all_gsea_results,
        file.path(output_dir, "all_gsea_results_BRCA.rds"))

cat("✓ GSEA tamamlandı.\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 6. GSVA — Hallmark pathway aktivasyon skorları (örnek düzeyinde)
# ─────────────────────────────────────────────────────────────────────────────
cat("GSVA Hallmark analizi başlıyor...\n")

# Symbol matrisine çevir
row_symbols <- ensembl_to_symbol(rownames(vst_matrix))

symbol_df <- data.frame(
  ensembl = rownames(vst_matrix),
  symbol  = row_symbols,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(symbol)) %>%
  group_by(symbol) %>%
  slice(1) %>%
  ungroup()

vst_symbol <- vst_matrix[symbol_df$ensembl, ]
rownames(vst_symbol) <- symbol_df$symbol

cat("  GSVA matrisi:", nrow(vst_symbol), "gen ×", ncol(vst_symbol), "örnek\n")

gsva_param  <- gsvaParam(vst_symbol, hallmark_list_symbol)
gsva_scores <- gsva(gsva_param, verbose = FALSE)

write.csv(gsva_scores,
          file.path(output_dir, "GSVA_Hallmark_scores_BRCA.csv"))
saveRDS(gsva_scores,
        file.path(output_dir, "gsva_hallmark_scores_BRCA.rds"))

cat("  GSVA tamamlandı:", nrow(gsva_scores), "pathway ×",
    ncol(gsva_scores), "örnek\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 7. GÖRSELLEŞTİRME
# ─────────────────────────────────────────────────────────────────────────────

# ── 7A. Her varyant için ayrı Hallmark GSEA dot plot ──────────────────────
# [L1-L4] literatür bulgularına uygun beklenen yönler:
#   H2AX      → pozitif NES: DNA damage repair, cell cycle, E2F targets
#   H2AZ1     → pozitif NES: proliferation, MYC targets, G2M checkpoint [L2]
#   MACROH2A1 → negatif NES (baskılama): EMT, inflammatory (bağlama göre) [L3]
#   MACROH2A2 → negatif NES: oncogenic programs, self-renewal pathways [L4]

gsea_colors <- c(Activated = "#B2182B", Suppressed = "#2166AC")

cat("Bireysel Hallmark dot plotlar oluşturuluyor...\n")

for (h_name in histone_symbols) {

  hallmark_obj <- all_gsea_results[[h_name]]$Hallmark
  if (is.null(hallmark_obj)) { cat("  ATLANADI:", h_name, "\n"); next }

  gsea_res <- hallmark_obj@result

  if (nrow(gsea_res) == 0) {
    cat("  Anlamlı Hallmark sonucu yok:", h_name, "\n")
    next
  }

  dot_df <- gsea_res %>%
    filter(p.adjust < 0.05) %>%
    mutate(
      Direction   = ifelse(NES > 0, "Activated", "Suppressed"),
      log10padj   = -log10(p.adjust),
      Pathway     = clean_pathway(ID),
      Pathway     = ifelse(nchar(Pathway) > 42,
                           paste0(substr(Pathway, 1, 40), ".."), Pathway),
      Pathway     = factor(Pathway, levels = Pathway[order(NES)])
    ) %>%
    arrange(NES)

  if (nrow(dot_df) == 0) next

  p_dot <- ggplot(dot_df,
                  aes(x = NES, y = Pathway,
                      size = setSize, color = Direction)) +
    geom_point(alpha = 0.88) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey55", linewidth = 0.4) +
    scale_color_manual(values = gsea_colors, name = "Direction") +
    scale_size_continuous(
      range = c(3, 10),
      name  = "Gene Set Size",
      guide = guide_legend(override.aes = list(color = "grey40"))
    ) +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    labs(
      title    = paste0("TCGA-BRCA — ", h_name,
                        ": Hallmark GSEA"),
      subtitle = paste0(
        "Spearman korelasyon bazlı sıralı liste | ",
        "Anlamlı pathway'ler (adj.p < 0.05), n = ", nrow(dot_df)
      ),
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(color = "grey40", size = 9),
      axis.text.y        = element_text(size = 10),
      axis.text.x        = element_text(size = 10),
      legend.position    = "right",
      panel.grid.major.y = element_line(color = "grey93"),
      panel.grid.minor   = element_blank()
    )

  plot_height <- max(4, nrow(dot_df) * 0.40 + 2.5)
  out_path    <- file.path(output_dir,
                           paste0(h_name, "_GSEA_Hallmark_dotplot_BRCA.jpeg"))
  ggsave(out_path, p_dot,
         width = 9, height = plot_height,
         dpi = 220, units = "in", limitsize = FALSE)
  cat("  →", basename(out_path), "kaydedildi.\n")
}

# ── 7B. ORTAK "MULTI-VARYANT" BUBBLE DOT PLOT ─────────────────────────────
# Yüklenen görselle aynı format:
#   x  = NES
#   y  = Pathway (>=2 varyantda ortak olanlar)
#   renk = Histone Variant
#   boyut = -log10(p.adjust)
# ─────────────────────────────────────────────────────────────────────────────
cat("\nOrtak pathway bubble dot plot oluşturuluyor...\n")

# Tüm Hallmark sonuçlarını birleştir
combined_df <- bind_rows(
  lapply(histone_symbols, function(h) {
    obj <- all_gsea_results[[h]]$Hallmark
    if (is.null(obj) || nrow(obj@result) == 0) return(NULL)
    obj@result %>%
      filter(p.adjust < 0.05) %>%
      mutate(
        Gene      = h,
        log10padj = -log10(p.adjust),
        Pathway   = clean_pathway(ID)
      ) %>%
      select(Pathway, NES, log10padj, setSize, Gene, p.adjust)
  })
)

if (!is.null(combined_df) && nrow(combined_df) > 0) {

  # Kaç varyantda ortak olduğunu say
  pathway_counts <- combined_df %>%
    group_by(Pathway) %>%
    summarise(n_variants = n_distinct(Gene), .groups = "drop")

  # Eşik: en az 2 varyantda geçen pathway'ler; en fazla TOP_N göster
  # En informatif olanları seç: önce 4 varyantda ortak olanlar,
  # sonra en yüksek ortalama |NES|'e göre sırala
  TOP_N_HALLMARK <- 30

  common_pathways <- pathway_counts %>%
    filter(n_variants >= 2) %>%
    left_join(
      combined_df %>%
        group_by(Pathway) %>%
        summarise(mean_abs_NES = mean(abs(NES)), .groups = "drop"),
      by = "Pathway"
    ) %>%
    arrange(desc(n_variants), desc(mean_abs_NES)) %>%
    slice_head(n = TOP_N_HALLMARK)

  cat("  Gösterilecek pathway sayısı:", nrow(common_pathways),
      "(toplam >=2 varyant:", nrow(pathway_counts %>% filter(n_variants >= 2)), ")\n")

  plot_df <- combined_df %>%
    filter(Pathway %in% common_pathways$Pathway)

  # Pathway sıralaması: ortalama NES'e göre
  pathway_order <- plot_df %>%
    group_by(Pathway) %>%
    summarise(mean_NES = mean(NES), .groups = "drop") %>%
    arrange(mean_NES) %>%
    pull(Pathway)

  plot_df <- plot_df %>%
    mutate(
      Pathway = factor(Pathway, levels = pathway_order),
      Gene    = factor(Gene, levels = histone_symbols)
    )

  n_pathways <- length(pathway_order)
  nes_lim    <- max(abs(plot_df$NES)) * 1.15

  p_common <- ggplot(plot_df,
                     aes(x = NES, y = Pathway,
                         size  = log10padj,
                         color = Gene)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey55", linewidth = 0.45) +
    geom_point(alpha = 0.82, stroke = 0.3) +
    scale_color_manual(values = VARIANT_COLORS, name = "Histone Variant") +
    scale_size_continuous(
      range  = c(2, 8),                          # max nokta boyutu küçültüldü
      name   = expression(-log[10](padj)),
      breaks = c(1, 2, 3, 4)
    ) +
    scale_x_continuous(
      limits = c(-nes_lim, nes_lim),
      expand = expansion(mult = 0.05)
    ) +
    labs(
      title    = "Common GSEA Pathways Across 4 Histone Variants",
      subtitle = paste0("TCGA-BRCA — top ", n_pathways,
                        " pathways enriched in >=2 genes"),
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    guides(
      color = guide_legend(title = "Histone Variant",
                           override.aes = list(size = 4)),
      size  = guide_legend(title = expression(-log[10](padj)))
    ) +
    theme_bw(base_size = 11) +                   # base font küçültüldü
    theme(
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(color = "grey40", size = 9),
      axis.text.y        = element_text(size = 8.5),   # y etiket küçültüldü
      axis.text.x        = element_text(size = 9),
      axis.title.x       = element_text(size = 10),
      legend.position    = "right",
      legend.title       = element_text(face = "bold", size = 9),
      legend.text        = element_text(size = 8.5),
      panel.grid.major.y = element_line(color = "grey93"),
      panel.grid.minor   = element_blank(),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = NA)
    )

  # Sabit genişlik + pathway başına 0.28 in yükseklik (önceki: 0.42)
  plot_h <- n_pathways * 0.28 + 2.8
  out_common <- file.path(output_dir,
                          "BRCA_CommonPathways_4variants_bubble_dotplot.jpeg")
  ggsave(out_common, p_common,
         width = 9, height = plot_h,
         dpi = 240, units = "in", limitsize = FALSE)
  cat("  →", basename(out_common),
      sprintf("(%.1f × %.1f in)\n", 9, plot_h))

  # ── 7C. GO-BP ortak bubble dot plot (aynı mantık) ──
  combined_go_df <- bind_rows(
    lapply(histone_symbols, function(h) {
      obj <- all_gsea_results[[h]]$GO
      if (is.null(obj) || nrow(obj@result) == 0) return(NULL)
      obj@result %>%
        filter(p.adjust < 0.05) %>%
        mutate(
          Gene      = h,
          log10padj = -log10(p.adjust),
          Pathway   = Description
        ) %>%
        select(Pathway, NES, log10padj, setSize, Gene, p.adjust)
    })
  )

  if (!is.null(combined_go_df) && nrow(combined_go_df) > 0) {

    go_counts <- combined_go_df %>%
      group_by(Pathway) %>%
      summarise(n_variants = n_distinct(Gene), .groups = "drop")

    TOP_N_GO <- 25   # GO-BP çok pathway üretir; en informatif 25'i göster

    common_go <- go_counts %>%
      filter(n_variants >= 2) %>%
      left_join(
        combined_go_df %>%
          group_by(Pathway) %>%
          summarise(mean_abs_NES = mean(abs(NES)), .groups = "drop"),
        by = "Pathway"
      ) %>%
      arrange(desc(n_variants), desc(mean_abs_NES)) %>%
      slice_head(n = TOP_N_GO)

    cat("  GO-BP gösterilecek pathway sayısı:", nrow(common_go),
        "(toplam >=2 varyant:", nrow(go_counts %>% filter(n_variants >= 2)), ")\n")

    go_plot_df <- combined_go_df %>%
      filter(Pathway %in% common_go$Pathway)

    go_order <- go_plot_df %>%
      group_by(Pathway) %>%
      summarise(mean_NES = mean(NES), .groups = "drop") %>%
      arrange(mean_NES) %>%
      pull(Pathway)

    go_plot_df <- go_plot_df %>%
      mutate(
        Pathway = factor(Pathway, levels = go_order),
        Gene    = factor(Gene, levels = histone_symbols)
      )

    # Uzun pathway isimlerini kırp
    go_plot_df <- go_plot_df %>%
      mutate(Pathway = factor(
        ifelse(nchar(as.character(Pathway)) > 50,
               paste0(substr(as.character(Pathway), 1, 48), ".."),
               as.character(Pathway)),
        levels = unique(ifelse(nchar(go_order) > 50,
                               paste0(substr(go_order, 1, 48), ".."),
                               go_order))
      ))

    nes_lim_go <- max(abs(go_plot_df$NES)) * 1.15
    n_go_paths  <- length(levels(go_plot_df$Pathway))

    p_go_common <- ggplot(go_plot_df,
                          aes(x = NES, y = Pathway,
                              size  = log10padj,
                              color = Gene)) +
      geom_vline(xintercept = 0, linetype = "dashed",
                 color = "grey55", linewidth = 0.45) +
      geom_point(alpha = 0.82, stroke = 0.3) +
      scale_color_manual(values = VARIANT_COLORS, name = "Histone Variant") +
      scale_size_continuous(
        range  = c(2, 8),                        # max nokta boyutu küçültüldü
        name   = expression(-log[10](padj))
      ) +
      scale_x_continuous(
        limits = c(-nes_lim_go, nes_lim_go),
        expand = expansion(mult = 0.05)
      ) +
      labs(
        title    = "Common GO-BP GSEA Pathways Across 4 Histone Variants",
        subtitle = paste0("TCGA-BRCA — top ", n_go_paths,
                          " pathways enriched in >=2 genes"),
        x = "Normalized Enrichment Score (NES)",
        y = NULL
      ) +
      guides(
        color = guide_legend(title = "Histone Variant",
                             override.aes = list(size = 4)),
        size  = guide_legend(title = expression(-log[10](padj)))
      ) +
      theme_bw(base_size = 11) +                 # base font küçültüldü
      theme(
        plot.title         = element_text(face = "bold", size = 13),
        plot.subtitle      = element_text(color = "grey40", size = 9),
        axis.text.y        = element_text(size = 8.5),  # y etiket küçültüldü
        axis.text.x        = element_text(size = 9),
        axis.title.x       = element_text(size = 10),
        legend.position    = "right",
        legend.title       = element_text(face = "bold", size = 9),
        legend.text        = element_text(size = 8.5),
        panel.grid.major.y = element_line(color = "grey93"),
        panel.grid.minor   = element_blank(),
        plot.background    = element_rect(fill = "white", color = NA),
        panel.background   = element_rect(fill = "white", color = NA)
      )

    # Sabit genişlik + pathway başına 0.28 in yükseklik
    go_plot_h <- n_go_paths * 0.28 + 2.8
    out_go <- file.path(output_dir,
                        "BRCA_CommonGO-BP_4variants_bubble_dotplot.jpeg")
    ggsave(out_go, p_go_common,
           width  = 9,
           height = go_plot_h,
           dpi    = 240, units = "in", limitsize = FALSE)
    cat("  →", basename(out_go),
        sprintf("(%.1f × %.1f in)\n", 9, go_plot_h))
  }

} else {
  cat("  UYARI: Ortak pathway bulunamadı veya GSEA sonuçları boş.\n")
}

# ── 7D. GSVA Bar Plot — Her varyant için High vs Low ───────────────────────
cat("\nGSVA Hallmark bar plotlar oluşturuluyor...\n")

for (h_name in histone_symbols) {

  if (!h_name %in% rownames(vst_symbol)) {
    cat("  ATLANADI (symbol matriste yok):", h_name, "\n"); next
  }

  h_expr <- vst_symbol[h_name, ]
  q_lo   <- quantile(h_expr, 0.33)
  q_hi   <- quantile(h_expr, 0.66)

  high_idx <- which(h_expr > q_hi)
  low_idx  <- which(h_expr < q_lo)

  med_high <- rowMedians(gsva_scores[, high_idx])
  med_low  <- rowMedians(gsva_scores[, low_idx])
  names(med_high) <- rownames(gsva_scores)
  names(med_low)  <- rownames(gsva_scores)

  delta <- med_high - med_low

  # Top 20 pathway by |delta|
  top20 <- names(sort(abs(delta), decreasing = TRUE)[1:min(20, length(delta))])

  bar_df <- data.frame(
    Pathway    = rep(clean_pathway(top20), 2),
    Group      = rep(c(paste0(h_name, " High"), paste0(h_name, " Low")),
                     each = length(top20)),
    MedianGSVA = c(med_high[top20], med_low[top20]),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      Pathway = factor(Pathway,
                       levels = clean_pathway(top20[order(delta[top20])])),
      Group   = factor(Group,
                       levels = c(paste0(h_name, " High"),
                                  paste0(h_name, " Low")))
    )

  p_bar <- ggplot(bar_df,
                  aes(x = MedianGSVA, y = Pathway, fill = Group)) +
    geom_col(position = position_dodge(width = 0.7),
             width = 0.64, alpha = 0.90) +
    geom_vline(xintercept = 0, color = "grey30",
               linewidth = 0.4, linetype = "solid") +
    scale_fill_manual(
      values = c(
        setNames("#B2182B", paste0(h_name, " High")),
        setNames("#4393C3", paste0(h_name, " Low"))
      ),
      name = paste0(h_name, " Expression")
    ) +
    scale_x_continuous(expand = expansion(mult = 0.12)) +
    labs(
      title    = paste0("TCGA-BRCA — ", h_name,
                        ": Hallmark GSVA Scores"),
      subtitle = paste0(
        "Median GSVA score — High vs Low expression groups",
        " (top 20 pathways by differential score)"
      ),
      x = "Median GSVA Score",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(color = "grey40", size = 9),
      axis.text.y        = element_text(size = 10),
      axis.text.x        = element_text(size = 10),
      legend.position    = "top",
      legend.title       = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey92"),
      panel.grid.minor   = element_blank()
    )

  out_bar <- file.path(output_dir,
                       paste0(h_name, "_GSVA_Hallmark_barplot_BRCA.jpeg"))
  ggsave(out_bar, p_bar, width = 9, height = 7, dpi = 220, units = "in")
  cat("  →", basename(out_bar), "kaydedildi.\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# 8. ÖZET TABLO
# ─────────────────────────────────────────────────────────────────────────────
cat("\n═══════════════════════════════════════════\n")
cat("ÖZET — BRCA GSEA Sonuç Sayıları\n")
cat("═══════════════════════════════════════════\n")

summary_df <- bind_rows(lapply(histone_symbols, function(h) {
  res <- all_gsea_results[[h]]
  if (is.null(res)) return(NULL)
  data.frame(
    Gene          = h,
    GO_BP         = if (!is.null(res$GO))        nrow(res$GO@result)        else 0,
    KEGG          = if (!is.null(res$KEGG))       nrow(res$KEGG@result)      else 0,
    Reactome      = if (!is.null(res$Reactome))   nrow(res$Reactome@result)  else 0,
    Hallmark_GSEA = if (!is.null(res$Hallmark))   nrow(res$Hallmark@result)  else 0,
    C6_Oncogenic  = if (!is.null(res$Oncogenic))  nrow(res$Oncogenic@result) else 0
  )
}))

print(summary_df)

write.csv(summary_df,
          file.path(output_dir, "BRCA_GSEA_summary.csv"),
          row.names = FALSE)

